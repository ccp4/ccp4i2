/**
 * Moorhen Scene parse/serialise/validate.
 *
 * Pure functions over the MoorhenScene type. No DOM, no Moorhen runtime,
 * no React — that keeps this trivially testable and reusable in Node
 * (e.g. CLI scene linter) as well as in the renderer.
 *
 * The lifter (Moorhen state → Scene) and the resolver (Scene → Moorhen
 * backupSession) live in separate files because they have very different
 * dependencies. This file is just the format itself.
 */

import * as YAML from "yaml";

import {
  MoorhenScene,
  SCENE_SCHEMA_VERSION,
  SceneDomain,
  SceneFileRef,
  SceneElement,
  SceneLsqMatch,
  SceneMap,
  SceneMapColumns,
  SceneRepresentation,
  SceneCentre,
  SceneColour,
  SceneColourSelection,
  SceneSuperpose,
  isSceneHexColour,
  isSceneNamedColour,
  isSceneRawColour,
} from "../types/moorhen-scene";

// --------------------------------------------------------------------------
// Errors
// --------------------------------------------------------------------------

/** A single validation problem, located by JSON-pointer-ish path. */
export interface SceneValidationError {
  /** Dotted path to the offending field, e.g. "files[1].fileId". */
  path: string;
  message: string;
}

export class SceneParseError extends Error {
  constructor(public readonly errors: SceneValidationError[]) {
    super(
      `Invalid Moorhen scene: ${errors.length} error(s)\n` +
        errors.map((e) => `  - ${e.path}: ${e.message}`).join("\n"),
    );
    this.name = "SceneParseError";
  }
}

// --------------------------------------------------------------------------
// Parse: YAML string → MoorhenScene
// --------------------------------------------------------------------------

/**
 * Parse a YAML scene document into a typed MoorhenScene.
 * Throws SceneParseError if validation fails.
 */
export function parseScene(yamlText: string): MoorhenScene {
  let raw: unknown;
  try {
    raw = YAML.parse(yamlText);
  } catch (e) {
    throw new SceneParseError([
      { path: "", message: `YAML parse error: ${(e as Error).message}` },
    ]);
  }
  const { scene, errors } = validateScene(raw);
  if (errors.length > 0) throw new SceneParseError(errors);
  return scene;
}

/**
 * Validate a parsed JS value against the scene schema.
 * Returns the (possibly partially-constructed) scene plus a list of errors.
 * Callers that need exceptions should use parseScene().
 */
export function validateScene(raw: unknown): {
  scene: MoorhenScene;
  errors: SceneValidationError[];
} {
  const errors: SceneValidationError[] = [];

  if (!isObject(raw)) {
    errors.push({ path: "", message: "scene must be a YAML mapping" });
    return { scene: { scene: "", version: SCENE_SCHEMA_VERSION }, errors };
  }

  const sceneName = strField(raw, "scene", "", errors);
  const version = numField(raw, "version", SCENE_SCHEMA_VERSION, errors);
  if (version !== SCENE_SCHEMA_VERSION) {
    errors.push({
      path: "version",
      message: `unsupported scene version ${version}; this reader handles version ${SCENE_SCHEMA_VERSION}`,
    });
  }

  const out: MoorhenScene = { scene: sceneName, version };

  if ("authoredIn" in raw) {
    const a = (raw as Record<string, unknown>).authoredIn;
    if (isObject(a)) {
      out.authoredIn = {
        projectId: optStr(a, "projectId", "authoredIn.projectId", errors),
        projectName: optStr(a, "projectName", "authoredIn.projectName", errors),
        createdAt: optStr(a, "createdAt", "authoredIn.createdAt", errors),
        createdBy: optStr(a, "createdBy", "authoredIn.createdBy", errors),
        ccp4i2Version: optStr(a, "ccp4i2Version", "authoredIn.ccp4i2Version", errors),
      };
    } else {
      errors.push({ path: "authoredIn", message: "must be a mapping" });
    }
  }

  if ("files" in raw) {
    out.files = validateFiles((raw as Record<string, unknown>).files, errors);
  }

  if ("superpose" in raw) {
    out.superpose = validateSuperpose(
      (raw as Record<string, unknown>).superpose,
      out.files ?? [],
      errors,
    );
  }

  if ("globalDictionaries" in raw) {
    out.globalDictionaries = validateDictionaryNameList(
      (raw as Record<string, unknown>).globalDictionaries,
      out.files ?? [],
      "globalDictionaries",
      errors,
    );
  }

  if ("domains" in raw) {
    out.domains = validateDomains((raw as Record<string, unknown>).domains, errors);
  }

  if ("elements" in raw) {
    out.elements = validateElements(
      (raw as Record<string, unknown>).elements,
      out.files ?? [],
      errors,
    );
  }

  if ("maps" in raw) {
    out.maps = validateMaps(
      (raw as Record<string, unknown>).maps,
      out.files ?? [],
      errors,
    );
  }

  if ("activeMap" in raw) {
    const am = optStr(raw, "activeMap", "activeMap", errors);
    if (am) {
      // Validate cross-reference against the maps block. The maps block
      // is optional, so we accept any string here when it isn't present
      // (the resolver will warn at apply time).
      if (out.maps && !out.maps.some((m) => m.name === am)) {
        errors.push({
          path: "activeMap",
          message: `"${am}" does not name any entry in maps:`,
        });
      } else {
        out.activeMap = am;
      }
    }
  }

  if ("view" in raw) {
    const v = (raw as Record<string, unknown>).view;
    if (isObject(v)) {
      out.view = {
        origin: tupleField<3>(v, "origin", 3, "view.origin", errors),
        centre: validateCentre(v.centre, out.files ?? [], errors),
        quat: tupleField<4>(v, "quat", 4, "view.quat", errors),
        zoom: optNum(v, "zoom", "view.zoom", errors),
        clipStart: optNum(v, "clipStart", "view.clipStart", errors),
        clipEnd: optNum(v, "clipEnd", "view.clipEnd", errors),
        fogStart: optNum(v, "fogStart", "view.fogStart", errors),
        fogEnd: optNum(v, "fogEnd", "view.fogEnd", errors),
        background: optStr(v, "background", "view.background", errors),
      };
    } else {
      errors.push({ path: "view", message: "must be a mapping" });
    }
  }

  if ("resolver" in raw) {
    const r = (raw as Record<string, unknown>).resolver;
    if (isObject(r)) {
      const policy = optStr(r, "onMissingResidues", "resolver.onMissingResidues", errors);
      if (policy && policy !== "clamp-and-log" && policy !== "strict") {
        errors.push({
          path: "resolver.onMissingResidues",
          message: `must be "clamp-and-log" or "strict", got "${policy}"`,
        });
      } else if (policy) {
        out.resolver = { onMissingResidues: policy as "clamp-and-log" | "strict" };
      }
    } else {
      errors.push({ path: "resolver", message: "must be a mapping" });
    }
  }

  return { scene: out, errors };
}

// --------------------------------------------------------------------------
// Validation helpers for sub-objects
// --------------------------------------------------------------------------

function validateFiles(
  raw: unknown,
  errors: SceneValidationError[],
): SceneFileRef[] | undefined {
  if (!Array.isArray(raw)) {
    errors.push({ path: "files", message: "must be a sequence" });
    return undefined;
  }
  const seen = new Set<string>();
  const out: SceneFileRef[] = [];
  raw.forEach((entry, i) => {
    const p = `files[${i}]`;
    if (!isObject(entry)) {
      errors.push({ path: p, message: "must be a mapping" });
      return;
    }
    const name = strField(entry, "name", "", errors, p);
    if (!name) {
      errors.push({ path: `${p}.name`, message: "required" });
    } else if (seen.has(name)) {
      errors.push({ path: `${p}.name`, message: `duplicate file name "${name}"` });
    } else {
      seen.add(name);
    }

    const ref: SceneFileRef = { name };
    if ("kind" in entry) {
      const kind = optStr(entry, "kind", `${p}.kind`, errors);
      if (kind && kind !== "coordinates" && kind !== "dictionary" && kind !== "mtz") {
        errors.push({
          path: `${p}.kind`,
          message: `must be "coordinates", "dictionary", or "mtz", got "${kind}"`,
        });
      } else if (kind === "coordinates" || kind === "dictionary" || kind === "mtz") {
        ref.kind = kind;
      }
    }
    if ("pdb" in entry) ref.pdb = optStr(entry, "pdb", `${p}.pdb`, errors);
    if ("cifText" in entry) ref.cifText = optStr(entry, "cifText", `${p}.cifText`, errors);
    if ("bundle" in entry) ref.bundle = optStr(entry, "bundle", `${p}.bundle`, errors);
    if ("url" in entry) ref.url = optStr(entry, "url", `${p}.url`, errors);
    if ("path" in entry) ref.path = optStr(entry, "path", `${p}.path`, errors);
    if ("projectId" in entry) ref.projectId = optStr(entry, "projectId", `${p}.projectId`, errors);
    if ("projectName" in entry) ref.projectName = optStr(entry, "projectName", `${p}.projectName`, errors);
    if ("fileId" in entry) ref.fileId = optNum(entry, "fileId", `${p}.fileId`, errors);
    if ("job" in entry) ref.job = optNum(entry, "job", `${p}.job`, errors);
    if ("param" in entry) ref.param = optStr(entry, "param", `${p}.param`, errors);

    // pdb id sanity check (4 char or extended 8 char). Keep it permissive
    // but flag obvious nonsense.
    if (ref.pdb && !/^[0-9A-Za-z]{4}$|^pdb_[0-9a-z]{8}$/.test(ref.pdb)) {
      errors.push({
        path: `${p}.pdb`,
        message: `does not look like a PDB ID (got "${ref.pdb}")`,
      });
    }

    // PDB IDs only make sense for coordinate refs.
    if (ref.pdb && ref.kind === "dictionary") {
      errors.push({
        path: `${p}.pdb`,
        message: "pdb: is only valid for coordinate refs, not dictionaries",
      });
    }

    // cifText only makes sense for dictionary refs (it's the inline-dict
    // shape the lifter writes when a loaded dict has no recoverable URL).
    if (ref.cifText && ref.kind !== "dictionary") {
      errors.push({
        path: `${p}.cifText`,
        message: "cifText: is only valid for dictionary refs (kind: dictionary)",
      });
    }

    // At least one resolvable reference must be set.
    const hasPdb = !!ref.pdb;
    const hasCifText = !!ref.cifText;
    const hasBundle = !!ref.bundle;
    const hasUrl = !!ref.url;
    const hasPath = !!ref.path;
    const hasFileId = ref.fileId !== undefined;
    const hasJobParam = ref.job !== undefined && !!ref.param;
    if (!hasPdb && !hasCifText && !hasBundle && !hasUrl && !hasPath && !hasFileId && !hasJobParam) {
      errors.push({
        path: p,
        message: "must set one of: pdb, url, path, bundle, fileId (+projectId), job+param (+projectId), or cifText (for dictionaries)",
      });
    }
    // Project-internal forms need a projectId.
    if ((hasFileId || hasJobParam) && !ref.projectId) {
      errors.push({
        path: `${p}.projectId`,
        message: "required when fileId or job+param is set",
      });
    }
    // job and param are paired.
    if ((ref.job !== undefined) !== !!ref.param) {
      errors.push({
        path: p,
        message: "job and param must be set together",
      });
    }

    out.push(ref);
  });
  return out;
}

const RANGE_RE = /^(-?\d+)-(-?\d+)$/;

function validateSuperpose(
  raw: unknown,
  files: SceneFileRef[],
  errors: SceneValidationError[],
): SceneSuperpose[] | undefined {
  if (!Array.isArray(raw)) {
    errors.push({ path: "superpose", message: "must be a sequence" });
    return undefined;
  }
  const fileNames = new Set(files.map((f) => f.name));
  const out: SceneSuperpose[] = [];

  raw.forEach((entry, i) => {
    const p = `superpose[${i}]`;
    if (!isObject(entry)) {
      errors.push({ path: p, message: "must be a mapping" });
      return;
    }
    const method = strField(entry, "method", "", errors, p);
    const move = strField(entry, "move", "", errors, p);
    const onto = strField(entry, "onto", "", errors, p);
    if (!method) {
      errors.push({ path: `${p}.method`, message: 'required ("ssm" or "lsq")' });
    } else if (method !== "ssm" && method !== "lsq") {
      errors.push({
        path: `${p}.method`,
        message: `must be "ssm" or "lsq", got "${method}"`,
      });
    }
    if (!move) errors.push({ path: `${p}.move`, message: "required" });
    if (!onto) errors.push({ path: `${p}.onto`, message: "required" });

    // Cross-reference the files block — surfaces typos early.
    if (move && files.length > 0 && !fileNames.has(move)) {
      errors.push({
        path: `${p}.move`,
        message: `unknown file "${move}" (not in top-level files block)`,
      });
    }
    if (onto && files.length > 0 && !fileNames.has(onto)) {
      errors.push({
        path: `${p}.onto`,
        message: `unknown file "${onto}" (not in top-level files block)`,
      });
    }
    if (move && onto && move === onto) {
      errors.push({
        path: p,
        message: `cannot superpose a file onto itself ("${move}")`,
      });
    }

    if (method === "ssm") {
      const movChain = strField(entry, "movChain", "", errors, p);
      const refChain = strField(entry, "refChain", "", errors, p);
      if (!movChain) errors.push({ path: `${p}.movChain`, message: "required for ssm" });
      if (!refChain) errors.push({ path: `${p}.refChain`, message: "required for ssm" });
      out.push({ method: "ssm", move, onto, movChain, refChain });
      return;
    }

    if (method === "lsq") {
      const hasMatches = "matches" in entry;
      const chain = optStr(entry, "chain", `${p}.chain`, errors);
      const range = rangeField(entry, "range", errors, p);
      const hasShorthand = !!chain || !!range;

      // Mutually exclusive — having both invites silent precedence bugs.
      if (hasMatches && hasShorthand) {
        errors.push({
          path: p,
          message: 'use either `matches` or the `chain`+`range` shorthand, not both',
        });
      }

      let matches: SceneLsqMatch[] = [];

      if (hasShorthand) {
        // chain + range together; both required.
        if (!chain) {
          errors.push({ path: `${p}.chain`, message: "required when range is set" });
        }
        if (!range) {
          errors.push({ path: `${p}.range`, message: "required when chain is set" });
        } else if (!RANGE_RE.test(range)) {
          errors.push({
            path: `${p}.range`,
            message: `must be "start-end", got "${range}"`,
          });
        } else {
          const m = RANGE_RE.exec(range)!;
          const s = parseInt(m[1], 10);
          const e = parseInt(m[2], 10);
          if (e < s) {
            errors.push({
              path: `${p}.range`,
              message: `end (${e}) must be >= start (${s})`,
            });
          }
        }
        // Even without `matches`, this is a valid LSQ entry — the resolver
        // expands the shorthand at apply-time.
      } else {
        const matchesRaw = (entry as Record<string, unknown>).matches;
        if (!Array.isArray(matchesRaw)) {
          errors.push({
            path: `${p}.matches`,
            message:
              "required for lsq (must be a sequence of range matches) — or use the `chain`+`range` shorthand",
          });
        } else if (matchesRaw.length === 0) {
          errors.push({
            path: `${p}.matches`,
            message: "must contain at least one match entry",
          });
        } else {
          matches = matchesRaw
            .map((m, j) => validateLsqMatch(m, `${p}.matches[${j}]`, errors))
            .filter((m): m is SceneLsqMatch => m !== null);
        }
      }

      const matchType = optStr(entry, "matchType", `${p}.matchType`, errors);
      if (matchType && matchType !== "all" && matchType !== "main" && matchType !== "ca") {
        errors.push({
          path: `${p}.matchType`,
          message: `must be "all", "main" or "ca", got "${matchType}"`,
        });
      }
      out.push({
        method: "lsq",
        move,
        onto,
        ...(hasShorthand ? { chain, range } : { matches }),
        matchType: (matchType as "all" | "main" | "ca" | undefined) ?? undefined,
      });
      return;
    }

    // Unknown method already errored above; skip emission.
  });

  return out;
}

function validateLsqMatch(
  raw: unknown,
  path: string,
  errors: SceneValidationError[],
): SceneLsqMatch | null {
  if (!isObject(raw)) {
    errors.push({ path, message: "must be a mapping" });
    return null;
  }
  const refChain = strField(raw, "refChain", "", errors, path);
  const refRange = rangeField(raw, "refRange", errors, path) ?? "";
  const movChain = strField(raw, "movChain", "", errors, path);
  const movRange = rangeField(raw, "movRange", errors, path) ?? "";
  if (!refChain) errors.push({ path: `${path}.refChain`, message: "required" });
  if (!movChain) errors.push({ path: `${path}.movChain`, message: "required" });
  for (const [field, val] of [["refRange", refRange], ["movRange", movRange]] as const) {
    if (!val) {
      errors.push({ path: `${path}.${field}`, message: "required" });
    } else if (!RANGE_RE.test(val)) {
      errors.push({
        path: `${path}.${field}`,
        message: `must be "start-end", got "${val}"`,
      });
    } else {
      const m = RANGE_RE.exec(val)!;
      const s = parseInt(m[1], 10);
      const e = parseInt(m[2], 10);
      if (e < s) {
        errors.push({
          path: `${path}.${field}`,
          message: `end (${e}) must be >= start (${s})`,
        });
      }
    }
  }
  return { refChain, refRange, movChain, movRange };
}

/**
 * Validate a list of dictionary file-name references — used for both
 * the top-level `globalDictionaries:` block and per-element
 * `dictionaries:` lists. Each entry must be a string that names an
 * existing file ref in the `files:` block, and that ref must have
 * `kind: dictionary`.
 */
function validateDictionaryNameList(
  raw: unknown,
  files: SceneFileRef[],
  path: string,
  errors: SceneValidationError[],
): string[] | undefined {
  if (!Array.isArray(raw)) {
    errors.push({ path, message: "must be a sequence of file names" });
    return undefined;
  }
  const dictNames = new Set(
    files.filter((f) => f.kind === "dictionary").map((f) => f.name),
  );
  const knownNames = new Set(files.map((f) => f.name));
  const out: string[] = [];
  raw.forEach((entry, i) => {
    const p = `${path}[${i}]`;
    if (typeof entry !== "string") {
      errors.push({ path: p, message: `must be a string, got ${typeof entry}` });
      return;
    }
    if (files.length === 0) {
      // No files block at all; can't cross-check. Leave as-is.
      out.push(entry);
      return;
    }
    if (!knownNames.has(entry)) {
      errors.push({
        path: p,
        message: `unknown file "${entry}" (not in top-level files block)`,
      });
      return;
    }
    if (!dictNames.has(entry)) {
      errors.push({
        path: p,
        message: `file "${entry}" is not a dictionary (missing \`kind: dictionary\`)`,
      });
      return;
    }
    out.push(entry);
  });
  return out;
}

function validateDomains(
  raw: unknown,
  errors: SceneValidationError[],
): SceneDomain[] | undefined {
  if (!Array.isArray(raw)) {
    errors.push({ path: "domains", message: "must be a sequence" });
    return undefined;
  }
  const seen = new Set<string>();
  const out: SceneDomain[] = [];
  raw.forEach((entry, i) => {
    const p = `domains[${i}]`;
    if (!isObject(entry)) {
      errors.push({ path: p, message: "must be a mapping" });
      return;
    }
    const e = entry as Record<string, unknown>;
    const name = strField(entry, "name", "", errors, p);
    const color = strField(entry, "color", "", errors, p);
    if (!name) errors.push({ path: `${p}.name`, message: "required" });
    if (!color) {
      errors.push({ path: `${p}.color`, message: "required" });
    } else if (!isHexColor(color)) {
      errors.push({
        path: `${p}.color`,
        message: `must be hex like "#rrggbb" or "#rrggbbaa", got "${color}"`,
      });
    }
    if (name) {
      if (seen.has(name)) {
        errors.push({ path: `${p}.name`, message: `duplicate domain name "${name}"` });
      }
      seen.add(name);
    }

    const domain: SceneDomain = { name, color };
    if ("selection" in e) {
      // Preferred CID form.
      const sel = strField(entry, "selection", "", errors, p);
      if (!sel) {
        errors.push({ path: `${p}.selection`, message: "must be a non-empty CID string" });
      } else {
        domain.selection = sel;
      }
    } else {
      // Legacy chain + optional range (omitted ⇒ whole chain).
      const chain = parseChainField(entry, `${p}.chain`, errors);
      if (chain === undefined) {
        errors.push({ path: `${p}.chain`, message: "required (or use `selection`)" });
      } else {
        domain.chain = chain;
      }
      if ("range" in e) {
        const range = rangeField(entry, "range", errors, p) ?? "";
        if (!range || !RANGE_RE.test(range)) {
          errors.push({ path: `${p}.range`, message: `must be "start-end", got "${range}"` });
        } else {
          const m = RANGE_RE.exec(range)!;
          const start = parseInt(m[1], 10);
          const end = parseInt(m[2], 10);
          if (end < start) {
            errors.push({ path: `${p}.range`, message: `end (${end}) must be >= start (${start})` });
          } else {
            domain.range = range;
          }
        }
      }
    }
    out.push(domain);
  });
  return out;
}

function validateElements(
  raw: unknown,
  files: SceneFileRef[],
  errors: SceneValidationError[],
): SceneElement[] | undefined {
  if (!Array.isArray(raw)) {
    errors.push({ path: "elements", message: "must be a sequence" });
    return undefined;
  }
  const fileNames = new Set(files.map((f) => f.name));
  const out: SceneElement[] = [];
  raw.forEach((entry, i) => {
    const p = `elements[${i}]`;
    if (!isObject(entry)) {
      errors.push({ path: p, message: "must be a mapping" });
      return;
    }
    const file = strField(entry, "file", "", errors, p);
    if (!file) {
      errors.push({ path: `${p}.file`, message: "required" });
    } else if (files.length > 0 && !fileNames.has(file)) {
      errors.push({
        path: `${p}.file`,
        message: `unknown file "${file}" (not in top-level files block)`,
      });
    }

    const el: SceneElement = { file };
    if ("dictionaries" in entry) {
      const dicts = validateDictionaryNameList(
        (entry as Record<string, unknown>).dictionaries,
        files,
        `${p}.dictionaries`,
        errors,
      );
      if (dicts) el.dictionaries = dicts;
    }
    if ("representations" in entry) {
      const reps = (entry as Record<string, unknown>).representations;
      if (!Array.isArray(reps)) {
        errors.push({ path: `${p}.representations`, message: "must be a sequence" });
      } else {
        el.representations = reps
          .map((r, j) => validateRepresentation(r, `${p}.representations[${j}]`, errors))
          .filter((r): r is SceneRepresentation => r !== null);
      }
    }
    out.push(el);
  });
  return out;
}

function validateMaps(
  raw: unknown,
  files: SceneFileRef[],
  errors: SceneValidationError[],
): SceneMap[] | undefined {
  if (!Array.isArray(raw)) {
    errors.push({ path: "maps", message: "must be a sequence" });
    return undefined;
  }
  const fileNames = new Set(files.map((f) => f.name));
  const kindByName = new Map(files.map((f) => [f.name, f.kind ?? "coordinates"]));
  const seen = new Set<string>();
  const out: SceneMap[] = [];
  raw.forEach((entry, i) => {
    const p = `maps[${i}]`;
    if (!isObject(entry)) {
      errors.push({ path: p, message: "must be a mapping" });
      return;
    }
    const name = strField(entry, "name", "", errors, p);
    if (!name) {
      errors.push({ path: `${p}.name`, message: "required" });
    } else if (seen.has(name)) {
      errors.push({ path: `${p}.name`, message: `duplicate map name "${name}"` });
    } else {
      seen.add(name);
    }

    const file = strField(entry, "file", "", errors, p);
    const fileKind = file ? kindByName.get(file) : undefined;
    if (!file) {
      errors.push({ path: `${p}.file`, message: "required" });
    } else if (files.length > 0 && !fileNames.has(file)) {
      errors.push({
        path: `${p}.file`,
        message: `unknown file "${file}" (not in top-level files block)`,
      });
    } else if (fileNames.has(file) && fileKind !== "mtz" && fileKind !== "map") {
      errors.push({
        path: `${p}.file`,
        message: `"${file}" must be a file with kind: "mtz" or "map"`,
      });
    }

    // `kind: "map"` files (real-space CCP4 maps / masks) take no columns —
    // they're read directly. Only `kind: "mtz"` files need a column spec.
    const isMapFile = fileKind === "map";
    const columnsRaw = (entry as Record<string, unknown>).columns;
    let columns: SceneMapColumns | undefined;
    if (isMapFile) {
      if (columnsRaw !== undefined) {
        errors.push({
          path: `${p}.columns`,
          message: `not allowed for kind: "map" file "${file}" (read directly, no columns)`,
        });
      }
    } else if (!isObject(columnsRaw)) {
      errors.push({ path: `${p}.columns`, message: "required mapping (F + PHI minimum)" });
    } else {
      columns = {
        F: optStr(columnsRaw, "F", `${p}.columns.F`, errors),
        PHI: optStr(columnsRaw, "PHI", `${p}.columns.PHI`, errors),
        Fobs: optStr(columnsRaw, "Fobs", `${p}.columns.Fobs`, errors),
        SigFobs: optStr(columnsRaw, "SigFobs", `${p}.columns.SigFobs`, errors),
        FreeR: optStr(columnsRaw, "FreeR", `${p}.columns.FreeR`, errors),
        useWeight: optBool(columnsRaw, "useWeight", `${p}.columns.useWeight`, errors),
        calcStructFact: optBool(columnsRaw, "calcStructFact", `${p}.columns.calcStructFact`, errors),
      };
      // Drop undefined fields so the round-trip stays clean.
      (Object.keys(columns) as (keyof SceneMapColumns)[]).forEach((k) => {
        if (columns![k] === undefined) delete columns![k];
      });
      // Need F + PHI unless calcStructFact (then Moorhen computes them
      // from Fobs/SigFobs/FreeR).
      const haveFP = !!columns.F && !!columns.PHI;
      const haveCalc = columns.calcStructFact && !!columns.Fobs && !!columns.SigFobs;
      if (!haveFP && !haveCalc) {
        errors.push({
          path: `${p}.columns`,
          message: "must set F + PHI (or calcStructFact + Fobs + SigFobs)",
        });
      }
    }

    const map: SceneMap = columns ? { name, file, columns } : { name, file };

    if ("isMask" in entry) {
      map.isMask = optBool(entry, "isMask", `${p}.isMask`, errors);
    }
    if ("isDifference" in entry) {
      map.isDifference = optBool(entry, "isDifference", `${p}.isDifference`, errors);
    }
    if ("contourLevel" in entry) {
      map.contourLevel = optNum(entry, "contourLevel", `${p}.contourLevel`, errors);
    }
    if ("radius" in entry) {
      map.radius = optNum(entry, "radius", `${p}.radius`, errors);
    }
    if ("alpha" in entry) {
      map.alpha = optNum(entry, "alpha", `${p}.alpha`, errors);
    }
    if ("style" in entry) {
      const s = optStr(entry, "style", `${p}.style`, errors);
      if (s && s !== "lines" && s !== "solid" && s !== "lit-lines") {
        errors.push({
          path: `${p}.style`,
          message: `must be "lines", "solid", or "lit-lines", got "${s}"`,
        });
      } else if (s) {
        map.style = s as SceneMap["style"];
      }
    }
    for (const k of ["colour", "positiveColour", "negativeColour"] as const) {
      if (k in entry) {
        const c = optStr(entry, k, `${p}.${k}`, errors);
        if (c) {
          if (!isHexColor(c)) {
            errors.push({
              path: `${p}.${k}`,
              message: `must be hex like "#rrggbb" or "#rrggbbaa", got "${c}"`,
            });
          } else {
            map[k] = c;
          }
        }
      }
    }
    if ("visible" in entry) {
      map.visible = optBool(entry, "visible", `${p}.visible`, errors);
    }

    out.push(map);
  });
  return out;
}

function validateRepresentation(
  raw: unknown,
  path: string,
  errors: SceneValidationError[],
): SceneRepresentation | null {
  if (!isObject(raw)) {
    errors.push({ path, message: "must be a mapping" });
    return null;
  }
  const style = strField(raw, "style", "", errors, path);
  if (!style) {
    errors.push({ path: `${path}.style`, message: "required" });
  }
  const rep: SceneRepresentation = { style };
  if ("selection" in raw) {
    rep.selection = optStr(raw, "selection", `${path}.selection`, errors);
  }
  if ("colour" in raw) {
    const c = (raw as Record<string, unknown>).colour;
    rep.colour = validateColour(c, `${path}.colour`, errors) ?? undefined;
  }
  if ("alpha" in raw) {
    const a = optNum(raw, "alpha", `${path}.alpha`, errors);
    if (a !== undefined) {
      if (a < 0 || a > 1) {
        errors.push({ path: `${path}.alpha`, message: "must be in [0, 1]" });
      } else {
        rep.alpha = a;
      }
    }
  }
  return rep;
}

function validateCentre(
  raw: unknown,
  files: SceneFileRef[],
  errors: SceneValidationError[],
): SceneCentre | undefined {
  if (raw === undefined) return undefined;
  if (!isObject(raw)) {
    errors.push({ path: "view.centre", message: "must be a mapping { file, selection? }" });
    return undefined;
  }
  // centre takes only file + selection. Reject anything else loudly — it's a
  // typo (e.g. "-selection"), and silently dropping it gives a wrong-but-quiet
  // centre (the whole molecule instead of the intended selection).
  for (const k of Object.keys(raw as Record<string, unknown>)) {
    if (k !== "file" && k !== "selection") {
      errors.push({
        path: `view.centre.${k}`,
        message: `unknown key "${k}" — centre takes only "file" and "selection"`,
      });
    }
  }
  const file = strField(raw, "file", "", errors, "view.centre");
  if (!file) {
    errors.push({ path: "view.centre.file", message: "required" });
    return undefined;
  }
  if (files.length > 0 && !files.some((f) => f.name === file)) {
    errors.push({
      path: "view.centre.file",
      message: `unknown file "${file}" (not in top-level files block)`,
    });
    return undefined;
  }
  const centre: SceneCentre = { file };
  if ("selection" in (raw as Record<string, unknown>)) {
    const sel = optStr(raw, "selection", "view.centre.selection", errors);
    if (sel) centre.selection = sel;
  }
  return centre;
}

function validateColour(
  raw: unknown,
  path: string,
  errors: SceneValidationError[],
): SceneColour | null {
  if (Array.isArray(raw)) {
    // Per-selection colour list: [{ selection, colour }, ...].
    const list: SceneColourSelection[] = [];
    raw.forEach((entry, i) => {
      const ep = `${path}[${i}]`;
      if (!isObject(entry)) {
        errors.push({ path: ep, message: "must be a mapping { selection, colour }" });
        return;
      }
      const selection = strField(entry, "selection", "", errors, ep);
      if (!selection) errors.push({ path: `${ep}.selection`, message: "required" });
      const colour = strField(entry, "colour", "", errors, ep);
      if (!colour) {
        errors.push({ path: `${ep}.colour`, message: "required" });
      } else if (!isHexColor(colour)) {
        errors.push({ path: `${ep}.colour`, message: "must be hex (#rrggbb or #rrggbbaa)" });
      }
      if (selection && colour && isHexColor(colour)) list.push({ selection, colour });
    });
    return list.length > 0 ? list : null;
  }
  if (typeof raw === "string") {
    if (isHexColor(raw)) return raw;
    if (isSceneNamedColour(raw as SceneColour)) return raw as SceneColour;
    errors.push({
      path,
      message: `unknown colour "${raw}" — expected hex (#rrggbb or #rrggbbaa) or named scheme`,
    });
    return null;
  }
  if (isObject(raw) && "raw" in raw) {
    const r = (raw as { raw: unknown }).raw;
    if (!isObject(r)) {
      errors.push({ path: `${path}.raw`, message: "must be a mapping" });
      return null;
    }
    const ruleType = strField(r, "ruleType", "", errors, `${path}.raw`);
    const args = (r as Record<string, unknown>).args;
    if (!ruleType) {
      errors.push({ path: `${path}.raw.ruleType`, message: "required" });
    }
    if (!Array.isArray(args)) {
      errors.push({ path: `${path}.raw.args`, message: "must be a sequence" });
      return null;
    }
    return {
      raw: {
        ruleType,
        args: args as (string | number)[],
        isMultiColourRule: optBool(r, "isMultiColourRule", `${path}.raw.isMultiColourRule`, errors),
        applyColourToNonCarbonAtoms: optBool(
          r,
          "applyColourToNonCarbonAtoms",
          `${path}.raw.applyColourToNonCarbonAtoms`,
          errors,
        ),
      },
    };
  }
  errors.push({
    path,
    message:
      "must be a hex string, named scheme, per-selection list, or { raw: ... }",
  });
  return null;
}

// --------------------------------------------------------------------------
// Serialise: MoorhenScene → YAML string
// --------------------------------------------------------------------------

/**
 * Serialise a scene back to YAML. Field order is canonical (top-level
 * keys appear in the order declared in the type), and empty/undefined
 * fields are omitted so round-trip is stable.
 */
export function serialiseScene(scene: MoorhenScene): string {
  return YAML.stringify(buildOrderedScene(scene), { lineWidth: 0 });
}

/**
 * Like serialiseScene, but attaches per-file comments above each entry
 * in the files: sequence. Used by the lifter to record source URLs
 * without polluting the schema.
 *
 * fileComments is keyed by SceneFileRef.name; missing entries produce
 * no comment.
 */
export function serialiseSceneWithComments(
  scene: MoorhenScene,
  fileComments: Record<string, string> = {},
): string {
  const ordered = buildOrderedScene(scene);
  const doc = new YAML.Document(ordered);

  if (scene.files && scene.files.length > 0) {
    const filesNode = doc.get("files") as YAML.YAMLSeq | undefined;
    if (filesNode && Array.isArray(filesNode.items)) {
      scene.files.forEach((f, i) => {
        const comment = fileComments[f.name];
        const item = filesNode.items[i];
        if (comment && item && typeof item === "object") {
          (item as { commentBefore?: string }).commentBefore = ` ${comment}`;
        }
      });
    }
  }

  return doc.toString({ lineWidth: 0 });
}

function buildOrderedScene(scene: MoorhenScene): Record<string, unknown> {
  const ordered: Record<string, unknown> = {};
  ordered.scene = scene.scene;
  ordered.version = scene.version;
  if (scene.authoredIn && hasAnyValue(scene.authoredIn)) {
    ordered.authoredIn = stripUndefined(scene.authoredIn);
  }
  if (scene.files && scene.files.length > 0) {
    ordered.files = scene.files.map((f) => stripUndefined(f));
  }
  if (scene.superpose && scene.superpose.length > 0) {
    ordered.superpose = scene.superpose.map((sp) => stripUndefined(sp));
  }
  if (scene.globalDictionaries && scene.globalDictionaries.length > 0) {
    ordered.globalDictionaries = scene.globalDictionaries;
  }
  if (scene.domains && scene.domains.length > 0) {
    ordered.domains = scene.domains;
  }
  if (scene.elements && scene.elements.length > 0) {
    ordered.elements = scene.elements.map((el) => {
      const out: Record<string, unknown> = { file: el.file };
      if (el.dictionaries && el.dictionaries.length > 0) {
        out.dictionaries = el.dictionaries;
      }
      if (el.representations && el.representations.length > 0) {
        out.representations = el.representations.map((r) => stripUndefined(r));
      }
      return out;
    });
  }
  if (scene.maps && scene.maps.length > 0) {
    ordered.maps = scene.maps.map((m) => {
      const out: Record<string, unknown> = {
        name: m.name,
        file: m.file,
      };
      // columns only for MTZ-backed maps; real-space map/mask refs omit it.
      if (m.columns !== undefined) out.columns = stripUndefined(m.columns);
      // Optional render-state fields: emit only the ones the lifter
      // (or hand-author) set, so the YAML stays small.
      const optional: (keyof SceneMap)[] = [
        "isMask", "isDifference", "contourLevel", "radius", "alpha", "style",
        "colour", "positiveColour", "negativeColour", "visible",
      ];
      for (const k of optional) {
        if (m[k] !== undefined) out[k] = m[k];
      }
      return out;
    });
  }
  if (scene.activeMap) {
    ordered.activeMap = scene.activeMap;
  }
  if (scene.view && hasAnyValue(scene.view)) {
    ordered.view = stripUndefined(scene.view);
  }
  if (scene.resolver && scene.resolver.onMissingResidues) {
    ordered.resolver = { onMissingResidues: scene.resolver.onMissingResidues };
  }
  return ordered;
}

// --------------------------------------------------------------------------
// Small utilities
// --------------------------------------------------------------------------

function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null && !Array.isArray(x);
}

function isHexColor(s: string): boolean {
  // 6-hex (#rrggbb) or 8-hex with alpha (#rrggbbaa). Moorhen's own
  // MoorhenColourRule.parseHexToRgba accepts both, and Moorhen's
  // default per-chain colour rules come out as 8-hex, so the scene
  // format mirrors that.
  return /^#[0-9a-fA-F]{6}([0-9a-fA-F]{2})?$/.test(s);
}

function hasAnyValue(o: object): boolean {
  return Object.values(o as Record<string, unknown>).some((v) => v !== undefined);
}

function stripUndefined<T extends object>(o: T): Partial<T> {
  const out: Record<string, unknown> = {};
  for (const [k, v] of Object.entries(o as Record<string, unknown>)) {
    if (v !== undefined) out[k] = v;
  }
  return out as Partial<T>;
}

function strField(
  obj: Record<string, unknown>,
  key: string,
  fallback: string,
  errors: SceneValidationError[],
  parent = "",
): string {
  const v = obj[key];
  if (v === undefined) return fallback;
  if (typeof v !== "string") {
    errors.push({
      path: parent ? `${parent}.${key}` : key,
      message: `must be a string, got ${typeof v}`,
    });
    return fallback;
  }
  return v;
}

/**
 * Read a residue-range field. Accepts either a "start-end" string or a
 * bare integer (treated as the single-residue range "N-N"). Returns the
 * normalised string form so downstream code only sees one shape. Pushes
 * an error and returns "" if the field is present but the wrong type.
 */
function rangeField(
  obj: Record<string, unknown>,
  key: string,
  errors: SceneValidationError[],
  parent = "",
): string | undefined {
  const v = obj[key];
  if (v === undefined) return undefined;
  if (typeof v === "number" && Number.isFinite(v) && Number.isInteger(v)) {
    return `${v}-${v}`;
  }
  if (typeof v === "string") return v;
  errors.push({
    path: parent ? `${parent}.${key}` : key,
    message: `must be a string or integer, got ${typeof v}`,
  });
  return undefined;
}

function numField(
  obj: Record<string, unknown>,
  key: string,
  fallback: number,
  errors: SceneValidationError[],
  parent = "",
): number {
  const v = obj[key];
  if (v === undefined) return fallback;
  if (typeof v !== "number") {
    errors.push({
      path: parent ? `${parent}.${key}` : key,
      message: `must be a number, got ${typeof v}`,
    });
    return fallback;
  }
  return v;
}

function optStr(
  obj: Record<string, unknown>,
  key: string,
  path: string,
  errors: SceneValidationError[],
): string | undefined {
  const v = obj[key];
  if (v === undefined) return undefined;
  if (typeof v !== "string") {
    errors.push({ path, message: `must be a string, got ${typeof v}` });
    return undefined;
  }
  return v;
}

function optNum(
  obj: Record<string, unknown>,
  key: string,
  path: string,
  errors: SceneValidationError[],
): number | undefined {
  const v = obj[key];
  if (v === undefined) return undefined;
  if (typeof v !== "number") {
    errors.push({ path, message: `must be a number, got ${typeof v}` });
    return undefined;
  }
  return v;
}

function optBool(
  obj: Record<string, unknown>,
  key: string,
  path: string,
  errors: SceneValidationError[],
): boolean | undefined {
  const v = obj[key];
  if (v === undefined) return undefined;
  if (typeof v !== "boolean") {
    errors.push({ path, message: `must be a boolean, got ${typeof v}` });
    return undefined;
  }
  return v;
}

function tupleField<N extends number>(
  obj: Record<string, unknown>,
  key: string,
  length: N,
  path: string,
  errors: SceneValidationError[],
): (N extends 3 ? [number, number, number] : [number, number, number, number]) | undefined {
  const v = obj[key];
  if (v === undefined) return undefined;
  if (!Array.isArray(v) || v.length !== length) {
    errors.push({ path, message: `must be an array of ${length} numbers` });
    return undefined;
  }
  if (!v.every((x) => typeof x === "number")) {
    errors.push({ path, message: `all elements must be numbers` });
    return undefined;
  }
  return v as never;
}

/**
 * Parse a `chain:` field accepting either a string (single chain id or
 * the "*" wildcard) or a sequence of strings.
 *
 * Returns `undefined` only when the field is missing — that lets the
 * caller emit its own "required" error so error messages don't double
 * up. Returns an empty string for typed-but-empty values, so the
 * required-check still fires uniformly.
 */
function parseChainField(
  obj: Record<string, unknown>,
  path: string,
  errors: SceneValidationError[],
): string | string[] | undefined {
  const v = obj["chain"];
  if (v === undefined) return undefined;
  if (typeof v === "string") return v;
  if (Array.isArray(v)) {
    if (v.length === 0) {
      errors.push({ path, message: "chain list must not be empty" });
      return "";
    }
    if (!v.every((x): x is string => typeof x === "string")) {
      errors.push({ path, message: "chain list entries must all be strings" });
      return "";
    }
    return v;
  }
  errors.push({ path, message: `must be a string or sequence of strings, got ${typeof v}` });
  return "";
}

// Re-export the guards so callers don't need two imports.
export { isSceneHexColour, isSceneNamedColour, isSceneRawColour };
