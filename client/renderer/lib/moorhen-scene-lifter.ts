/**
 * Moorhen Scene lifter: Moorhen state → MoorhenScene.
 *
 * Captures what's currently shown in the viewer as a scene the user can
 * download, hand-edit, and re-apply. v0 scope mirrors the user's actual
 * authoring workflow:
 *
 *  - Always lifted: files, camera, per-molecule representation list with
 *    their CID selections.
 *  - Lifted if recognisable: colour rules — only named multi-rule schemes
 *    (b-factor, af2-plddt, secondary-structure, ...), single-colour hex,
 *    and the by-domain pipe-delimited shape that parses back cleanly.
 *  - Unknown colour rules: emitted under the `{ raw: { ruleType, args } }`
 *    escape hatch so they round-trip losslessly even when unreadable.
 *  - Dropped: bond options, m2t params, lighting, SSAO. The user can add
 *    these by hand in the yaml if they want them.
 *
 * No serialisation in this file — it returns a MoorhenScene. The
 * companion serialiser (moorhen-scene.ts) handles YAML conversion.
 * The "download scene" UI button glues the two together.
 */

import type { moorhen } from "moorhen/types/moorhen";

import { extractFileIdFromUniqueId } from "./moorhen-view-state";
import {
  MoorhenScene,
  SCENE_SCHEMA_VERSION,
  SceneColour,
  SceneDomain,
  SceneElement,
  SceneFileRef,
  SceneRepresentation,
  SceneView,
} from "../types/moorhen-scene";

// --------------------------------------------------------------------------
// Public API
// --------------------------------------------------------------------------

export interface LiftCtx {
  /** The molecules currently in Moorhen. */
  molecules: moorhen.Molecule[];
  /** glRef state (origin, quat, zoom, clip*, fog*). Pass `store.getState().glRef`. */
  glRef: {
    origin: number[] | Float32Array;
    quat: number[] | Float32Array;
    zoom: number;
    clipStart?: number;
    clipEnd?: number;
    fogStart?: number;
    fogEnd?: number;
  };
  /** Optional: the ccp4i2 project UUID. If set, file refs get both
   *  projectId and (derivable) fileId. */
  projectId?: string;
  /** Optional: human-readable project name (advisory). */
  projectName?: string;
  /** Optional: scene name. Default "scene-<isoDate>". */
  sceneName?: string;
  /** Optional: who is authoring (email/username). */
  createdBy?: string;
}

/**
 * Capture the current Moorhen state as a MoorhenScene. The returned
 * object is independently serialisable by serialiseScene().
 *
 * Per-rep `uniqueId` strings (the loader URL) survive as carried-over
 * properties on the returned SceneFileRef under a non-enumerable internal
 * symbol, used by the serialiser to emit a YAML comment naming the
 * original source. Callers that go through serialiseSceneWithLiftHints
 * see those comments; everyone else sees a clean MoorhenScene.
 */
export function liftScene(ctx: LiftCtx): MoorhenScene {
  const scene: MoorhenScene = {
    scene: ctx.sceneName ?? `scene-${new Date().toISOString().slice(0, 19).replace(/[T:]/g, "-")}`,
    version: SCENE_SCHEMA_VERSION,
  };

  scene.authoredIn = {
    projectId: ctx.projectId,
    projectName: ctx.projectName,
    createdAt: new Date().toISOString(),
    createdBy: ctx.createdBy,
  };
  if (!hasAnyValue(scene.authoredIn as Record<string, unknown>)) {
    delete scene.authoredIn;
  }

  scene.files = ctx.molecules.map((mol) => liftFileRef(mol, ctx));
  if (scene.files.length === 0) delete scene.files;

  // Lift dictionaries scoped to each molecule. For each loaded molecule
  // we enumerate non-standard residue types and ask Moorhen for the
  // stored dict text. We emit them as `kind: dictionary` file refs with
  // an inline `cifText:`, keyed by molecule + comp_id so the same dict
  // applied to two molecules ends up as two separate entries (matches
  // the per-molecule scoping semantics).
  const liftedDicts: { mol: moorhen.Molecule; molFileName: string; refs: { name: string; comp_id: string }[] }[] = [];
  ctx.molecules.forEach((mol, i) => {
    const molFileName = scene.files?.[i]?.name ?? `mol${i}`;
    const liftedRefs = liftDictionariesForMolecule(mol, molFileName, scene.files!);
    if (liftedRefs.length > 0) {
      liftedDicts.push({ mol, molFileName, refs: liftedRefs });
    }
  });

  const view = liftView(ctx.glRef);
  if (view) scene.view = view;

  const elements = ctx.molecules
    .map((mol, i) => liftElement(mol, scene.files?.[i]?.name ?? `mol${i}`))
    .filter((el): el is SceneElement => el !== null);

  // Attach the lifted dict refs to the element for their molecule.
  for (const { molFileName, refs } of liftedDicts) {
    const el = elements.find((e) => e.file === molFileName);
    if (el) el.dictionaries = refs.map((r) => r.name);
  }

  if (elements.length > 0) scene.elements = elements;

  return scene;
}

/**
 * Enumerate non-standard residues in a molecule and emit a SceneFileRef
 * per (molecule, comp_id) pair where Moorhen has a stored dictionary
 * for that residue. Appends the new refs to `allFiles` (the scene's
 * `files:` block) and returns the (name, comp_id) list for the caller
 * to use when wiring them into the element.
 *
 * Naming: dict refs get names like `dict-<molFileName>-<comp_id>` to
 * stay readable while staying unique across multiple molecules with
 * same-named ligands.
 */
function liftDictionariesForMolecule(
  mol: moorhen.Molecule,
  molFileName: string,
  allFiles: SceneFileRef[],
): { name: string; comp_id: string }[] {
  const ligands = mol.ligands ?? [];
  // Dedupe by comp_id — multiple instances of the same ligand share a dict.
  const compIds = Array.from(new Set(ligands.map((l) => l.resName).filter(Boolean)));
  const out: { name: string; comp_id: string }[] = [];

  for (const comp_id of compIds) {
    if (STANDARD_MONOMERS.has(comp_id)) continue;
    let dictText = "";
    try {
      dictText = mol.getDict(comp_id) || "";
    } catch {
      // getDict may throw for comp_ids without a stored dict — that's
      // fine, just skip them.
      continue;
    }
    if (!dictText) continue;

    const name = `dict-${molFileName}-${comp_id}`;
    // Skip if a ref with this name already exists (shouldn't happen but
    // defensive against pathological dupe-suppression).
    if (allFiles.some((f) => f.name === name)) continue;

    allFiles.push({
      name,
      kind: "dictionary",
      cifText: dictText,
    });
    out.push({ name, comp_id });
  }
  return out;
}

/**
 * Residue types Coot already knows about — we shouldn't emit dicts for
 * these because they're built-in. Conservative list: 20 amino acids,
 * standard nucleotides, water/ions/cofactors that come with the monomer
 * library. If a structure carries a non-standard variant (e.g. SEP,
 * PTR) we want to capture *that* dict, so those are intentionally absent.
 */
const STANDARD_MONOMERS: ReadonlySet<string> = new Set([
  // Standard amino acids
  "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
  // Standard nucleotides
  "A", "C", "G", "T", "U", "DA", "DC", "DG", "DT", "DU",
  // Waters and very common ions/cofactors
  "HOH", "WAT", "H2O", "OH2",
]);

/**
 * Per-file annotations the lifter can produce that don't fit in the
 * MoorhenScene type itself (e.g. "the loader URL this file came from",
 * which we want emitted as a YAML comment by the serialiser).
 *
 * Returned as a parallel structure keyed by file name. The download-UI
 * glue passes these to the serialiser; pure parse/serialise flows ignore
 * them.
 */
export interface SceneLiftHints {
  /** Comments to attach above each `files[i]` entry, keyed by SceneFileRef.name. */
  fileComments: Record<string, string>;
}

/**
 * Same as liftScene, but also returns hints for richer YAML output
 * (currently: a comment above each file entry naming its uniqueId).
 */
export function liftSceneWithHints(ctx: LiftCtx): {
  scene: MoorhenScene;
  hints: SceneLiftHints;
} {
  const scene = liftScene(ctx);
  const fileComments: Record<string, string> = {};
  ctx.molecules.forEach((mol, i) => {
    const fr = scene.files?.[i];
    if (fr && mol.uniqueId) {
      fileComments[fr.name] = `loaded from ${mol.uniqueId}`;
    }
  });
  return { scene, hints: { fileComments } };
}

/**
 * Lift a scene AND a parallel asset map suitable for bundling into a
 * `.scene.zip`. Compared to liftScene/liftSceneWithHints, this:
 *
 *  - For each loaded coordinate molecule, asks Moorhen for the current
 *    coord text via mol.getAtoms() and emits a `bundle:` file ref. The
 *    bytes go into the asset map under `assets/<safeName>.<ext>`.
 *  - Already-emitted `cifText:` dictionary refs are moved into bundle
 *    assets too, because once we have a zip there's no reason to inline
 *    them in the YAML.
 *
 * The returned scene is self-contained: serialised to YAML and zipped
 * up alongside its asset map, it round-trips on any machine.
 */
export async function liftSceneToBundle(ctx: LiftCtx): Promise<{
  scene: MoorhenScene;
  hints: SceneLiftHints;
  assets: Map<string, ArrayBuffer>;
}> {
  const { scene, hints } = liftSceneWithHints(ctx);
  const assets = new Map<string, ArrayBuffer>();

  // 1. Replace each coord file ref with a bundle: ref carrying the
  //    current coords. Skip refs the lifter wrote as PDB IDs / project
  //    fileIds — those are already self-resolving and don't need
  //    inlining.
  const filesByName = new Map<string, SceneFileRef>();
  for (const f of scene.files ?? []) filesByName.set(f.name, f);

  for (const mol of ctx.molecules) {
    const ref = scene.files?.find(
      (f) => f.name === sanitiseName(mol.name) || f.name === `mol${ctx.molecules.indexOf(mol)}`,
    );
    if (!ref) continue;
    if (ref.kind === "dictionary") continue; // handled below
    // Only inline coords for refs that aren't already pdb/fileId/url —
    // those are portable on their own. Bundle is for "this only exists
    // in memory or on local disk".
    const alreadyPortable = !!ref.pdb || ref.fileId !== undefined || !!ref.url;
    if (alreadyPortable) continue;

    try {
      const coordText = await mol.getAtoms();
      if (!coordText) continue;
      const ext = (mol.coordsFormat === "mmcif" || mol.coordsFormat === "mmjson") ? "cif" : "pdb";
      const assetPath = `assets/${ref.name}.${ext}`;
      assets.set(assetPath, new TextEncoder().encode(coordText).buffer);
      // Replace the unresolvable path: form with the portable bundle: form.
      ref.bundle = assetPath;
      delete ref.path;
      // Preserve the original path as a comment so the round-tripped
      // YAML still records where the user originally pulled this from.
      if (mol.uniqueId && !hints.fileComments[ref.name]) {
        hints.fileComments[ref.name] = `originally from ${mol.uniqueId}`;
      }
    } catch (err) {
      console.warn(`[lift] could not bundle coords for ${ref.name}:`, err);
    }
  }

  // 2. Move cifText dicts into bundle assets. Cleaner YAML + smaller
  //    when the same dict appears via multiple scenes.
  for (const ref of scene.files ?? []) {
    if (ref.kind !== "dictionary" || !ref.cifText) continue;
    const assetPath = `assets/${ref.name}.cif`;
    assets.set(assetPath, new TextEncoder().encode(ref.cifText).buffer);
    ref.bundle = assetPath;
    delete ref.cifText;
  }

  return { scene, hints, assets };
}

// --------------------------------------------------------------------------
// Files
// --------------------------------------------------------------------------

function liftFileRef(mol: moorhen.Molecule, ctx: LiftCtx): SceneFileRef {
  const name = sanitiseName(mol.name) || `mol${mol.molNo ?? "?"}`;

  // Project-internal: ccp4i2 file id embedded in the loader URL pattern.
  const fileId = extractFileIdFromUniqueId(mol.uniqueId ?? "");
  if (fileId !== null && ctx.projectId) {
    return {
      name,
      projectId: ctx.projectId,
      projectName: ctx.projectName,
      fileId,
    };
  }

  // Plain URL fallback for non-ccp4i2 loaders (e.g. PDBe direct fetches).
  if (mol.uniqueId && /^https?:\/\//.test(mol.uniqueId)) {
    return { name, url: mol.uniqueId };
  }

  // Path fallback for local fetches that don't look like URLs.
  if (mol.uniqueId) {
    return { name, path: mol.uniqueId };
  }

  // Last resort: a placeholder url that the user will edit.
  return { name, url: "TODO: source URL for this molecule" };
}

// --------------------------------------------------------------------------
// View
// --------------------------------------------------------------------------

function liftView(glRef: LiftCtx["glRef"]): SceneView | undefined {
  if (!glRef) return undefined;
  // glRef.origin / glRef.quat can be Float32Array — coerce to arrays so
  // they JSON/YAML-stringify as sequences not as objects.
  const origin = toArray3(glRef.origin);
  const quat = toArray4(glRef.quat);
  const view: SceneView = {};
  if (origin) view.origin = origin;
  if (quat) view.quat = quat;
  if (typeof glRef.zoom === "number") view.zoom = round(glRef.zoom, 4);
  // clip/fog: only include if set; the apply path will use Moorhen
  // defaults for anything omitted.
  if (typeof glRef.clipStart === "number") view.clipStart = glRef.clipStart;
  if (typeof glRef.clipEnd === "number") view.clipEnd = glRef.clipEnd;
  if (typeof glRef.fogStart === "number") view.fogStart = glRef.fogStart;
  if (typeof glRef.fogEnd === "number") view.fogEnd = glRef.fogEnd;
  return Object.keys(view).length > 0 ? view : undefined;
}

// --------------------------------------------------------------------------
// Elements + representations
// --------------------------------------------------------------------------

function liftElement(mol: moorhen.Molecule, fileName: string): SceneElement | null {
  const reps = (mol.representations ?? []).filter((r) => r && r.style);
  // Only lift visible reps; hidden ones are usually leftovers the user
  // wouldn't expect in a captured scene.
  const visible = reps.filter((r) => r.visible !== false);
  if (visible.length === 0) return null;

  const out: SceneRepresentation[] = visible.map(liftRepresentation);
  return { file: fileName, representations: out };
}

function liftRepresentation(rep: moorhen.MoleculeRepresentation): SceneRepresentation {
  const out: SceneRepresentation = { style: rep.style };
  if (rep.cid && rep.cid !== "/*/*/*/*") out.selection = rep.cid;

  const colour = liftColour(rep.colourRules ?? []);
  if (colour !== undefined) out.colour = colour;

  return out;
}

// --------------------------------------------------------------------------
// Colours — the only bit that does pattern matching
// --------------------------------------------------------------------------

const NAMED_MULTI_RULES = new Set([
  "b-factor",
  "b-factor-norm",
  "af2-plddt",
  "secondary-structure",
  "jones-rainbow",
  "mol-symm",
]);

function liftColour(rules: moorhen.ColourRule[]): SceneColour | undefined {
  if (rules.length === 0) return undefined;

  // We only attempt to lift the *first* rule. Multiple rules per
  // representation is rare in normal Moorhen UI flows; if we see it,
  // the safer thing is to emit the first one and drop the rest rather
  // than invent some lossy merge. (A future v2 could lift a sequence.)
  const r = rules[0];

  // 1. Named multi-rule scheme (b-factor, etc.).
  if (NAMED_MULTI_RULES.has(r.ruleType)) {
    return r.ruleType as SceneColour;
  }

  // 2. Single-colour rule with a hex colour (Moorhen's "molecule" ruleType
  //    or anything with isMultiColourRule=false).
  if (!r.isMultiColourRule && typeof r.color === "string" && /^#[0-9a-fA-F]{6}$/.test(r.color)) {
    return r.color;
  }

  // 3. by-domain shape: a single multi-rule whose args[0] is a string
  //    of `//chain/start-end^#rrggbb|...` segments. Recognise and emit
  //    as `by-domain` (the calling code can reconstitute domains from
  //    the top-level domains: block; we don't try to lift those here
  //    because we don't know what the segments *mean*).
  if (
    r.isMultiColourRule &&
    typeof r.args?.[0] === "string" &&
    /^\/\/[^/]+\/-?\d+--?\d+\^#[0-9a-fA-F]{6}(\|\/\/[^/]+\/-?\d+--?\d+\^#[0-9a-fA-F]{6})*$/.test(
      r.args[0] as string,
    )
  ) {
    return "by-domain";
  }

  // 4. Escape hatch: keep the rule verbatim. Lossless, ugly.
  return {
    raw: {
      ruleType: r.ruleType,
      args: r.args ?? [],
      isMultiColourRule: r.isMultiColourRule,
      applyColourToNonCarbonAtoms: r.applyColourToNonCarbonAtoms,
    },
  };
}

// --------------------------------------------------------------------------
// Small utilities
// --------------------------------------------------------------------------

function sanitiseName(name?: string): string {
  if (!name) return "";
  // Strip URL-ish punctuation; YAML keys are happier with quiet identifiers.
  return name.replace(/[^A-Za-z0-9_.-]+/g, "_");
}

function toArray3(v: number[] | Float32Array | undefined): [number, number, number] | undefined {
  if (!v || v.length < 3) return undefined;
  return [round(v[0], 3), round(v[1], 3), round(v[2], 3)];
}

function toArray4(
  v: number[] | Float32Array | undefined,
): [number, number, number, number] | undefined {
  if (!v || v.length < 4) return undefined;
  return [round(v[0], 4), round(v[1], 4), round(v[2], 4), round(v[3], 4)];
}

function round(n: number, dp: number): number {
  const f = Math.pow(10, dp);
  return Math.round(n * f) / f;
}

function hasAnyValue(o: Record<string, unknown>): boolean {
  return Object.values(o).some((v) => v !== undefined);
}
