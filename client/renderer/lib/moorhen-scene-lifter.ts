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
  SceneColourSelection,
  SceneDomain,
  SceneElement,
  SceneFileRef,
  SceneMap,
  SceneMapColumns,
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
  /** Optional: Moorhen's monomer library URL root (e.g. "/baby-gru/monomers").
   *  When set, the lifter HEAD-probes `${monomerLibraryPath}/<l>/<UPPER>.cif`
   *  per ligand and skips dicts that already resolve there, so common
   *  cofactors (ATP, HEM, NAD, …) don't bloat the captured YAML.
   *  Without it, only the built-in STANDARD_MONOMERS set is skipped. */
  monomerLibraryPath?: string;
  /** Moorhen maps currently in the store. Each must carry a populated
   *  uniqueId (the MTZ loader URL) plus selectedColumns + isDifference
   *  for the lifter to emit a SceneMap. */
  maps?: moorhen.Map[];
  /** Per-map render state from the mapContourSettings slice, keyed by
   *  map.molNo. The wrapper flattens this from the store before calling
   *  lift — keeps the lifter store-agnostic and easier to unit-test. */
  mapState?: Record<number, MapRenderState>;
  /** molNo of the currently-active map (state.generalStates.activeMap.molNo).
   *  Drives scene.activeMap. */
  activeMapMolNo?: number;
  /** Per-molecule dictionary provenance: molNo → (comp_id → source project
   *  file). When a ligand's dict came from a project file, the lifter emits a
   *  `kind: dictionary` ref by `fileId` (terse, re-fetchable) instead of inlining
   *  cifText. Keyed by molNo (NOT comp_id alone) so two molecules that both call
   *  their ligand LIG, with different chemistry from different dict files, stay
   *  distinct and per-molecule-scoped. */
  dictSources?: Map<number, Map<string, { fileId: number; projectId?: string }>>;
}

/**
 * Per-map render state pulled from Moorhen's redux slices and handed to
 * the lifter in LiftCtx.mapState. All fields optional — the lifter only
 * emits the ones whose captured value differs from Moorhen's defaults
 * so the YAML stays terse.
 */
export interface MapRenderState {
  contourLevel?: number;
  radius?: number;
  alpha?: number;
  style?: "lines" | "solid" | "lit-lines";
  /** Hex (#rrggbb) — for non-difference maps. */
  colour?: string;
  /** Hex — for the positive lobe of a difference map. */
  positiveColour?: string;
  /** Hex — for the negative lobe of a difference map. */
  negativeColour?: string;
  visible?: boolean;
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
    const liftedRefs = liftDictionariesForMolecule(
      mol, molFileName, scene.files!,
      mol.molNo != null ? ctx.dictSources?.get(mol.molNo) : undefined,
      ctx.projectId,
    );
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

  // Maps. Each map gets a SceneFileRef (kind: "mtz") appended to
  // scene.files plus a SceneMap entry referencing it by name. The
  // wrapper flattens map render state from the contour-settings slice
  // into LiftCtx.mapState; activeMap is the molNo whose name we
  // promote to scene.activeMap.
  if (ctx.maps && ctx.maps.length > 0) {
    const maps: SceneMap[] = [];
    ctx.maps.forEach((map, i) => {
      const fileRef = liftMapFileRef(map, ctx, i);
      // De-dupe file refs by name. (Possible if a map and a molecule
      // share a sanitised loader URL — rare but defensive.)
      if (!(scene.files ?? []).some((f) => f.name === fileRef.name)) {
        scene.files = scene.files ?? [];
        scene.files.push(fileRef);
      }
      const mapName = sanitiseName(map.name) || `map${i}`;
      const sceneMap = liftSceneMap(map, fileRef.name, mapName, ctx.mapState);
      maps.push(sceneMap);
      if (ctx.activeMapMolNo !== undefined && map.molNo === ctx.activeMapMolNo) {
        scene.activeMap = mapName;
      }
    });
    if (maps.length > 0) scene.maps = maps;
  }

  // Fold any whole-chain colour lists into domains: + colour: by-domain, so a
  // molecule's per-chain colouring is stated once rather than inline per rep.
  hoistPerChainColours(scene);

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
  sources?: Map<string, { fileId: number; projectId?: string }>,
  fallbackProjectId?: string,
): { name: string; comp_id: string }[] {
  const ligands = mol.ligands ?? [];
  // Dedupe by comp_id — multiple instances of the same ligand share a dict.
  const compIds = Array.from(new Set(ligands.map((l) => l.resName).filter(Boolean)));
  const out: { name: string; comp_id: string }[] = [];

  for (const comp_id of compIds) {
    if (STANDARD_MONOMERS.has(comp_id)) continue;

    // Preferred: the dict came from a project file → emit a terse, re-fetchable
    // `fileId` ref instead of inlining cifText. Deduped by fileId (a multi-comp
    // dict shared by several ligands becomes ONE ref); per-molecule scoping is
    // preserved by listing it in THIS molecule's `dictionaries:`.
    const source = sources?.get(comp_id);
    if (source) {
      const name = `dict-file-${source.fileId}`;
      if (!allFiles.some((f) => f.name === name)) {
        allFiles.push({
          name,
          kind: "dictionary",
          fileId: source.fileId,
          projectId: source.projectId ?? fallbackProjectId,
        });
      }
      out.push({ name, comp_id });
      continue;
    }

    // No project source → inline cifText. dropLibraryDicts() may later remove
    // it if the comp_id resolves in the receiver's monomer library.
    let dictText = "";
    try {
      dictText = mol.getDict(comp_id) || "";
    } catch {
      // getDict may throw for comp_ids without a stored dict — skip them.
      continue;
    }
    if (!dictText) continue;

    const name = `dict-${molFileName}-${comp_id}`;
    if (allFiles.some((f) => f.name === name)) continue;
    allFiles.push({ name, kind: "dictionary", cifText: dictText });
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
 * STRAIGHT lift for Capture: real provenance + hints, NO bundling — but it DOES
 * drop dict refs the receiver's monomer library already has, so library monomers
 * (CL, ATP, …) aren't carried; the receiving Moorhen resolves those itself. Custom
 * ligand dicts stay — as terse `fileId` refs when the project supplied them, or
 * inline `cifText` otherwise. Async because the library check is a HEAD probe.
 */
export async function liftSceneStraight(ctx: LiftCtx): Promise<{
  scene: MoorhenScene;
  hints: SceneLiftHints;
}> {
  const { scene, hints } = liftSceneWithHints(ctx);
  if (ctx.monomerLibraryPath) {
    await dropLibraryDicts(scene, ctx.molecules, ctx.monomerLibraryPath);
  }
  return { scene, hints };
}

/**
 * Remove INLINE (cifText) dictionary refs whose comp_id resolves in the
 * receiver's monomer library — the receiving Moorhen fetches the same dict from
 * the same library when it loads the coords, so carrying our copy only bloats the
 * YAML (and risks a stale variant). `fileId`/`url`/`bundle` dict refs are NEVER
 * dropped: a project explicitly supplied that dict, so it must travel and win
 * even if its comp_id name collides with a library entry.
 */
async function dropLibraryDicts(
  scene: MoorhenScene,
  molecules: moorhen.Molecule[],
  monomerLibraryPath: string,
): Promise<void> {
  const libraryHas = await probeMonomerLibrary(molecules, monomerLibraryPath);
  if (libraryHas.size === 0) return;
  const dropNames = new Set<string>();
  for (const ref of scene.files ?? []) {
    if (ref.kind !== "dictionary" || !ref.cifText) continue; // only inline dicts
    // Names are `dict-<molFileName>-<COMP_ID>`; take the trailing comp_id (the
    // molFileName may itself contain hyphens).
    const m = /^dict-.*-([A-Za-z0-9]+)$/.exec(ref.name);
    if (m && libraryHas.has(m[1].toUpperCase())) dropNames.add(ref.name);
  }
  if (dropNames.size === 0) return;
  scene.files = (scene.files ?? []).filter((f) => !dropNames.has(f.name));
  for (const el of scene.elements ?? []) {
    if (!el.dictionaries) continue;
    el.dictionaries = el.dictionaries.filter((n) => !dropNames.has(n));
    if (el.dictionaries.length === 0) delete el.dictionaries;
  }
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

  // 0. Drop inline dict refs the receiver's monomer library already has (shared
  //    with the straight-lift path).
  if (ctx.monomerLibraryPath) {
    await dropLibraryDicts(scene, ctx.molecules, ctx.monomerLibraryPath);
  }

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
      // Replace the unresolvable relativeUrl: form with the portable bundle: form.
      ref.bundle = assetPath;
      delete ref.relativeUrl;
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

/**
 * Resolve a SceneFileRef to a URL the browser can `fetch()`, or null
 * if the ref doesn't carry enough information (e.g. a bare `path:` or
 * an already-bundled ref). Supplied by the wrapper, which knows the
 * site-specific proxy routes.
 */
export type SceneRefUrlResolver = (ref: SceneFileRef) => string | null;

export interface PromoteCtx {
  /** Parsed scene to promote. Mutated in place; pass a clone if you
   *  need the original. */
  scene: MoorhenScene;
  /** Asset bytes already present (e.g. from a .scene.zip the user opened).
   *  New assets are merged into a fresh map so the caller can choose to
   *  keep, discard, or diff. */
  existingAssets: Map<string, ArrayBuffer>;
  /** Function that yields a fetchable URL for any URL-shaped ref kind
   *  (pdb / url / fileId+projectId / path-as-URL). */
  resolveUrl: SceneRefUrlResolver;
  /** Live molecules — used to re-collect library-resolvable dicts that
   *  the lifter omitted from the YAML. Optional: if not provided, no
   *  library-dict re-collection happens. */
  molecules?: moorhen.Molecule[];
  /** Optional monomer library URL root for library-dict re-collection. */
  monomerLibraryPath?: string;
}

/**
 * Walk a scene and gather every external dependency into a single
 * self-contained bundle. For each URL-resolvable ref (pdb / url /
 * fileId+projectId / path), fetch the bytes and rewrite the ref to
 * `bundle: assets/<unique>.<ext>`. For each ligand whose dict was
 * skipped at lift time (because it lives in the monomer library),
 * fetch it from the library and add a fresh dict ref.
 *
 * Asset paths are name-mangled on collision (`name.cif`, `name-2.cif`).
 *
 * Returns the in-place-mutated scene and a fresh assets map (merging
 * `existingAssets` with all newly-fetched bytes). The caller is
 * responsible for serialising the scene to YAML and zipping the assets.
 */
export async function promoteSceneToPortable(ctx: PromoteCtx): Promise<{
  scene: MoorhenScene;
  assets: Map<string, ArrayBuffer>;
  warnings: string[];
}> {
  const { scene, existingAssets, resolveUrl, molecules, monomerLibraryPath } = ctx;
  const assets = new Map(existingAssets);
  const warnings: string[] = [];

  // Names already taken in the asset map — used to mangle on collision.
  const usedPaths = new Set<string>(assets.keys());
  const claim = (preferred: string): string => {
    if (!usedPaths.has(preferred)) {
      usedPaths.add(preferred);
      return preferred;
    }
    const dot = preferred.lastIndexOf(".");
    const base = dot > 0 ? preferred.slice(0, dot) : preferred;
    const ext = dot > 0 ? preferred.slice(dot) : "";
    for (let i = 2; ; i++) {
      const candidate = `${base}-${i}${ext}`;
      if (!usedPaths.has(candidate)) {
        usedPaths.add(candidate);
        return candidate;
      }
    }
  };

  // 1. URL-resolvable refs → bundle. Coord and dict refs alike: anything
  //    with a fetchable URL gets pulled and stashed under assets/.
  for (const ref of scene.files ?? []) {
    if (ref.bundle) continue; // already bundled
    if (ref.kind === "dictionary" && ref.cifText) {
      // Inline cifText: move it to a bundle asset (mirrors liftSceneToBundle
      // step 2). Cheaper than refetching and works without a URL.
      const assetPath = claim(`assets/${ref.name}.cif`);
      assets.set(assetPath, new TextEncoder().encode(ref.cifText).buffer);
      ref.bundle = assetPath;
      delete ref.cifText;
      continue;
    }
    const url = resolveUrl(ref);
    if (!url) {
      warnings.push(`No URL for ref "${ref.name}"; left as-is`);
      continue;
    }
    try {
      const res = await fetch(url);
      if (!res.ok) {
        warnings.push(`Fetch ${url} for "${ref.name}" returned ${res.status}; left as-is`);
        continue;
      }
      const buf = await res.arrayBuffer();
      const ext = extForRef(ref, url);
      const assetPath = claim(`assets/${ref.name}.${ext}`);
      assets.set(assetPath, buf);
      // Replace symbolic source fields with the bundle reference. Keep
      // `kind` since downstream resolver cares.
      delete ref.pdb;
      delete ref.url;
      delete ref.relativeUrl;
      delete ref.fileId;
      delete ref.projectId;
      delete ref.projectName;
      ref.bundle = assetPath;
    } catch (err) {
      warnings.push(
        `Fetch ${url} for "${ref.name}" failed (${err instanceof Error ? err.message : "unknown"}); left as-is`,
      );
    }
  }

  // 2. Library-dict re-collection. The lifter strips dicts whose
  //    comp_id resolves to the monomer library; self-contained mode
  //    wants them back. Walk live molecules; for each non-standard
  //    comp_id that the library has AND isn't already in the scene's
  //    dict refs, fetch it and append a fresh ref.
  if (molecules && monomerLibraryPath) {
    const existingDictNames = new Set(
      (scene.files ?? []).filter((f) => f.kind === "dictionary").map((f) => f.name),
    );
    for (const mol of molecules) {
      const molFileName = findFileNameForMol(scene, mol);
      if (!molFileName) continue;
      const ligands = mol.ligands ?? [];
      const compIds = Array.from(new Set(ligands.map((l) => l.resName).filter(Boolean)));
      for (const compId of compIds) {
        if (STANDARD_MONOMERS.has(compId)) continue; // never travelled
        const dictName = `dict-${molFileName}-${compId}`;
        if (existingDictNames.has(dictName)) continue;
        const url = `${monomerLibraryPath}/${compId[0].toLowerCase()}/${compId.toUpperCase()}.cif`;
        try {
          const res = await fetch(url);
          if (!res.ok) continue;
          const buf = await res.arrayBuffer();
          const assetPath = claim(`assets/${dictName}.cif`);
          assets.set(assetPath, buf);
          const newRef: SceneFileRef = {
            name: dictName,
            kind: "dictionary",
            bundle: assetPath,
          };
          (scene.files ??= []).push(newRef);
          existingDictNames.add(dictName);
          const el = (scene.elements ?? []).find((e) => e.file === molFileName);
          if (el) {
            el.dictionaries = el.dictionaries ?? [];
            if (!el.dictionaries.includes(dictName)) el.dictionaries.push(dictName);
          }
        } catch {
          // network errors → leave the dict out; receiver's own Moorhen
          // will fetch from its own library at apply time.
        }
      }
    }
  }

  return { scene, assets, warnings };
}

/**
 * Best-effort extension picker for a fetched asset. Dict refs are CIF;
 * coord refs fall back to the URL's trailing extension, else "cif".
 */
function extForRef(ref: SceneFileRef, url: string): string {
  if (ref.kind === "dictionary") return "cif";
  if (ref.kind === "mtz") return "mtz";
  const m = /\.([A-Za-z0-9]{2,4})(?:\?|$)/.exec(url);
  if (m) return m[1].toLowerCase();
  return "cif";
}

/**
 * Find which scene.files entry corresponds to a loaded molecule, using
 * the same heuristics as liftScene: sanitised name or mol-index. Used
 * by promote to attach re-collected dict refs to the right element.
 */
function findFileNameForMol(
  scene: MoorhenScene,
  mol: moorhen.Molecule,
): string | null {
  const sanitised = sanitiseName(mol.name);
  for (const f of scene.files ?? []) {
    if (f.name === sanitised) return f.name;
  }
  // Project-internal: match by fileId in uniqueId.
  const fid = extractFileIdFromUniqueId(mol.uniqueId ?? "");
  if (fid !== null) {
    for (const f of scene.files ?? []) {
      if (f.fileId === fid) return f.name;
    }
  }
  // Bundle round-trip: uniqueId is `bundle:<path>`.
  if (mol.uniqueId?.startsWith("bundle:")) {
    const path = mol.uniqueId.slice("bundle:".length);
    for (const f of scene.files ?? []) {
      if (f.bundle === path) return f.name;
    }
  }
  return null;
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

  // Fallback: a uniqueId that isn't an absolute URL is an origin-relative
  // loader URL (e.g. /api/proxy/pdbe/…) — emit it as relativeUrl.
  if (mol.uniqueId) {
    return { name, relativeUrl: mol.uniqueId };
  }

  // Last resort: a placeholder url that the user will edit.
  return { name, url: "TODO: source URL for this molecule" };
}

// --------------------------------------------------------------------------
// Maps
// --------------------------------------------------------------------------

/** True if a loaded MoorhenMap came from a real-space CCP4 map file (incl. a
 *  mask), tagged at load by the wrapper map loaders, rather than from MTZ. */
function isCcp4MapFile(map: moorhen.Map): boolean {
  return !!(map as unknown as { isCcp4MapFile?: boolean }).isCcp4MapFile;
}

/**
 * Lift a SceneFileRef for a loaded map. MTZ maps get kind: "mtz"; real-space
 * CCP4 map files (incl. masks) get kind: "map". The name is disambiguated
 * from any coord ref by a kind suffix.
 */
function liftMapFileRef(map: moorhen.Map, ctx: LiftCtx, index: number): SceneFileRef {
  const kind: SceneFileRef["kind"] = isCcp4MapFile(map) ? "map" : "mtz";
  const baseName = sanitiseName(map.name) || `map${index}`;
  const name = `${baseName}__${kind}`;

  const fileId = extractFileIdFromUniqueId(map.uniqueId ?? "");
  if (fileId !== null && ctx.projectId) {
    return {
      name,
      kind,
      projectId: ctx.projectId,
      projectName: ctx.projectName,
      fileId,
    };
  }
  if (map.uniqueId && /^https?:\/\//.test(map.uniqueId)) {
    return { name, kind, url: map.uniqueId };
  }
  if (map.uniqueId) {
    return { name, kind, relativeUrl: map.uniqueId };
  }
  return { name, kind, url: "TODO: URL for this map" };
}

/**
 * Build a SceneMap entry for a loaded MoorhenMap. Pulls the column
 * spec straight from `map.selectedColumns` and reads optional render
 * state from the wrapper-supplied `mapState` keyed by molNo. Difference
 * maps only emit positive/negative colours; non-difference maps emit
 * a single `colour` field.
 */
function liftSceneMap(
  map: moorhen.Map,
  fileName: string,
  mapName: string,
  mapState?: Record<number, MapRenderState>,
): SceneMap {
  // Real-space CCP4 map files (incl. masks) carry no columns; MTZ maps do.
  const mapFile = isCcp4MapFile(map);
  const out: SceneMap = mapFile
    ? { name: mapName, file: fileName }
    : { name: mapName, file: fileName, columns: stripUndefinedColumns(map.selectedColumns) };
  if (mapFile && (map as unknown as { isCcp4Mask?: boolean }).isCcp4Mask) {
    out.isMask = true;
  }
  if (map.isDifference) out.isDifference = true;

  const s = mapState?.[map.molNo];
  if (!s) return out;

  if (s.contourLevel !== undefined) out.contourLevel = round(s.contourLevel, 4);
  if (s.radius !== undefined) out.radius = round(s.radius, 2);
  if (s.alpha !== undefined && s.alpha !== 1) out.alpha = round(s.alpha, 3);
  if (s.style && s.style !== "lit-lines") out.style = s.style;
  if (s.visible === false) out.visible = false;

  // Colour fields: difference maps split into positive/negative;
  // non-difference maps carry a single colour.
  if (map.isDifference) {
    if (s.positiveColour) out.positiveColour = s.positiveColour;
    if (s.negativeColour) out.negativeColour = s.negativeColour;
  } else if (s.colour) {
    out.colour = s.colour;
  }
  return out;
}

/**
 * Drop undefined / empty-string fields from a column spec so the YAML
 * doesn't carry no-op entries like `Fobs:` (with no value).
 */
function stripUndefinedColumns(c: Partial<SceneMapColumns> | undefined): SceneMapColumns {
  const out: SceneMapColumns = {};
  if (!c) return out;
  const strKeys: (keyof SceneMapColumns)[] = ["F", "PHI", "Fobs", "SigFobs", "FreeR"];
  for (const k of strKeys) {
    const v = c[k];
    if (typeof v === "string" && v.length > 0) (out as Record<string, unknown>)[k] = v;
  }
  if (c.useWeight) out.useWeight = true;
  if (c.calcStructFact) out.calcStructFact = true;
  return out;
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

  // Per-representation opacity (Moorhen `nonCustomOpacity`); only emit a
  // non-default value (default 1 = opaque ⇒ omit, matching the map `alpha`).
  const a = (rep as { nonCustomOpacity?: number }).nonCustomOpacity;
  if (typeof a === "number" && a < 1) out.alpha = round(a, 3);

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

const HEX_RE = /^#[0-9a-fA-F]{6}([0-9a-fA-F]{2})?$/;

/** A single-colour rule with a hex colour (Moorhen "molecule" ruleType etc.). */
function isSingleHexRule(r: moorhen.ColourRule): boolean {
  return !r.isMultiColourRule && typeof r.color === "string" && HEX_RE.test(r.color);
}

function liftColour(rules: moorhen.ColourRule[]): SceneColour | undefined {
  if (rules.length === 0) return undefined;

  // All single-colour rules → a single hex (one rule) or a per-selection list
  // (many). The list is the faithful capture of coot's default per-chain
  // colouring: one entry per chain with the colour coot assigned — instead of
  // mistaking the first chain's hex for the whole representation. Whole-chain
  // and residue-range CIDs are the same shape, so by-domain re-applies through
  // the same path.
  if (rules.every(isSingleHexRule)) {
    if (rules.length === 1) return rules[0].color;
    return rules.map((r) => ({ selection: r.cid, colour: r.color }));
  }

  const r = rules[0];

  // 1. Named multi-rule scheme (b-factor, etc.).
  if (NAMED_MULTI_RULES.has(r.ruleType)) {
    return r.ruleType as SceneColour;
  }

  // 2. A compiled WHOLE-CHAIN per-selection rule: one multi-rule whose args[0]
  //    is `//chain^#hex|...` (no residue range). This is what applying a
  //    hoisted per-chain colouring (domains: + by-domain) produces, so decompose
  //    it back to the list — hoistPerChainColours then re-folds it into
  //    domains: + by-domain, keeping capture→apply→capture stable.
  if (
    r.isMultiColourRule &&
    typeof r.args?.[0] === "string" &&
    /^\/\/[^/^|]+\^#[0-9a-fA-F]{6}(?:[0-9a-fA-F]{2})?(\|\/\/[^/^|]+\^#[0-9a-fA-F]{6}(?:[0-9a-fA-F]{2})?)*$/.test(
      r.args[0] as string,
    )
  ) {
    return (r.args[0] as string).split("|").map((seg) => {
      const ix = seg.lastIndexOf("^");
      return { selection: seg.slice(0, ix), colour: seg.slice(ix + 1) };
    });
  }

  // 3. by-domain shape: a single multi-rule whose args[0] is a string
  //    of `//chain/start-end^#rrggbb|...` segments. Recognise and emit
  //    as `by-domain` (the calling code can reconstitute domains from
  //    the top-level domains: block; we don't try to lift those here
  //    because we don't know what the segments *mean*).
  if (
    r.isMultiColourRule &&
    typeof r.args?.[0] === "string" &&
    /^\/\/[^/]+\/-?\d+--?\d+\^#[0-9a-fA-F]{6}(?:[0-9a-fA-F]{2})?(\|\/\/[^/]+\/-?\d+--?\d+\^#[0-9a-fA-F]{6}(?:[0-9a-fA-F]{2})?)*$/.test(
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

/** Chain id of a whole-chain CID ("//A", "/1/A"), else null (has residue/atom
 *  parts, or is a wildcard). */
function chainOfWholeChainCid(cid: string): string | null {
  const parts = cid.split("/"); // "//A" -> ["","","A"]; "/1/A" -> ["","1","A"]
  if (parts.length === 3 && parts[2] !== "" && parts[2] !== "*") return parts[2];
  return null;
}

/**
 * Hoist whole-chain colour lists into the shared `domains:` block + `colour:
 * by-domain` on the representations, so a molecule's per-chain colouring is
 * stated ONCE instead of repeated inline on every representation — the same
 * "define a colour map once, representations adopt it" pattern by-domain already
 * uses, at whole-chain granularity.
 *
 * Conservative: only fires when every adopting rep's colour list is entirely
 * whole-chain CIDs, the chain→colour mapping is consistent across reps, and the
 * `domains:` block is empty (so authored domains are never disturbed). Anything
 * else (residue-range lists, conflicting per-chain colours) is left inline.
 */
function hoistPerChainColours(scene: MoorhenScene): void {
  if (!scene.elements || (scene.domains && scene.domains.length > 0)) return;
  const chainColour = new Map<string, string>();
  const adopters: SceneRepresentation[] = [];
  let conflict = false;
  for (const el of scene.elements) {
    for (const rep of el.representations ?? []) {
      const c = rep.colour;
      if (!Array.isArray(c) || c.length === 0) continue;
      const chains = c.map((e) => chainOfWholeChainCid(e.selection));
      if (chains.some((ch) => ch === null)) continue; // not all whole-chain
      adopters.push(rep);
      c.forEach((e, i) => {
        const ch = chains[i] as string;
        const prev = chainColour.get(ch);
        if (prev !== undefined && prev !== e.colour) conflict = true;
        chainColour.set(ch, e.colour);
      });
    }
  }
  if (conflict || chainColour.size === 0 || adopters.length === 0) return;
  scene.domains = [...chainColour.entries()].map(([chain, color]) => ({
    name: chain,
    selection: `//${chain}`, // whole chain, CID form
    color,
  }));
  for (const rep of adopters) rep.colour = "by-domain";
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

/**
 * HEAD-probe each candidate comp_id against the monomer library and
 * return the set that resolves there (uppercase). A non-2xx (or thrown
 * fetch) is treated as "not present" — conservative: we'd rather ship
 * a dict the receiver might already have than drop one they don't.
 *
 * `STANDARD_MONOMERS` are short-circuited (always treated as present)
 * so we don't issue probes for the 20 amino acids on every capture.
 */
async function probeMonomerLibrary(
  molecules: moorhen.Molecule[],
  monomerLibraryPath: string,
): Promise<Set<string>> {
  const compIds = new Set<string>();
  for (const mol of molecules) {
    for (const lig of mol.ligands ?? []) {
      if (lig.resName) compIds.add(lig.resName.toUpperCase());
    }
  }
  const present = new Set<string>();
  await Promise.all(
    [...compIds].map(async (cid) => {
      if (STANDARD_MONOMERS.has(cid)) {
        present.add(cid);
        return;
      }
      const url = `${monomerLibraryPath}/${cid[0].toLowerCase()}/${cid}.cif`;
      try {
        const res = await fetch(url, { method: "HEAD" });
        if (res.ok) present.add(cid);
      } catch {
        // network errors → treat as not-present; the dict will travel.
      }
    }),
  );
  return present;
}
