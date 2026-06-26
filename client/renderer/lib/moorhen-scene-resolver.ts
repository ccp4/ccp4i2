/**
 * Moorhen Scene resolver.
 *
 * Takes a parsed MoorhenScene plus the molecules currently loaded in
 * Moorhen and applies the scene: clears existing representations on each
 * targeted molecule, adds new ones per the scene, with colour rules
 * compiled from the scene's domain definitions.
 *
 * v0 scope: only the simplest path. Files in the scene must match
 * structures already loaded in Moorhen (matched by file id extracted
 * from the molecule's uniqueId). Fetching new structures from `files:`
 * entries (url/path/job+param) is deferred — drop the file on Moorhen
 * first, then apply the scene.
 *
 * Camera (view block) is applied by dispatching the same Redux actions
 * the existing PasteViewLinkField uses.
 */

import type { Dispatch } from "redux";
import {
  setOrigin,
  setQuat,
  setZoom,
  setClipStart,
  setClipEnd,
  setFogStart,
  setFogEnd,
  setResetClippingFogging,
  setBackgroundColor,
  setRequestDrawScene,
  addCustomRepresentation,
  removeCustomRepresentation,
  setContourLevel,
  setMapRadius,
  setMapAlpha,
  setMapStyle,
  setMapColours,
  setPositiveMapColours,
  setNegativeMapColours,
  showMap,
  hideMap,
  setActiveMap,
} from "moorhen/react-lib";
import type { moorhen } from "moorhen/types/moorhen";

// We deliberately don't import MoorhenColourRule from "moorhen/react-lib". Although
// the class is exported in Moorhen's source, the installed package's
// bundled moorhen.js does not re-expose it on the package export object
// (runtime "is not a constructor"). Instead, we add colour rules via
// MoleculeRepresentation.addColourRule, which constructs the rule
// internally with the right commandCentre + parent molecule wiring.

import { extractFileIdFromUniqueId } from "./moorhen-view-state";
import {
  MoorhenScene,
  SceneColour,
  SceneColourSelection,
  SceneDomain,
  SceneFileRef,
  SceneMap,
  SceneRepresentation,
  SceneLsqMatch,
  SceneSuperpose,
  isSceneHexColour,
  isSceneNamedColour,
  isSceneRawColour,
} from "../types/moorhen-scene";

// --------------------------------------------------------------------------
// Result + log entries
// --------------------------------------------------------------------------

/** One modification noted by the clamp-and-log resolver policy. */
export interface SceneResolveLogEntry {
  file: string;            // SceneFileRef.name
  domain: string;          // SceneDomain.name
  message: string;         // human-readable
}

export interface SceneResolveResult {
  /** Files in the scene that couldn't be matched to a loaded molecule. */
  unresolvedFiles: string[];
  /** Modifications made by clamp-and-log (or anything else worth noting). */
  log: SceneResolveLogEntry[];
  /** Per-file representation counts actually applied. */
  applied: { file: string; molNo: number; representations: number }[];
}

// --------------------------------------------------------------------------
// Entry point
// --------------------------------------------------------------------------

/**
 * Loads a structure described by a SceneFileRef and returns the
 * MoorhenMolecule that ends up in the store, or null if the ref is not
 * fetchable on its own (e.g. a bare `path:` with no matching loaded mol).
 *
 * The resolver doesn't know how to talk to Coot or the store directly —
 * it just calls this. The wrapper provides the implementation in terms
 * of its existing fetchMolecule / fetchFile / apiText machinery.
 */
export type SceneFileFetcher = (
  ref: SceneFileRef,
) => Promise<moorhen.Molecule | null>;

/**
 * Fetches the raw text of a dictionary CIF and returns it. The resolver
 * is then responsible for loading it into Coot (globally and/or scoped
 * to specific molecule molNos). A single dictionary file may contain
 * multiple `data_comp_*` blocks — Coot loads them all in one call, so
 * the resolver doesn't need to split.
 */
export type SceneDictionaryFetcher = (
  ref: SceneFileRef,
) => Promise<string | null>;

/**
 * Loads a dictionary CIF string into Coot's dictionary store for a
 * specific molecule (or globally, when molNo is -999999). Mirrors the
 * existing wrapper pattern: load globally first so coords containing the
 * monomers parse, then re-associate per-molecule so structures using
 * differently-defined ligands with the same residue name don't collide.
 */
export type SceneDictionaryLoader = (
  dictText: string,
  molNo: number,
) => Promise<void>;

/**
 * Fetch an MTZ file ref's bytes and load it into Coot as a Moorhen map
 * with the given column spec. Owned by the wrapper because it
 * instantiates MoorhenMap and dispatches addMap — neither of which
 * belong in the resolver.
 *
 * Returns the resulting moorhen.Map (so the resolver can pick up its
 * molNo for redux state updates), or null on fetch / parse failure.
 */
export type SceneMapFetcher = (
  ref: SceneFileRef,
  sceneMap: SceneMap,
) => Promise<moorhen.Map | null>;

interface ResolveCtx {
  scene: MoorhenScene;
  molecules: moorhen.Molecule[];
  /** Maps currently loaded in Moorhen. Used to match scene.maps entries
   *  to existing maps before falling back to mapFetcher. */
  maps?: moorhen.Map[];
  dispatch: Dispatch;
  /** Optional. When provided, the resolver fetches structures that aren't
   *  already loaded but whose file ref carries enough info (pdb, url,
   *  fileId+projectId). Without it, unmatched files just go to
   *  unresolvedFiles. */
  fetcher?: SceneFileFetcher;
  /** Optional. When provided alongside `dictionaryLoader`, dict file
   *  refs (`kind: dictionary`) are fetched into raw CIF text. */
  dictionaryFetcher?: SceneDictionaryFetcher;
  /** Optional. Required for any dictionary scoping to take effect. The
   *  resolver calls this once per (dict, molecule) pair. */
  dictionaryLoader?: SceneDictionaryLoader;
  /** Optional. Required for scene.maps[] to be applied. Without it,
   *  map entries are dropped with a log entry. */
  mapFetcher?: SceneMapFetcher;
  /** Live glRef snapshot (zoom, fogClipOffset). Used by `view.clip:
   *  { front, back }` to derive clip/fog the way coot does (clip = zoom*depth,
   *  fog offset by fogClipOffset). Falls back to sane defaults if absent. */
  glRef?: { zoom: number; fogClipOffset: number };
}

// --------------------------------------------------------------------------
// Geometry: bounding sphere of a selection (for view.slab)
// --------------------------------------------------------------------------

interface GemmiVec<T> { size(): number; get(i: number): T; delete(): void; }
interface GemmiAtom { pos: { x: number; y: number; z: number }; }
interface CCP4GemmiModule {
  Selection: new (cid: string) => { delete: () => void };
  selection_get_models(sel: unknown, struct: unknown): GemmiVec<unknown>;
  selection_get_chains(sel: unknown, model: unknown): GemmiVec<unknown>;
  selection_get_residues(sel: unknown, chain: unknown): GemmiVec<unknown>;
  selection_get_atoms(sel: unknown, residue: unknown): GemmiVec<GemmiAtom>;
}

/**
 * Centroid + bounding radius (Å) of a CID selection, walked over the molecule's
 * cached gemmi structure via the gemmi WASM module Moorhen exposes as
 * `window.CCP4Module`. Orientation-independent. Returns null if gemmi / the
 * structure is unavailable or the selection matched no atoms. emscripten objects
 * are .delete()d so the WASM heap doesn't grow on repeated applies.
 */
function selectionBoundingSphere(
  mol: moorhen.Molecule,
  cid: string,
): { centre: [number, number, number]; radius: number } | null {
  const M =
    typeof window !== "undefined"
      ? (window as unknown as { CCP4Module?: CCP4GemmiModule }).CCP4Module
      : undefined;
  const struct = (mol as unknown as { gemmiStructure?: unknown }).gemmiStructure;
  if (!M || !struct || typeof M.Selection !== "function") return null;

  const del = (x: unknown) => {
    try { (x as { delete?: () => void } | null)?.delete?.(); } catch { /* ignore */ }
  };
  const xs: number[] = [], ys: number[] = [], zs: number[] = [];
  // Accept "||"-separated CIDs (the representation-path union operator) as well
  // as primitive gemmi CIDs, so one selection syntax works everywhere — a scene
  // author shouldn't have to know that view directives feed gemmi directly while
  // representations split on "||". Walk each sub-selection and union its atoms; a
  // malformed part is skipped, not fatal.
  const parts = cid.split("||").map((s) => s.trim()).filter(Boolean);
  for (const part of parts.length ? parts : [cid]) {
    let sel: { delete: () => void } | null = null;
    try {
      sel = new M.Selection(part);
      const models = M.selection_get_models(sel, struct);
      for (let i = 0; i < models.size(); i++) {
        const model = models.get(i);
        const chains = M.selection_get_chains(sel, model);
        for (let j = 0; j < chains.size(); j++) {
          const chain = chains.get(j);
          const residues = M.selection_get_residues(sel, chain);
          for (let k = 0; k < residues.size(); k++) {
            const res = residues.get(k);
            const atoms = M.selection_get_atoms(sel, res);
            for (let l = 0; l < atoms.size(); l++) {
              const a = atoms.get(l);
              xs.push(a.pos.x); ys.push(a.pos.y); zs.push(a.pos.z);
              del(a);
            }
            atoms.delete(); del(res);
          }
          residues.delete(); del(chain);
        }
        chains.delete(); del(model);
      }
      models.delete();
    } catch {
      // skip this sub-selection; other parts may still contribute atoms
    } finally {
      del(sel);
    }
  }

  const n = xs.length;
  if (n === 0) return null;
  const cx = xs.reduce((s, v) => s + v, 0) / n;
  const cy = ys.reduce((s, v) => s + v, 0) / n;
  const cz = zs.reduce((s, v) => s + v, 0) / n;
  let r = 0;
  for (let i = 0; i < n; i++) {
    const d = Math.hypot(xs[i] - cx, ys[i] - cy, zs[i] - cz);
    if (d > r) r = d;
  }
  return { centre: [cx, cy, cz], radius: r };
}

/**
 * Resolve the molecule a view directive (centre/slab) acts on. With an explicit
 * `file`, look it up. Without one, default to the sole loaded molecule — the
 * common "only one structure open" case — and report an error otherwise (none
 * loaded, or ambiguous with several). `ref` is a display string for logging.
 */
function resolveViewMolecule(
  file: string | undefined,
  fileBindings: Map<string, moorhen.Molecule>,
): { mol?: moorhen.Molecule; ref: string; error?: string } {
  if (file !== undefined) {
    const mol = fileBindings.get(file);
    return mol
      ? { mol, ref: file }
      : { ref: file, error: `file "${file}" not bound to a loaded molecule` };
  }
  if (fileBindings.size === 1) {
    return { mol: [...fileBindings.values()][0], ref: "(sole molecule)" };
  }
  return {
    ref: "(unspecified)",
    error:
      fileBindings.size === 0
        ? "no file specified and no molecule loaded"
        : `no file specified but ${fileBindings.size} molecules loaded — add a file:`,
  };
}

/**
 * Apply a scene to the currently-loaded molecules.
 *
 * Returns a SceneResolveResult describing what was applied and what
 * couldn't be resolved. Throws only on programmer errors; missing
 * structures and domain clamping are reported in the result, not raised.
 */
export async function applyScene(ctx: ResolveCtx): Promise<SceneResolveResult> {
  const { scene, molecules, dispatch, fetcher, dictionaryFetcher, dictionaryLoader, mapFetcher } = ctx;
  const maps = ctx.maps ?? [];
  const result: SceneResolveResult = {
    unresolvedFiles: [],
    log: [],
    applied: [],
  };

  const policy = scene.resolver?.onMissingResidues ?? "clamp-and-log";

  // Split file refs by kind. Dictionaries are fetched and loaded
  // GLOBALLY first, so coords containing those monomers parse correctly
  // when they're loaded immediately after. Per-element re-association
  // (the scoping step) happens later, once coord molNos are known.
  const allFiles = scene.files ?? [];
  const dictRefs = allFiles.filter((f) => f.kind === "dictionary");
  const coordRefs = allFiles.filter((f) => f.kind !== "dictionary");

  // 1a. Fetch and globally-load dictionary text. Keep the raw text
  //     keyed by name so we can re-load per-molecule later.
  const dictTexts = new Map<string, string>();
  for (const fr of dictRefs) {
    if (!dictionaryFetcher || !dictionaryLoader) {
      result.unresolvedFiles.push(fr.name);
      continue;
    }
    if (!isFetchable(fr)) {
      result.unresolvedFiles.push(fr.name);
      continue;
    }
    try {
      const text = await dictionaryFetcher(fr);
      if (!text) {
        result.unresolvedFiles.push(fr.name);
        continue;
      }
      dictTexts.set(fr.name, text);
      // Global association: imol = -999999 (Coot's "any" sentinel).
      // Coots parses every `data_comp_*` block in a single call, so
      // multi-comp dicts are handled in one shot.
      await dictionaryLoader(text, -999999);
    } catch (e) {
      console.warn(`[scene] dictionary fetch failed for ${fr.name}:`, e);
      result.unresolvedFiles.push(fr.name);
    }
  }

  // 1b. Match coord scene files to loaded molecules. For any file that
  //     doesn't match, fall back to the fetcher (if provided + ref is
  //     fetchable). Newly-fetched molecules are accumulated into a
  //     local list so later files can also bind against them.
  const livePool: moorhen.Molecule[] = [...molecules];
  const fileBindings = new Map<string, moorhen.Molecule>();
  for (const fr of coordRefs) {
    const existing = matchOneFile(fr, livePool);
    if (existing) {
      fileBindings.set(fr.name, existing);
      continue;
    }
    if (fetcher && isFetchable(fr)) {
      try {
        const fetched = await fetcher(fr);
        if (fetched) {
          livePool.push(fetched);
          fileBindings.set(fr.name, fetched);
          continue;
        }
      } catch (e) {
        console.warn(`[scene] fetch failed for ${fr.name}:`, e);
      }
    }
    result.unresolvedFiles.push(fr.name);
  }

  // 1c. Scope dictionaries to specific molecules. Re-association overrides
  //     the earlier global load for the named molecule, so two molecules
  //     with same-named ligands can carry different chemistry. Order:
  //     globalDictionaries first (apply to every coord molecule), then
  //     per-element dictionaries (which take precedence on their molecule
  //     because they're loaded last).
  if (dictionaryLoader) {
    const globalDictNames = scene.globalDictionaries ?? [];
    for (const dictName of globalDictNames) {
      const text = dictTexts.get(dictName);
      if (!text) continue; // unresolved noted above
      for (const mol of fileBindings.values()) {
        if (mol.molNo === undefined || mol.molNo === null) continue;
        try {
          await dictionaryLoader(text, mol.molNo);
        } catch (e) {
          console.warn(`[scene] global dict ${dictName} failed on molNo ${mol.molNo}:`, e);
        }
      }
    }
    const moleculesNeedingRedraw = new Set<moorhen.Molecule>();
    for (const element of scene.elements ?? []) {
      const elDicts = element.dictionaries ?? [];
      if (elDicts.length === 0) continue;
      const mol = fileBindings.get(element.file);
      if (!mol || mol.molNo === undefined || mol.molNo === null) continue;
      for (const dictName of elDicts) {
        const text = dictTexts.get(dictName);
        if (!text) continue;
        try {
          await dictionaryLoader(text, mol.molNo);
          moleculesNeedingRedraw.add(mol);
        } catch (e) {
          console.warn(`[scene] dict ${dictName} failed on ${element.file}:`, e);
        }
      }
    }
    // Redraw molecules that got per-element dict re-associations, so
    // bond rendering picks up the new geometry before reps are applied.
    for (const mol of moleculesNeedingRedraw) {
      try {
        mol.setAtomsDirty(true);
        await mol.redraw();
      } catch (e) {
        console.warn(`[scene] redraw after dict scope failed on ${mol.name}:`, e);
      }
    }
  }

  // 2. Run superpositions before rendering. The view's camera is meaningless
  //    against un-aligned coords, so any saved quaternion expects this
  //    step to have happened. Each entry mutates the *moving* molecule's
  //    display transform in place; the reference is untouched.
  for (const sp of scene.superpose ?? []) {
    const mov = fileBindings.get(sp.move);
    const ref = fileBindings.get(sp.onto);
    if (!mov || !ref) {
      const missing = !mov ? sp.move : sp.onto;
      result.log.push({
        file: missing,
        domain: `superpose ${sp.move}→${sp.onto}`,
        message: `cannot superpose: file "${missing}" not bound`,
      });
      continue;
    }
    if (ref.molNo === undefined || ref.molNo === null) {
      result.log.push({
        file: sp.onto,
        domain: `superpose ${sp.move}→${sp.onto}`,
        message: `reference molecule has no molNo yet; skipped`,
      });
      continue;
    }
    try {
      await runSuperpose(sp, mov, ref);
    } catch (e) {
      console.warn(
        `[scene] superpose failed (${sp.method} ${sp.move}→${sp.onto}):`,
        e,
      );
      result.log.push({
        file: sp.move,
        domain: `superpose ${sp.move}→${sp.onto}`,
        message: `${sp.method} superpose failed: ${e instanceof Error ? e.message : "unknown error"}`,
      });
    }
  }

  // 3. Apply each element to its bound molecule.
  for (const element of scene.elements ?? []) {
    const mol = fileBindings.get(element.file);
    if (!mol) continue; // already noted in unresolvedFiles

    // Clear all existing custom representations on this molecule so the
    // scene fully controls the look. Built-in style toggles aren't custom
    // reps so this is non-destructive for things like default ribbons
    // added by the loader. Symmetrically tell the Redux store so the
    // Moorhen Models drawer drops them from its "Custom representations"
    // list.
    try {
      for (const rep of mol.representations ?? []) {
        if (rep.isCustom) {
          await mol.removeRepresentation(rep.uniqueId);
          dispatch(removeCustomRepresentation(rep));
        }
      }
    } catch (e) {
      console.warn(`[scene] failed to clear representations on ${mol.name}:`, e);
    }

    const reps = element.representations ?? [];

    // The scene fully owns the look of any molecule it touches: hide
    // every loader-default (non-custom) representation so only what the
    // scene asks for is drawn. Without this, every molecule that has a
    // scene element but no CRs/MolecularSurface in it (e.g. a fragment
    // showing only `style: ligands`) keeps its loader-default ribbon
    // visible — producing one ghost ribbon per loaded molecule.
    //
    // We hide rather than remove so the panel's built-in B/R/S toggles
    // can bring them back if the user wants.
    try {
      for (const rep of mol.representations ?? []) {
        if (!rep.isCustom && rep.visible) {
          rep.hide();
        }
      }
    } catch (e) {
      console.warn(`[scene] failed to hide overshadowed defaults on ${mol.name}:`, e);
    }

    let added = 0;
    for (const rep of reps) {
      const ok = await applyRepresentation({
        molecule: mol,
        rep,
        domains: scene.domains ?? [],
        fileName: element.file,
        log: result.log,
        policy,
        dispatch,
      });
      if (ok) added++;
    }
    result.applied.push({
      file: element.file,
      molNo: mol.molNo ?? -1,
      representations: added,
    });
  }

  // 2.5 Maps. For each scene.maps entry: bind to a loaded map if one
  //     matches the file ref; otherwise fetch+load via mapFetcher.
  //     Then dispatch contour/colour/style/visibility updates, and
  //     promote the named activeMap (if any).
  const liveMapPool: moorhen.Map[] = [...maps];
  const mapBindings = new Map<string, moorhen.Map>();
  for (const sceneMap of scene.maps ?? []) {
    const fileRef = (scene.files ?? []).find((f) => f.name === sceneMap.file);
    if (!fileRef) {
      result.log.push({
        file: sceneMap.file,
        domain: `map ${sceneMap.name}`,
        message: `map references unknown file "${sceneMap.file}"`,
      });
      continue;
    }
    const existing = matchOneMap(fileRef, liveMapPool);
    let map: moorhen.Map | null = existing;
    if (!map && mapFetcher && isFetchable(fileRef)) {
      try {
        map = await mapFetcher(fileRef, sceneMap);
        if (map) liveMapPool.push(map);
      } catch (e) {
        result.log.push({
          file: sceneMap.file,
          domain: `map ${sceneMap.name}`,
          message: `MTZ fetch failed: ${e instanceof Error ? e.message : "unknown"}`,
        });
        continue;
      }
    }
    if (!map) {
      result.log.push({
        file: sceneMap.file,
        domain: `map ${sceneMap.name}`,
        message: mapFetcher
          ? "no matching loaded map and ref is not fetchable"
          : "no matching loaded map (no mapFetcher provided)",
      });
      continue;
    }
    mapBindings.set(sceneMap.name, map);
    applyMapState(map, sceneMap, dispatch);
  }
  if (scene.activeMap) {
    const active = mapBindings.get(scene.activeMap);
    if (active) {
      dispatch(setActiveMap(active));
    } else {
      result.log.push({
        file: "",
        domain: `activeMap ${scene.activeMap}`,
        message: `activeMap "${scene.activeMap}" not bound to any loaded map`,
      });
    }
  }

  // 3. Camera. Mirror PasteViewLinkField exactly so behaviour matches.
  if (scene.view) {
    // Camera target — where the view points. `centre` centres on a selection's
    // centroid; `origin` is an explicit point. This is the ONLY thing that moves
    // the camera; the slab below is independent (z-depth only).
    if (scene.view.centre) {
      const { selection } = scene.view.centre;
      const cid = selection ?? "/*/*/*/*";
      const { mol, ref, error } = resolveViewMolecule(scene.view.centre.file, fileBindings);
      if (error || !mol) {
        result.log.push({
          file: ref,
          domain: "view.centre",
          message: `cannot centre: ${error}`,
        });
      } else if (cid.includes("||")) {
        // mol.centreOn feeds gemmi a single CID; for a "||"-union compute the
        // centroid ourselves and set the origin to it.
        const sphere = selectionBoundingSphere(mol, cid);
        if (sphere) {
          dispatch(setOrigin(sphere.centre));
        } else {
          result.log.push({
            file: ref,
            domain: "view.centre",
            message: `centre: selection "${cid}" matched no atoms (or gemmi unavailable)`,
          });
        }
      } else {
        try {
          await mol.centreOn(cid, false, false); // no animate, don't touch zoom
        } catch (e) {
          result.log.push({
            file: ref,
            domain: "view.centre",
            message: `centre on "${cid}" failed: ${e instanceof Error ? e.message : "unknown error"}`,
          });
        }
      }
    } else if (scene.view.origin) {
      dispatch(setOrigin(scene.view.origin));
    }

    // Slab: a z-depth clip/fog window ONLY (see the clip block below). It sizes
    // the depth to the selection's bounding radius but does NOT move the camera —
    // the clip brackets the current origin in depth, so a slab only "contains"
    // its selection when `centre` has put the origin on that selection. centre
    // and slab are independent and both settable; pair them to frame a selection.
    let slabDepth: number | undefined;
    if (scene.view.slab) {
      const { selection, pad } = scene.view.slab;
      const cid = selection ?? "/*/*/*/*";
      const { mol, ref, error } = resolveViewMolecule(scene.view.slab.file, fileBindings);
      if (error || !mol) {
        result.log.push({
          file: ref, domain: "view.slab",
          message: `cannot slab: ${error}`,
        });
      } else {
        const sphere = selectionBoundingSphere(mol, cid);
        if (!sphere) {
          result.log.push({
            file: ref, domain: "view.slab",
            message: `slab: selection "${cid}" matched no atoms (or gemmi unavailable)`,
          });
        } else {
          slabDepth = sphere.radius + (pad ?? 0);
        }
      }
    }
    if (scene.view.quat) dispatch(setQuat(scene.view.quat));
    if (scene.view.zoom !== undefined) dispatch(setZoom(scene.view.zoom));

    // Clip & fog. Coot derives both from zoom and a shared pair of field depths
    // (Å in front of / behind the centre, default 8/21) and recomputes them on
    // zoom UNLESS resetClippingFogging is off. So a scene that sets clip/fog also
    // pins that flag, and `clip` gives intent-level control over the depths.
    const clip = scene.view.clip;
    const zoom = scene.view.zoom ?? ctx.glRef?.zoom ?? 1;
    const fco = ctx.glRef?.fogClipOffset ?? 250;
    if (slabDepth !== undefined) {
      // Slab: clip/fog the selection's bounding sphere. clipStart/clipEnd are an
      // ABSOLUTE half-thickness in world Å about the fogClipOffset plane — the
      // eye-space z is never scaled by zoom (only x/y are, in the projection), so
      // the depth is the bounding radius (+pad) directly, NOT zoom*radius. (The
      // field-depth `clip` form below DOES multiply by zoom because there the
      // author gives coot's pre-zoom field depths, not an absolute distance.)
      dispatch(setClipStart(slabDepth));
      dispatch(setClipEnd(slabDepth));
      dispatch(setFogStart(fco - slabDepth));
      dispatch(setFogEnd(fco + slabDepth));
      dispatch(setResetClippingFogging(false));
    } else if (clip === "auto") {
      dispatch(setResetClippingFogging(true)); // let coot recompute from zoom
    } else if (clip && typeof clip === "object") {
      // Field depths in Å; clip = zoom*depth, fog offset by fogClipOffset —
      // coot's own formula with author depths. Drives clip AND fog together.
      dispatch(setClipStart(zoom * clip.front));
      dispatch(setClipEnd(zoom * clip.back));
      dispatch(setFogStart(fco - zoom * clip.front));
      dispatch(setFogEnd(fco + zoom * clip.back));
      dispatch(setResetClippingFogging(false));
    } else {
      // Explicit numbers (or `clip: "lock"`): set what's given and freeze, so
      // coot doesn't recompute them away on the next zoom.
      const setAny =
        scene.view.clipStart !== undefined || scene.view.clipEnd !== undefined ||
        scene.view.fogStart !== undefined || scene.view.fogEnd !== undefined;
      if (scene.view.clipStart !== undefined) dispatch(setClipStart(scene.view.clipStart));
      if (scene.view.clipEnd !== undefined) dispatch(setClipEnd(scene.view.clipEnd));
      if (scene.view.fogStart !== undefined) dispatch(setFogStart(scene.view.fogStart));
      if (scene.view.fogEnd !== undefined) dispatch(setFogEnd(scene.view.fogEnd));
      if (clip === "lock" || setAny) dispatch(setResetClippingFogging(false));
    }

    if (scene.view.background) {
      const rgba = hexToRgba01(scene.view.background);
      if (rgba) dispatch(setBackgroundColor(rgba));
    }
  }
  dispatch(setRequestDrawScene(true));

  return result;
}

// --------------------------------------------------------------------------
// File matching
// --------------------------------------------------------------------------

/**
 * True iff the resolver could ask the fetcher to load this ref on its
 * own — i.e. the ref carries enough info to know where to fetch from.
 * `path:` alone is not fetchable (we don't read arbitrary local paths
 * from the browser). `job+param` IS fetchable: the host fetcher resolves
 * the job number + output param to a project file via the ccp4i2 REST API
 * (jobs → files) and loads it through the same proxy URL as fileId refs.
 */
export function isFetchable(fr: SceneFileRef): boolean {
  if (fr.pdb) return true;
  if (fr.url) return true;
  if (fr.fileId !== undefined && fr.projectId) return true;
  if (fr.job !== undefined && fr.param) return true;
  // Inline dict text: trivially "fetchable" — the fetcher just returns
  // the text. Only valid on dictionary refs (validator enforces this).
  if (fr.cifText && fr.kind === "dictionary") return true;
  // Bundle ref: fetcher looks up bytes/text in an in-memory asset map.
  // Always considered fetchable; if the map doesn't have the key the
  // fetcher returns null and the resolver records it as unresolved.
  if (fr.bundle) return true;
  return false;
}

/** Match a file ref against the currently-loaded molecules. */
function matchOneFile(
  fr: SceneFileRef,
  molecules: moorhen.Molecule[],
): moorhen.Molecule | null {
  // 1. Match by ccp4i2 fileId extracted from uniqueId.
  if (fr.fileId !== undefined) {
    for (const mol of molecules) {
      const fid = extractFileIdFromUniqueId(mol.uniqueId || "");
      if (fid === fr.fileId) return mol;
    }
  }

  // 2. Match by PDB id embedded in the proxied download URL.
  //    Fetched PDB structures get uniqueId =
  //    /api/proxy/pdbe/entry-files/download/<pdbid>.cif
  if (fr.pdb) {
    const pdbLower = fr.pdb.toLowerCase();
    for (const mol of molecules) {
      const u = (mol.uniqueId || "").toLowerCase();
      if (u.includes(`/pdbe/entry-files/download/${pdbLower}`)) return mol;
    }
  }

  // 3. Match a previously-fetched bundle ref. The wrapper sets the
  //    molecule's uniqueId to `bundle:<assetPath>` so re-applying the
  //    same scene against the same loaded bundle skips the re-fetch.
  if (fr.bundle) {
    const sentinel = `bundle:${fr.bundle}`;
    for (const mol of molecules) {
      if (mol.uniqueId === sentinel) return mol;
    }
  }

  // 4. Match by URL or path (uniqueId is set to the loader's URL/path).
  const candidates = [fr.url, fr.path].filter(Boolean) as string[];
  for (const c of candidates) {
    for (const mol of molecules) {
      if (mol.uniqueId === c) return mol;
    }
  }

  return null;
}

/**
 * Same shape as matchOneFile but for maps. Maps are keyed by molNo via
 * extractFileIdFromUniqueId (when loaded from a ccp4i2 file), by URL,
 * or by the bundle: sentinel for re-applied bundle scenes.
 */
function matchOneMap(fr: SceneFileRef, maps: moorhen.Map[]): moorhen.Map | null {
  if (fr.fileId !== undefined) {
    for (const m of maps) {
      const fid = extractFileIdFromUniqueId(m.uniqueId || "");
      if (fid === fr.fileId) return m;
    }
  }
  if (fr.bundle) {
    const sentinel = `bundle:${fr.bundle}`;
    for (const m of maps) {
      if (m.uniqueId === sentinel) return m;
    }
  }
  const candidates = [fr.url, fr.path].filter(Boolean) as string[];
  for (const c of candidates) {
    for (const m of maps) {
      if (m.uniqueId === c) return m;
    }
  }
  return null;
}

/**
 * Apply a SceneMap's render state to a loaded map by dispatching the
 * relevant Moorhen actions. Difference maps split colour into positive
 * + negative lobes; non-difference use the single mapColours slot.
 *
 * No-op when a field is undefined — Moorhen defaults stand.
 */
function applyMapState(
  map: moorhen.Map,
  sceneMap: SceneMap,
  dispatch: Dispatch,
): void {
  const molNo = map.molNo;
  // Mask styling (incl. Coot's suggested contour level/radius) is applied by
  // applyMaskDefaults at load time in the wrapper's map fetcher; here we apply
  // only the scene's explicit overrides, for masks and ordinary maps alike.
  if (sceneMap.contourLevel !== undefined) {
    // setContourLevel's typings have drifted across Moorhen versions;
    // cast keeps the resolver compatible with the installed shape.
    dispatch(setContourLevel({ molNo, contourLevel: sceneMap.contourLevel } as never));
  }
  if (sceneMap.radius !== undefined) {
    dispatch(setMapRadius({ molNo, radius: sceneMap.radius } as never));
  }
  if (sceneMap.alpha !== undefined) {
    dispatch(setMapAlpha({ molNo, alpha: sceneMap.alpha } as never));
  }
  if (sceneMap.style) {
    dispatch(setMapStyle({ molNo, style: sceneMap.style } as never));
  }
  if (sceneMap.isDifference) {
    if (sceneMap.positiveColour) {
      const rgb = hexToRgb01(sceneMap.positiveColour);
      if (rgb) dispatch(setPositiveMapColours({ molNo, rgb } as never));
    }
    if (sceneMap.negativeColour) {
      const rgb = hexToRgb01(sceneMap.negativeColour);
      if (rgb) dispatch(setNegativeMapColours({ molNo, rgb } as never));
    }
  } else if (sceneMap.colour) {
    const rgb = hexToRgb01(sceneMap.colour);
    if (rgb) dispatch(setMapColours({ molNo, rgb } as never));
  }
  if (sceneMap.visible === false) {
    dispatch(hideMap(map as never));
  } else if (sceneMap.visible === true) {
    dispatch(showMap(map as never));
  }
}

/** Parse #rrggbb / #rrggbbaa into the {r,g,b} 0-1 shape Moorhen wants. */
function hexToRgb01(hex: string): { r: number; g: number; b: number } | null {
  const m = /^#([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})(?:[0-9a-fA-F]{2})?$/.exec(hex);
  if (!m) return null;
  return {
    r: parseInt(m[1], 16) / 255,
    g: parseInt(m[2], 16) / 255,
    b: parseInt(m[3], 16) / 255,
  };
}

// --------------------------------------------------------------------------
// Superposition
// --------------------------------------------------------------------------

/**
 * Apply one SceneSuperpose entry to a pair of already-loaded molecules.
 * Mutates `mov` in place; `ref` is untouched. Throws on coot-side
 * errors so the caller can record them in the resolver log.
 */
async function runSuperpose(
  sp: SceneSuperpose,
  mov: moorhen.Molecule,
  ref: moorhen.Molecule,
): Promise<void> {
  if (sp.method === "ssm") {
    await mov.SSMSuperpose(sp.movChain, ref.molNo as number, sp.refChain, true);
    return;
  }
  // LSQ: translate the scene's matches (or chain+range shorthand) into
  // Moorhen's lskqbResidueRangeMatch shape.
  const matches = expandLsqMatches(sp);
  const residueMatches: moorhen.lskqbResidueRangeMatch[] = [];
  for (const m of matches) {
    const refMatch = RANGE_RE.exec(m.refRange);
    const movMatch = RANGE_RE.exec(m.movRange);
    if (!refMatch || !movMatch) continue; // validator caught it; defensive skip
    residueMatches.push({
      refChainId: m.refChain,
      movChainId: m.movChain,
      refResidueRange: [parseInt(refMatch[1], 10), parseInt(refMatch[2], 10)],
      movResidueRange: [parseInt(movMatch[1], 10), parseInt(movMatch[2], 10)],
    });
  }
  const matchTypeMap = { all: 0, main: 1, ca: 2 } as const;
  const matchType = matchTypeMap[sp.matchType ?? "main"];
  await mov.lsqkbSuperpose(ref.molNo as number, residueMatches, matchType, true);
}

/**
 * Resolve an LSQ entry's matches: if the scene used the `chain`+`range`
 * shorthand, expand it into a single SceneLsqMatch with the same chain
 * and range on both sides. Otherwise use the explicit `matches` block.
 *
 * Exported for unit testing — the schema validation accepts the
 * shorthand, but the expansion logic is the bit most likely to drift.
 */
export function expandLsqMatches(sp: {
  matches?: SceneLsqMatch[];
  chain?: string;
  range?: string;
}): SceneLsqMatch[] {
  if (sp.chain && sp.range) {
    return [{
      refChain: sp.chain,
      movChain: sp.chain,
      refRange: sp.range,
      movRange: sp.range,
    }];
  }
  return sp.matches ?? [];
}

const RANGE_RE = /^(-?\d+)-(-?\d+)$/;

/**
 * Split a Coot-style multi-CID selection into single CIDs. Coot uses
 * "||" as the multi-CID separator in some APIs (rigid_body_fit,
 * copy_fragment_using_cid, etc.) — but MoorhenMoleculeRepresentation
 * does *not* split them, so the resolver has to do it before calling
 * addRepresentation, otherwise the whole "||"-joined string is treated
 * as a single (invalid) CID and silently selects nothing.
 *
 * Empty chunks (from leading/trailing/repeated "||") are dropped, and
 * each chunk is trimmed. A string with no "||" comes back as a single-
 * element array.
 */
export function splitMultiCid(cid: string): string[] {
  if (!cid.includes("||")) return [cid];
  return cid
    .split("||")
    .map((s) => s.trim())
    .filter((s) => s.length > 0);
}

// --------------------------------------------------------------------------
// Representations
// --------------------------------------------------------------------------

interface ApplyRepCtx {
  molecule: moorhen.Molecule;
  rep: SceneRepresentation;
  domains: SceneDomain[];
  fileName: string;
  log: SceneResolveLogEntry[];
  policy: "clamp-and-log" | "strict";
  dispatch: Dispatch;
}

/**
 * A colour rule expressed as the arguments to MoleculeRepresentation.addColourRule.
 * Lets us stage the rule independently of how Moorhen wires it together.
 */
interface PendingRule {
  ruleType: string;
  cid: string;
  color: string;
  args: (string | number)[];
  isMultiColourRule: boolean;
  applyColourToNonCarbonAtoms?: boolean;
}

async function applyRepresentation(ctx: ApplyRepCtx): Promise<boolean> {
  const { molecule, rep, dispatch, fileName, log } = ctx;

  const cid = rep.selection ?? "/*/*/*/*";

  // Coot supports "||"-joined multi-CID strings as a convention in *some*
  // APIs (copy_fragment, rigid_body_fit), but MoorhenMoleculeRepresentation
  // does NOT split them — it passes the whole string through as one CID,
  // which silently selects nothing. We split here and emit one
  // representation per chunk, all sharing the same style and colour rules.
  const cids = splitMultiCid(cid);

  const pendingRules = buildPendingRules(ctx, cid);

  let anySucceeded = false;
  for (const subCid of cids) {
    try {
      const created = await molecule.addRepresentation(
        rep.style as moorhen.RepresentationStyles,
        subCid,
        true, // isCustom — keeps it under our control to clear later
      );
      if (!created) continue;
      const hasAlpha = typeof rep.alpha === "number" && rep.alpha < 1;
      if (pendingRules.length > 0) {
        for (const r of pendingRules) {
          // Colour rule CID stays as authored — the rule's CID and the
          // representation's CID don't have to match (e.g. by-domain
          // rules target absolute residue numbers regardless of the
          // representation's selection).
          created.addColourRule(
            r.ruleType,
            r.cid,
            r.color,
            r.args,
            r.isMultiColourRule,
            r.applyColourToNonCarbonAtoms ?? false,
          );
        }
        // Rebuild the buffers with the colour rules applied.
        await molecule.redrawRepresentation(created.uniqueId);
      }
      // Opacity (Moorhen `nonCustomOpacity`, 0..1, 1=opaque; surfaces included)
      // MUST be applied LAST, after any redraw. setNonCustomOpacity rewrites the
      // alpha of the *current* buffers, but Moorhen's draw()/redraw() rebuild
      // buffers opaque and do NOT re-apply nonCustomOpacity — so setting it
      // before the redraw above is silently clobbered when the rebuild lands
      // (the "transparent while loading, opaque once drawn" flip). Setting it
      // here, on the final buffers, makes it stick; it requests its own redraw.
      if (hasAlpha) created.setNonCustomOpacity(rep.alpha as number);
      dispatch(addCustomRepresentation(created));
      anySucceeded = true;
    } catch (e) {
      const detail = extractRepError(e);
      console.warn(
        `[scene] failed to add ${rep.style} on ${molecule.name} (cid=${subCid}):`,
        e,
      );
      // Surface in the result log so the user sees what went wrong
      // without having to open DevTools. Coot's `Invalid selection
      // syntax` errors are the most common cause and the most
      // actionable for the author.
      log.push({
        file: fileName,
        domain: `${rep.style} ${subCid}`,
        message: detail,
      });
    }
  }
  return anySucceeded;
}

/**
 * Pull a useful one-line message out of whatever Moorhen/Coot threw at
 * us. The exception shape varies (sometimes a plain Error with stack,
 * sometimes a Coot-side Exception with an `Array` message), so we
 * sniff and reach for the most informative string we can find.
 */
function extractRepError(e: unknown): string {
  if (!e) return "unknown error";
  const anyE = e as { message?: unknown; stack?: unknown };
  // Coot exceptions surface their message as an array like
  // ["std::runtime_error", "Invalid selection syntax: //A/115,116"].
  if (Array.isArray(anyE.message)) {
    return anyE.message.filter((s) => typeof s === "string").join(": ");
  }
  if (typeof anyE.message === "string") return anyE.message;
  if (typeof anyE.stack === "string") {
    // Coot's "Error: std::runtime_error,Invalid selection syntax..." pattern.
    const m = anyE.stack.match(/Error:\s*(.+?)(?:\n|$)/);
    if (m) return m[1].slice(0, 200);
  }
  try {
    return String(e).slice(0, 200);
  } catch {
    return "unknown error";
  }
}

function buildPendingRules(ctx: ApplyRepCtx, defaultCid: string): PendingRule[] {
  const { molecule, rep, domains, fileName, log, policy } = ctx;
  if (!rep.colour) return [];

  if (Array.isArray(rep.colour)) {
    // Per-selection colour list: one single-colour rule per entry. A whole-chain
    // CID ("//A") and a residue range ("//A/121-130") apply identically here —
    // it's the general form by-domain compiles to, and what coot's default
    // per-chain colouring round-trips through.
    return rep.colour.map((c) => ({
      ruleType: "molecule",
      cid: c.selection,
      color: c.colour,
      args: [c.selection, c.colour],
      isMultiColourRule: false,
    }));
  }

  if (isSceneHexColour(rep.colour)) {
    // libcoot's add_colour_rule reads cid+colour from args, not from
    // this.cid/this.color (which are only consulted by the bond-style
    // shim_set_bond_colours path). Without [cid, colour] in args,
    // ribbons / MolecularSurface / etc. silently no-op.
    return [
      {
        ruleType: "molecule",
        cid: defaultCid,
        color: rep.colour,
        args: [defaultCid, rep.colour],
        isMultiColourRule: false,
      },
    ];
  }

  if (isSceneNamedColour(rep.colour)) {
    if (rep.colour === "by-domain") {
      return buildByDomainPendingRule(molecule, domains, fileName, log, policy);
    }
    // Named schemes (b-factor, af2-plddt, etc.) are Moorhen multi-rules
    // whose args are filled in by Moorhen at apply-time. We pass an empty
    // args array; Moorhen's internal getMultiColourRuleArgs supplies them.
    return [
      {
        ruleType: rep.colour,
        cid: defaultCid,
        color: "#ffffff",
        args: [],
        isMultiColourRule: true,
      },
    ];
  }

  if (isSceneRawColour(rep.colour)) {
    const raw = rep.colour.raw;
    return [
      {
        ruleType: raw.ruleType,
        cid: defaultCid,
        color: "#ffffff",
        args: raw.args,
        isMultiColourRule: raw.isMultiColourRule ?? true,
        applyColourToNonCarbonAtoms: raw.applyColourToNonCarbonAtoms,
      },
    ];
  }

  return [];
}

// --------------------------------------------------------------------------
// by-domain compilation: clamp ranges, build pipe-delimited args
// --------------------------------------------------------------------------

/**
 * Turn a domain CID `selection` into colour-rule segments. When the CID is the
 * `//chain/start-end` shape (a concrete chain + numeric range) it is clamped to
 * the residues present and warned about — diagnostic parity with the legacy
 * chain+range form. Any other CID (whole chain `//F`, residue names, atoms,
 * wildcard chains) passes straight through to Coot.
 */
function selectionToSegments(
  selection: string,
  color: string,
  domainName: string,
  presentByChain: Map<string, Set<number>>,
  log: SceneResolveLogEntry[],
  fileName: string,
  policy: "clamp-and-log" | "strict",
): string[] {
  const m = /^(\/[^/]*\/([^/*]+))\/(-?\d+)-(-?\d+)$/.exec(selection);
  if (!m) return [`${selection}^${color}`];
  const prefix = m[1];
  const chainId = m[2];
  const start = parseInt(m[3], 10);
  const end = parseInt(m[4], 10);
  const present = presentByChain.get(chainId);
  if (!present) {
    log.push({ file: fileName, domain: domainName, message: `chain ${chainId} not present in molecule; skipped` });
    return [];
  }
  const subRanges = clampRangeToPresent(start, end, present);
  if (subRanges.length === 0) {
    log.push({ file: fileName, domain: domainName, message: `range ${start}-${end} has no present residues in chain ${chainId}; skipped` });
    return [];
  }
  if (subRanges.length > 1 || subRanges[0][0] !== start || subRanges[0][1] !== end) {
    log.push({
      file: fileName,
      domain: domainName,
      message: `chain ${chainId}: range ${start}-${end} resolved to ${subRanges
        .map(([s, e]) => `${s}-${e}`)
        .join(", ")} after clamping to present residues`,
    });
    if (policy === "strict") {
      throw new Error(
        `Scene resolver in strict mode: domain "${domainName}" range ${start}-${end} not fully present in chain ${chainId}`,
      );
    }
  }
  return subRanges.map(([s, e]) => `${prefix}/${s}-${e}^${color}`);
}

function buildByDomainPendingRule(
  molecule: moorhen.Molecule,
  domains: SceneDomain[],
  fileName: string,
  log: SceneResolveLogEntry[],
  policy: "clamp-and-log" | "strict",
): PendingRule[] {
  if (domains.length === 0) return [];

  // Index present residues per chain from molecule.sequences.
  const presentByChain = new Map<string, Set<number>>();
  for (const seq of molecule.sequences ?? []) {
    const present = new Set<number>();
    for (const r of seq.sequence) present.add(r.resNum);
    presentByChain.set(seq.chain, present);
  }

  const segments: string[] = [];
  for (const d of domains) {
    // Preferred CID form: clamp the //chain/start-end shape (diagnostic parity
    // with chain+range), pass any other CID straight through.
    if (d.selection) {
      segments.push(
        ...selectionToSegments(
          d.selection, d.color, d.name, presentByChain, log, fileName, policy,
        ),
      );
      continue;
    }
    if (!d.chain) continue; // neither selection nor chain (invalid; validated upstream)

    // Legacy: resolve the chain selector into concrete chain ids.
    // - "*"      → every chain present in the structure
    // - "A"      → exactly chain A (legacy single-chain form)
    // - ["A","B"] → exactly those chains
    const targetChains = resolveChainSelector(d.chain, presentByChain);
    if (targetChains.length === 0) {
      log.push({
        file: fileName,
        domain: d.name,
        message:
          typeof d.chain === "string"
            ? `chain ${d.chain} not present in molecule; skipped`
            : `no chains from [${(d.chain as string[]).join(", ")}] present; skipped`,
      });
      continue;
    }

    // Range-less domain ⇒ the WHOLE chain: one `//chain` segment per chain, no
    // residue range (so a molecule's per-chain colouring is a set of whole-chain
    // "domains" adopted via colour: by-domain, the same path as range domains).
    if (!d.range) {
      for (const chainId of targetChains) segments.push(`//${chainId}^${d.color}`);
      continue;
    }

    const m = /^(-?\d+)-(-?\d+)$/.exec(d.range);
    if (!m) continue;
    const start = parseInt(m[1], 10);
    const end = parseInt(m[2], 10);

    for (const chainId of targetChains) {
      const present = presentByChain.get(chainId);
      if (!present) {
        // Should only happen with an explicit list naming an absent
        // chain; "*" already excludes absent ones, single-chain hit the
        // outer check.
        log.push({
          file: fileName,
          domain: d.name,
          message: `chain ${chainId} not present in molecule; skipped`,
        });
        continue;
      }

      const subRanges = clampRangeToPresent(start, end, present);
      if (subRanges.length === 0) {
        log.push({
          file: fileName,
          domain: d.name,
          message: `range ${start}-${end} has no present residues in chain ${chainId}; skipped`,
        });
        continue;
      }
      if (subRanges.length > 1 || subRanges[0][0] !== start || subRanges[0][1] !== end) {
        log.push({
          file: fileName,
          domain: d.name,
          message:
            `chain ${chainId}: range ${start}-${end} resolved to ${subRanges
              .map(([s, e]) => `${s}-${e}`)
              .join(", ")} after clamping to present residues`,
        });
      }
      if (
        policy === "strict" &&
        (subRanges.length > 1 || subRanges[0][0] !== start || subRanges[0][1] !== end)
      ) {
        throw new Error(
          `Scene resolver in strict mode: domain "${d.name}" range ${start}-${end} not fully present in chain ${chainId}`,
        );
      }

      for (const [s, e] of subRanges) {
        segments.push(`//${chainId}/${s}-${e}^${d.color}`);
      }
    }
  }

  if (segments.length === 0) return [];

  // One multi-rule, args = pipe-joined segments. Matches Moorhen's own
  // internal shape for multi-residue colouring (see secondary-structure
  // colouring in baby-gru/src/utils/utils.ts).
  return [
    {
      ruleType: "by-domain", // label only; not a built-in scheme
      cid: "/*/*/*/*",
      color: "#ffffff",
      args: [segments.join("|")],
      isMultiColourRule: true,
    },
  ];
}

/**
 * Resolve a SceneDomain.chain value into a concrete list of chain ids
 * for the loaded structure.
 *
 *   - "*"        → every chain that has at least one residue present
 *                  (sorted for stable output)
 *   - "A" etc.   → just ["A"]; the chain might be absent — the caller
 *                  is responsible for handling that case
 *   - ["A","B"]  → that list verbatim (in declared order)
 *
 * For "*", returns an empty list when no chains are present (e.g. an
 * un-loaded molecule) so callers can short-circuit with a single log
 * entry rather than emitting a misleading per-chain error.
 */
export function resolveChainSelector(
  selector: string | string[],
  presentByChain: Map<string, Set<number>>,
): string[] {
  if (typeof selector === "string") {
    if (selector === "*") {
      return [...presentByChain.keys()].sort();
    }
    return [selector];
  }
  return selector;
}

/**
 * Given an inclusive range [start, end] and the set of residue numbers
 * actually present, return one or more sub-ranges that:
 *   1. Are inside [start, end]
 *   2. Are contiguous in the present set
 *   3. Don't include any residues not present
 *
 * Sub-ranges are returned in ascending order.
 */
export function clampRangeToPresent(
  start: number,
  end: number,
  present: Set<number>,
): [number, number][] {
  const out: [number, number][] = [];
  let runStart: number | null = null;
  let runEnd: number | null = null;
  for (let n = start; n <= end; n++) {
    if (present.has(n)) {
      if (runStart === null) runStart = n;
      runEnd = n;
    } else if (runStart !== null) {
      out.push([runStart, runEnd!]);
      runStart = null;
      runEnd = null;
    }
  }
  if (runStart !== null) out.push([runStart, runEnd!]);
  return out;
}

// --------------------------------------------------------------------------
// Small utilities
// --------------------------------------------------------------------------

function hexToRgba01(hex: string): [number, number, number, number] | null {
  if (!/^#[0-9a-fA-F]{6}$/.test(hex)) return null;
  const r = parseInt(hex.slice(1, 3), 16) / 255;
  const g = parseInt(hex.slice(3, 5), 16) / 255;
  const b = parseInt(hex.slice(5, 7), 16) / 255;
  return [r, g, b, 1];
}
