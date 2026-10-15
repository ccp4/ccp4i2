"use client";

/**
 * Loading things into a CCP4i2 Moorhen viewer, and applying and capturing
 * scenes there: coordinates with their dictionaries, maps from MTZ
 * coefficients or map files, a job's outputs, a file by its database id,
 * and every scene file-ref kind. One implementation for every viewer page,
 * on top of use-moorhen-viewer-core.
 *
 * The rules that live here and must not be re-learnt per page:
 * - A molecule's dictionaries are attached to that molecule alone, before
 *   Moorhen looks for missing monomers, and never to Coot's global store
 *   (lib/moorhen-dictionaries). Which file defined each ligand is remembered
 *   so a captured scene can say so.
 * - A difference map is never the active map (Moorhen refines against it).
 * - A crystal's map file written as P1 (PanDDA's are) takes a loaded model's
 *   space group, or coot boxes it into the unit cell (lib/moorhen-map-file).
 * - A job that wrote a scene is the authority on how its outputs look.
 */
import {
  addMap,
  addMolecule,
  MoorhenMap,
  MoorhenMolecule,
  setActiveMap,
  setContourLevel,
  setMapColours,
  setMapRadius,
  setRequestDrawScene,
  showMolecule,
} from "moorhen/react-lib";
import { useCallback, useEffect, useRef } from "react";
import { moorhen } from "moorhen/types/moorhen";
import { apiArrayBuffer, apiGet, apiText } from "../../api-fetch";
import {
  buildableCompIds,
  COORDINATE_TYPES,
  DICTIONARY_TYPE,
  fetchCompanionDictionaryFiles,
  fetchDictionaryTexts,
  fetchJobDictionaryFiles,
  fetchSubjobDictionaryFiles,
  loadWithDictionaries,
  privateDictionaryNumber,
  provenanceOf,
  type DictionaryFile,
  type DictionaryToAttach,
} from "../../lib/moorhen-dictionaries";
import { baseName } from "../../lib/moorhen-path-load";
import { hexToRgb01 } from "../../lib/event-evidence";
import { parseScene, serialiseScene } from "../../lib/scene";
import {
  applyScene,
  SceneDictionaryFetcher,
  SceneDictionaryLoader,
  SceneFileFetcher,
  SceneMapFetcher,
  SceneMapMasker,
  SceneResolveResult,
} from "../../lib/moorhen-scene-resolver";
import type {
  MaskMap,
  MoorhenScene,
  SceneDomain,
  SceneFileRef,
  SceneSuperpose,
} from "../../types/moorhen-scene";
import type { SceneBundleAssets } from "./moorhen-scenes-panel";
import { readCameraState } from "../../lib/moorhen-view-state";
import { collectMapRenderState } from "../../lib/moorhen-map-render-state";
import {
  applyMaskDefaults,
  ccp4DodgeEmClamp,
  ccp4Mode0ToFloat,
  ccp4RestoreSpaceGroup,
  isMaskSubType,
  markMaskMap,
  primeEmMapHeaderInfo,
  primeXtalMapContourStats,
  requireMoorhenInstance,
  type CrystalSymmetry,
} from "../../lib/moorhen-map-file";
import {
  liftSceneStraight,
  MASK_PROVENANCE_PREFIX,
  promoteSceneToPortable,
  SceneLiftHints,
  SceneRefUrlResolver,
} from "../../lib/moorhen-scene-lifter";
import {
  EMPTY_PROJECT_CONTEXT,
  ProjectContext,
  SCENE_FILE_TYPE,
  sceneRefUrl,
} from "../../lib/scene-ref-urls";
import type { MoorhenViewerCore } from "./use-moorhen-viewer-core";

/** A file record as the REST API lists it, as far as the loaders read it. */
export interface FileRecord {
  id: number;
  name?: string;
  annotation?: string;
  job_param_name?: string;
  type: string;
  sub_type?: number | null;
  directory?: number;
}

export interface LoadStructureOptions {
  /** The molecule's own dictionaries (its job's). Attached to THIS molecule
   *  only, never to Coot's global store: a global entry is inherited by every
   *  molecule with none of its own, which is how another job's ligand of the
   *  same name ends up with this chemistry. */
  dictionaries?: DictionaryToAttach[];
  /** Centre the view on the molecule (default). A second molecule brought
   *  into a live view must not move the camera. */
  centre?: boolean;
}

export interface MapFileOptions {
  isMask?: boolean;
  /** The file's sub_type (1 normal, 2 difference, 3 anomalous, 4 mask):
   *  labels the map and, for 2 and 3, contours it as a difference map. */
  subType?: number;
  /** An absolute level to open at, e.g. PanDDA's optimal contour for an
   *  event map. Omitted: 3 rmsd for a difference map, else Moorhen's own. */
  contourLevel?: number | null;
  /** The file's own annotation, kept for the panel's tooltip: the contour
   *  rows are labelled by map type ("2Fo-Fc", "Fo-Fc"), which with several
   *  datasets loaded gives several identically-named rows. */
  description?: string;
  /** The map's row label in the control panel, where its type would mislead
   *  (a PanDDA event map is neither 2Fo-Fc nor Fo-Fc). */
  shortLabel?: string;
  /** "#rrggbb" */
  colour?: string;
  radius?: number;
}

/** What a job's load found, for a page that keeps state about it (the
 *  campaign viewer's ligand panel). */
export interface JobLoadInfo {
  dictionaryFiles: DictionaryFile[];
  dictionaries: DictionaryToAttach[];
  coordFile: FileRecord | null;
  molecule: moorhen.Molecule | null;
}

export interface JobLoadOptions {
  /** Centre the view on the job's model. Off for a job added to a live view. */
  centre?: boolean;
  /** Told what the load found, before the maps load. Not called when the job
   *  authored its own scene (applied instead). */
  onLoaded?: (info: JobLoadInfo) => void;
}

export interface MoorhenLoadersOptions {
  core: MoorhenViewerCore;
  /** The project the page is scoped to, for scene refs that need one
   *  (job+param, projectPath) and name none of their own. A ref: the page
   *  fills it as it learns the project, without re-binding the loaders. */
  projectContext?: React.RefObject<ProjectContext>;
  /** What a captured scene records as its project. */
  captureProject?: { projectId?: string; projectName?: string };
}

export type MoorhenLoaders = ReturnType<typeof useMoorhenLoaders>;

const fileUrl = (fileId: number) => `/api/proxy/ccp4i2/files/${fileId}/download/`;

export function useMoorhenLoaders({ core, projectContext, captureProject }: MoorhenLoadersOptions) {
  const {
    commandCentre,
    moorhenInstanceRef,
    store,
    dispatch,
    molecules,
    maps,
    backgroundColor,
    defaultBondSmoothness,
    monomerLibraryPath,
    getOrigin,
  } = core;

  const fallbackContext = useRef<ProjectContext>(EMPTY_PROJECT_CONTEXT);
  const contextRef = projectContext ?? fallbackContext;
  const currentContext = useCallback(
    (): ProjectContext => contextRef.current ?? EMPTY_PROJECT_CONTEXT,
    [contextRef],
  );
  // Read at capture time, so the handler does not re-bind per render.
  const captureProjectRef = useRef(captureProject);
  captureProjectRef.current = captureProject;

  // Per-molecule dictionary provenance: molNo → (comp_id → source project file).
  // Populated when a job's dicts + coords load together; read by Capture so the
  // lifter can emit terse fileId dict refs instead of inlining cifText.
  const dictSourcesRef = useRef<Map<number, Map<string, { fileId: number; projectId?: string }>>>(new Map());

  // -------------------------------------------------------------------------
  // Coordinates
  // -------------------------------------------------------------------------

  /**
   * Build a MoorhenMolecule from already-fetched coord text and register
   * it in the store. Use when bytes are in memory (e.g. from a scene
   * bundle); skips the apiText round-trip that loadStructure does.
   * Sets uniqueId from the supplied marker so callers can later identify
   * the molecule (e.g. `bundle:assets/foo.pdb`).
   */
  const loadStructureFromText = useCallback(async (
    coordText: string,
    molName: string,
    uniqueIdMarker: string,
    opts: LoadStructureOptions = {},
  ): Promise<moorhen.Molecule | null> => {
    if (!commandCentre.current) return null;
    const newMolecule = new MoorhenMolecule(requireMoorhenInstance(moorhenInstanceRef));
    newMolecule.setBackgroundColour(backgroundColor);
    newMolecule.defaultBondOptions.smoothness = defaultBondSmoothness;
    try {
      const dictionaries = opts.dictionaries ?? [];
      // The molecule's own dictionaries go on BEFORE Moorhen goes looking for
      // missing monomers; otherwise the ligand is bonded from whatever the
      // monomer library returns for its code (see lib/moorhen-dictionaries).
      await loadWithDictionaries(newMolecule as any, dictionaries, () =>
        newMolecule.loadToCootFromString(coordText, molName),
      );
      if (newMolecule.molNo === -1) {
        throw new Error("Cannot read the fetched molecule...");
      }
      newMolecule.uniqueId = uniqueIdMarker;
      // Remember which file defines each of this molecule's ligands, so a
      // captured scene lists them under this molecule's element.
      const provenance = provenanceOf(dictionaries);
      if (newMolecule.molNo != null && provenance.size > 0) {
        dictSourcesRef.current.set(newMolecule.molNo, provenance);
      }

      // Try ribbon representation first (better for protein overview)
      // Fall back to CBs if ribbons fail (e.g., no protein backbone)
      try {
        await newMolecule.addRepresentation("CRs", "/*/*/*/*");
      } catch {
        await newMolecule.addRepresentation("CBs", "/*/*/*/*");
      }

      // Always try to add ligand representation
      try {
        await newMolecule.addRepresentation("ligands", "/*/*/*/*");
      } catch {
        console.warn("[loadStructureFromText] Ligands representation failed");
      }

      if (opts.centre !== false) await newMolecule.centreOn("/*/*/*/*", false, true);
      dispatch(addMolecule(newMolecule));
      dispatch(showMolecule({ molNo: newMolecule.molNo } as any));
      return newMolecule;
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot parse coords for ${molName} (${uniqueIdMarker})`);
      return null;
    }
  }, [commandCentre, moorhenInstanceRef, store, monomerLibraryPath, backgroundColor, defaultBondSmoothness, dispatch]);

  const loadStructure = useCallback(async (
    url: string,
    molName: string,
    opts: LoadStructureOptions = {},
  ): Promise<moorhen.Molecule | null> => {
    try {
      const pdbData = await apiText(url);
      return loadStructureFromText(pdbData, molName, url, opts);
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch PDB entry from ${url}, doing nothing...`);
      return null;
    }
  }, [loadStructureFromText]);

  /** loadStructure under its older name; the same thing. */
  const fetchMolecule = loadStructure;

  // -------------------------------------------------------------------------
  // Maps
  // -------------------------------------------------------------------------

  /** A map from MTZ coefficients (F/PHI), as a job writes them. */
  const fetchMap = useCallback(async (
    url: string,
    mapName: string,
    mapSubType: number = 1,
    description?: string,
  ) => {
    if (!commandCentre.current) return;
    let loadedMap: moorhen.Map | undefined;
    // subType: 1=normal, 2=difference, 3=anomalous difference
    // Both difference and anomalous maps use isDifference=true for contouring
    const isDiffMap = mapSubType === 2 || mapSubType === 3;
    try {
      const mtzData = await apiArrayBuffer(url);
      const newMap = await MoorhenMap.loadToCootFromMtzData(
        new Uint8Array(mtzData),
        mapName,
        {
          F: "F",
          PHI: "PHI",
          useWeight: false,
          isDifference: isDiffMap,
        } as moorhen.selectedMtzColumns,
        requireMoorhenInstance(moorhenInstanceRef),
      );
      loadedMap = newMap;
      newMap.uniqueId = url;
      (newMap as any).mapSubType = mapSubType;
      if (description) (newMap as any).ccp4i2Description = description;
      // Before addMap: an EM-flagged MTZ map crashes the viewer otherwise.
      primeEmMapHeaderInfo(newMap);
      if (mapSubType === 3) {
        newMap.defaultPositiveMapColour = { r: 1.0, g: 0.65, b: 0.0 };
        newMap.defaultNegativeMapColour = { r: 0.6, g: 0.3, b: 0.8 };
      }
      if (newMap.molNo === -1)
        throw new Error("Cannot read the fetched map...");
      dispatch(addMap(newMap));
      // Only set as active map for non-difference maps (subType 1).
      // Difference maps (Fo-Fc = 2, anomalous = 3) must never be the
      // active map because Moorhen refines against the active map.
      if (!isDiffMap) {
        dispatch(setActiveMap(newMap));
      }
      // Reduce initial contour level to 0.8x for more sensitive screening
      const state = store.getState() as any;
      const contourLevels = state.mapContourSettings?.contourLevels || [];
      const entry = contourLevels.find((c: any) => c.molNo === newMap.molNo);
      const currentLevel = entry?.contourLevel;
      if (currentLevel != null && !Number.isNaN(currentLevel)) {
        const reducedLevel = currentLevel * 0.8;
        dispatch(setContourLevel({ molNo: newMap.molNo, contourLevel: reducedLevel } as any));
        newMap.drawMapContour().catch((err: Error) => {
          console.error("Failed to redraw map contour after level adjustment:", err);
        });
      }
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch map from ${url}`);
    }
    return loadedMap;
  }, [commandCentre, moorhenInstanceRef, store, dispatch]);

  /** The crystals of the loaded models (space group and cell), read from
   *  each molecule's gemmi structure; P1 and unreadable ones are skipped. */
  const loadedCrystals = useCallback((): CrystalSymmetry[] => {
    // Moorhen 1.0 exposes gemmi as window.gemmiModule; older builds as CCP4Module.
    type Gemmi = { get_spacegroup_by_name?: (n: string) => { ccp4: number; delete?: () => void } };
    const w = window as unknown as { gemmiModule?: Gemmi; CCP4Module?: Gemmi };
    const gemmi = w.gemmiModule ?? w.CCP4Module;
    const list = (store.getState() as moorhen.State).molecules.moleculeList ?? [];
    const out: CrystalSymmetry[] = [];
    for (const mol of list) {
      try {
        const st = (mol as unknown as {
          gemmiStructure?: { spacegroup_hm: string; cell: { a: number; b: number; c: number; alpha: number; beta: number; gamma: number } };
        }).gemmiStructure;
        if (!st?.spacegroup_hm || !gemmi?.get_spacegroup_by_name) continue;
        const sg = gemmi.get_spacegroup_by_name(st.spacegroup_hm);
        const ccp4 = sg?.ccp4;
        sg?.delete?.();
        if (!ccp4 || ccp4 <= 1) continue;
        const c = st.cell;
        out.push({ ccp4, cell: [c.a, c.b, c.c, c.alpha, c.beta, c.gamma] });
      } catch {
        // a molecule without a readable structure contributes nothing
      }
    }
    return out;
  }, [store]);

  /** A real-space CCP4 map file (application/CCP4-map): a mask, a PanDDA
   *  event or Z-map. Unlike fetchMap this uses loadToCootFromMapData. Masks
   *  get the shared translucent-solid defaults and are never the active map. */
  const fetchMapFile = useCallback(async (
    url: string,
    mapName: string,
    opts: MapFileOptions = {},
  ): Promise<moorhen.Map | undefined> => {
    if (!commandCentre.current) return undefined;
    let newMap: moorhen.Map | undefined;
    const isDiffMap = opts.subType === 2 || opts.subType === 3;
    try {
      // Convert mode-0 (int8) CCP4 maps to float so Moorhen reads sane stats
      // (no-op if already float). For masks, also nudge the P1/orthogonal cell
      // off 90° so coot contours periodically instead of clamping to the cell box.
      let mapData = ccp4Mode0ToFloat(await apiArrayBuffer(url));
      // A crystal's map written as P1 (PanDDA's are) takes the space group of
      // a loaded model on the same cell, or coot boxes it into the unit cell.
      if (!opts.isMask) mapData = ccp4RestoreSpaceGroup(mapData, loadedCrystals());
      if (opts.isMask) mapData = ccp4DodgeEmClamp(mapData);
      newMap = await MoorhenMap.loadToCootFromMapData(
        new Uint8Array(mapData),
        mapName,
        isDiffMap,
        requireMoorhenInstance(moorhenInstanceRef),
      );
      if (newMap.molNo === -1) throw new Error("Cannot read the fetched map file...");
      newMap.uniqueId = url;
      // Tag so the lifter captures it as a kind: "map" ref (not MTZ).
      (newMap as any).isCcp4MapFile = true;
      if (opts.subType != null) (newMap as any).mapSubType = opts.subType;
      if (opts.description) (newMap as any).ccp4i2Description = opts.description;
      if (opts.shortLabel) (newMap as any).ccp4i2ShortLabel = opts.shortLabel;
      if (isDiffMap) newMap.isDifference = true;
      if (opts.subType === 3) {
        newMap.defaultPositiveMapColour = { r: 1.0, g: 0.65, b: 0.0 };
        newMap.defaultNegativeMapColour = { r: 0.6, g: 0.3, b: 0.8 };
      }
      if (opts.isMask) {
        markMaskMap(newMap);
      }
      dispatch(addMap(newMap));
      if (opts.isMask) {
        await applyMaskDefaults(dispatch, newMap as any);
      }
      const molNo = newMap.molNo;
      if (opts.radius != null) {
        dispatch(setMapRadius({ molNo, radius: opts.radius } as never));
      }
      const rgb = opts.colour ? hexToRgb01(opts.colour) : null;
      if (rgb) dispatch(setMapColours({ molNo, rgb } as never));
      // Where to open: a level the producer recorded (an event map's optimal
      // contour, in absolute map units), else +/-3 rmsd for a difference map
      // (a Z-map's rmsd is ~1, so that is z = 3). A normal map keeps Moorhen's
      // own default.
      let level = opts.contourLevel ?? undefined;
      if (level == null && isDiffMap && !opts.isMask) {
        try {
          const rmsd = await newMap.fetchMapRmsd();
          if (Number.isFinite(rmsd) && rmsd > 0) level = 3 * rmsd;
        } catch (err) {
          console.warn("map rmsd unavailable:", err);
        }
      }
      if (level != null && Number.isFinite(level) && level > 0) {
        dispatch(setContourLevel({ molNo, contourLevel: level } as any));
        newMap.drawMapContour().catch((err: Error) => {
          console.error("Failed to redraw map contour at its starting level:", err);
        });
      }
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch map file from ${url}`);
      return undefined;
    }
    return newMap;
  }, [commandCentre, moorhenInstanceRef, store, dispatch, loadedCrystals]);

  // -------------------------------------------------------------------------
  // Dictionaries on their own
  // -------------------------------------------------------------------------

  /**
   * Build a molecule of each monomer a dictionary describes, at the centre
   * of the view, with that dictionary attached to it and to nothing else.
   *
   * Moorhen's own import reads the dictionary into Coot's global store and
   * builds from there; that entry would then be inherited by every molecule
   * with a residue of the same name. Here it is read under a private number
   * (see privateDictionaryNumber), and the dictionary is attached to the new
   * molecule once it has a number of its own. Returns the molecules built.
   */
  const loadMoleculesFromDictionary = useCallback(async (
    dict: DictionaryToAttach,
    uniqueIdBase: string,
    opts: { centre?: boolean } = {},
  ): Promise<moorhen.Molecule[]> => {
    const built: moorhen.Molecule[] = [];
    if (!commandCentre.current) return built;
    const compIds = buildableCompIds(dict.text);
    if (compIds.length === 0) return built;
    const privateNumber = privateDictionaryNumber();
    await commandCentre.current.cootCommand(
      {
        returnType: "status",
        command: "read_dictionary_string",
        commandArgs: [dict.text, privateNumber],
        changesMolecules: [],
      },
      false
    );
    const originCoords = getOrigin().map((coord: number) => -coord);
    const provenance = provenanceOf([dict]);
    for (const code of compIds) {
      const result = (await commandCentre.current.cootCommand(
        {
          returnType: "status",
          command: "get_monomer_and_position_at",
          commandArgs: [code, privateNumber, ...originCoords],
        },
        true
      )) as moorhen.WorkerResponse<number>;
      const molNo = result.data.result.result;
      if (result.data.result.status !== "Completed" || molNo == null || molNo < 0) {
        console.warn(`[dictionary] Coot could not build ${code} from its dictionary`);
        continue;
      }
      const newMolecule = new MoorhenMolecule(requireMoorhenInstance(moorhenInstanceRef));
      newMolecule.uniqueId = `${uniqueIdBase}#${code}`;
      newMolecule.molNo = molNo;
      newMolecule.name = code;
      newMolecule.setBackgroundColour(backgroundColor);
      newMolecule.defaultBondOptions.smoothness = defaultBondSmoothness;
      newMolecule.coordsFormat = "mmcif";
      await Promise.all([
        newMolecule.fetchDefaultColourRules(),
        newMolecule.addDict(dict.text),
      ]);
      if (provenance.size > 0) dictSourcesRef.current.set(molNo, provenance);
      if (opts.centre !== false && built.length === 0) {
        // Pass "" (not "/*/*/*/*"): centreAndAlignViewOn appends "*" and
        // "CA" to the CID internally, so passing "/*/*/*/*" produces the
        // malformed "/*/*/*/**" selector and a WebAssembly.Exception from
        // Coot. The empty-string branch in moorhen's implementation falls
        // back to the correct "/*/*/*/*" wildcard internally.
        newMolecule.centreAndAlignViewOn("", false, 100);
      }
      await newMolecule.fetchIfDirtyAndDraw("ligands");
      dispatch(addMolecule(newMolecule));
      dispatch(showMolecule({ molNo } as any));
      built.push(newMolecule);
    }
    return built;
  }, [commandCentre, moorhenInstanceRef, backgroundColor, defaultBondSmoothness, getOrigin, dispatch]);

  /**
   * A job that wrote dictionaries and no coordinates (ImportDictionary,
   * MakeLink run without a model): draw the monomers its dictionaries
   * describe. When they describe none -- MakeLink's is modifications and a
   * link -- try the dictionaries its subjobs wrote, where AcedrgLink leaves
   * the linked pair.
   */
  const loadDictionaryOnlyJob = useCallback(async (
    jobId: number,
    ownDictionaries: DictionaryFile[],
    opts: { centre?: boolean } = {},
  ): Promise<number> => {
    let built = 0;
    for (const files of [ownDictionaries, null]) {
      const candidates = files ?? (await fetchSubjobDictionaryFiles(jobId));
      for (const dict of await fetchDictionaryTexts(candidates, currentContext().uuid ?? undefined)) {
        const url = fileUrl(dict.fileId as number);
        const molecules = await loadMoleculesFromDictionary(dict, url, {
          centre: opts.centre !== false && built === 0,
        });
        built += molecules.length;
      }
      if (built > 0) break;
    }
    if (built > 0) dispatch(setRequestDrawScene(true));
    return built;
  }, [loadMoleculesFromDictionary, dispatch, currentContext]);

  // -------------------------------------------------------------------------
  // Scenes: apply
  // -------------------------------------------------------------------------

  // handleApplyScene is defined below (it depends on the scene fetchers);
  // the job loaders reach it through this ref.
  const applySceneRef = useRef<((yamlText: string) => Promise<unknown>) | null>(null);

  // Bundle assets (from an opened .scene.zip in the Scenes panel) live in a
  // ref so the stable-identity fetchers below can reach this apply's assets;
  // refreshed at the top of each apply.
  const bundleAssetsRef = useRef<SceneBundleAssets>(new Map());
  // What the last applied scene asked for that cannot be read back off the
  // live view: a superposition leaves no trace on a molecule (SSM/LSQ move
  // coords inside coot), a masked map is bytes with no operands, and
  // `colour: by-domain` compiles to a rule that keeps neither a domain's name
  // nor its authored range. Capture re-emits them (see LiftCtx).
  const lastAppliedSuperposeRef = useRef<SceneSuperpose[] | undefined>(undefined);
  const lastAppliedMaskMapsRef = useRef<MaskMap[] | undefined>(undefined);
  const lastAppliedDomainsRef = useRef<SceneDomain[] | undefined>(undefined);

  // Scene fetcher: invoked by the resolver for file refs that aren't
  // already loaded. Bundle assets short-circuit the network; every other
  // ref kind goes through sceneRefUrl. Returns the molecule that ends up in
  // the store, or null if the fetch fails (the resolver then leaves it in
  // unresolvedFiles).
  const handleFetchSceneFile: SceneFileFetcher = useCallback(
    async (ref: SceneFileRef, fetchOpts) => {
      // The element's own dictionaries, attached as the molecule loads.
      const loadOpts = { dictionaries: fetchOpts?.dictionaries ?? [] };
      if (ref.bundle) {
        const buf = bundleAssetsRef.current.get(ref.bundle);
        if (!buf) {
          console.warn(
            `[scene] bundle lookup miss: wanted "${ref.bundle}"; keys in map (${bundleAssetsRef.current.size}):`,
            Array.from(bundleAssetsRef.current.keys()).slice(0, 10),
          );
          return null;
        }
        try {
          const coordText = new TextDecoder("utf-8").decode(buf);
          return await loadStructureFromText(
            coordText,
            ref.name || ref.bundle,
            // Sentinel so matchOneFile recognises re-applies of the
            // same bundle ref instead of re-fetching.
            `bundle:${ref.bundle}`,
            loadOpts,
          );
        } catch (err) {
          console.warn(`[scene] failed to load bundle coord ${ref.bundle}:`, err);
          return null;
        }
      }
      const url = await sceneRefUrl(ref, currentContext());
      if (!url) return null;
      const fallbackName =
        ref.pdb ? ref.pdb.toLowerCase()
        : ref.fileId !== undefined ? `file_${ref.fileId}`
        : ref.job !== undefined ? `job_${ref.job}_${ref.param}`
        : ref.projectPath ? baseName(ref.projectPath)
        : url;
      return loadStructure(url, ref.name || fallbackName, loadOpts);
    },
    [loadStructure, loadStructureFromText, currentContext],
  );

  // Dictionary fetcher: returns raw CIF text for a `kind: dictionary`
  // ref. PDB id form is not allowed (validator should reject) — only
  // url, fileId, job+param, projectPath, cifText, and bundle are valid.
  // Dicts may contain multiple `data_comp_*` blocks; the resolver hands
  // the whole text to Coot in one call so all blocks get parsed in one shot.
  const handleFetchSceneDictionary: SceneDictionaryFetcher = useCallback(
    async (ref: SceneFileRef): Promise<string | null> => {
      // Inline text: no network needed.
      if (ref.cifText) return ref.cifText;
      if (ref.bundle) {
        const buf = bundleAssetsRef.current.get(ref.bundle);
        if (!buf) return null;
        try {
          return new TextDecoder("utf-8").decode(buf);
        } catch (err) {
          console.warn(`[scene] dictionary ${ref.name} bundle decode failed:`, err);
          return null;
        }
      }
      const url = await sceneRefUrl(ref, currentContext());
      if (!url) return null;
      try {
        return await apiText(url);
      } catch (err) {
        console.warn(`[scene] failed to fetch dictionary ${ref.name}:`, err);
        return null;
      }
    },
    [currentContext],
  );

  // Dictionary loader: tells Coot to parse the given CIF and associate
  // it with the given molNo (use -999999 for global).
  const handleLoadSceneDictionary: SceneDictionaryLoader = useCallback(
    async (dictText: string, molNo: number): Promise<void> => {
      if (!commandCentre.current) return;
      await commandCentre.current.cootCommand(
        {
          returnType: "status",
          command: "read_dictionary_string",
          commandArgs: [dictText, molNo],
          changesMolecules: molNo >= 0 ? [molNo] : [],
        },
        false,
      );
    },
    [commandCentre],
  );

  // Where a scene ref's bytes live. Shared by the map fetcher (apply path)
  // and the promoter (export path) so they agree on what counts as resolvable.
  const resolveSceneRefUrl: SceneRefUrlResolver = useCallback(
    async (ref) => sceneRefUrl(ref, currentContext()),
    [currentContext],
  );

  // Map fetcher: produces a loaded MoorhenMap for a scene `maps[]` entry.
  // Real-space CCP4 map files (kind: "map", incl. masks) load via
  // loadToCootFromMapData; MTZ refs via loadToCootFromMtzData with the column
  // spec. Bundle-asset bytes short-circuit the network.
  const handleFetchSceneMap: SceneMapFetcher = useCallback(
    async (ref, sceneMap) => {
      if (!commandCentre.current) return null;
      let bytes: ArrayBuffer | null = null;
      let uniqueId: string | null = null;
      if (ref.bundle) {
        const buf = bundleAssetsRef.current.get(ref.bundle);
        if (!buf) {
          console.warn(`[scene] map bundle miss: ${ref.bundle}`);
          return null;
        }
        bytes = buf;
        uniqueId = `bundle:${ref.bundle}`;
      } else {
        const url = await resolveSceneRefUrl(ref);
        if (!url) return null;
        try {
          bytes = await apiArrayBuffer(url);
          uniqueId = url;
        } catch (err) {
          console.warn(`[scene] map fetch failed for ${ref.name}:`, err);
          return null;
        }
      }
      try {
        const mapInstance = requireMoorhenInstance(moorhenInstanceRef);
        let newMap: moorhen.Map;
        if (ref.kind === "map") {
          // mode-0 -> float (sane stats); masks also dodge coot's EM cell-clamp.
          // A crystal's map written as P1 takes the space group of a loaded
          // model on its cell (scene models load before maps).
          const mapBytes = ccp4Mode0ToFloat(bytes as ArrayBuffer);
          newMap = await MoorhenMap.loadToCootFromMapData(
            new Uint8Array(sceneMap.isMask
              ? ccp4DodgeEmClamp(mapBytes)
              : ccp4RestoreSpaceGroup(mapBytes, loadedCrystals())),
            sceneMap.name,
            !!sceneMap.isDifference,
            mapInstance,
          );
          (newMap as any).isCcp4MapFile = true;
          if (sceneMap.isMask) markMaskMap(newMap);
        } else {
          const cols = sceneMap.columns ?? {};
          newMap = await MoorhenMap.loadToCootFromMtzData(
            new Uint8Array(bytes as ArrayBuffer),
            sceneMap.name,
            {
              F: cols.F,
              PHI: cols.PHI,
              Fobs: cols.Fobs,
              SigFobs: cols.SigFobs,
              FreeR: cols.FreeR,
              useWeight: !!cols.useWeight,
              calcStructFact: !!cols.calcStructFact,
              isDifference: !!sceneMap.isDifference,
            } as moorhen.selectedMtzColumns,
            mapInstance,
          );
        }
        if (newMap.molNo === -1) return null;
        if (uniqueId) newMap.uniqueId = uniqueId;
        // The panels label a map by its type ("2Fo-Fc", "Fo-Fc"); a scene's
        // map is called what the scene calls it ("Z-map", "Event 1 map").
        (newMap as any).ccp4i2ShortLabel = sceneMap.name;
        (newMap as any).ccp4i2Description = sceneMap.name;
        primeEmMapHeaderInfo(newMap);
        dispatch(addMap(newMap));
        if (ref.kind === "map" && sceneMap.isMask) {
          await applyMaskDefaults(dispatch, newMap as any);
        }
        return newMap;
      } catch (err) {
        console.warn(`[scene] failed to load map for ${ref.name}:`, err);
        return null;
      }
    },
    [commandCentre, moorhenInstanceRef, dispatch, resolveSceneRefUrl, loadedCrystals],
  );

  // Map masker: mask a source map by a model's atom selection, producing a
  // NEW MoorhenMap. Mirrors Moorhen's MapMasking menu item
  // (mask_map_by_atom_selection → new molNo), plus the addMap/settings
  // wiring the resolver can't do itself. radius omitted → -1 (Coot default).
  const handleMaskSceneMap: SceneMapMasker = useCallback(
    async ({ sourceMap, model, selection, radius, invert, name }) => {
      if (!commandCentre.current) return null;
      try {
        const result = (await commandCentre.current.cootCommand(
          {
            returnType: "status",
            command: "mask_map_by_atom_selection",
            commandArgs: [
              model.molNo,
              sourceMap.molNo,
              selection,
              radius ?? -1,
              invert,
            ],
          },
          false,
        )) as moorhen.WorkerResponse<number>;
        const newMolNo = result?.data?.result?.result;
        if (newMolNo == null || newMolNo === -1) return null;
        const newMap = new MoorhenMap(requireMoorhenInstance(moorhenInstanceRef));
        newMap.molNo = newMolNo;
        newMap.name = name;
        // Inherit the difference-map flag from the source; a mask of a
        // difference map is still a difference map.
        newMap.isDifference = !!(sourceMap as any).isDifference;
        // beta.1 replaced getSuggestedSettings(); populate rmsd/suggested level
        // directly (masked map has no fileHeader, so full initialise() can't run).
        await primeXtalMapContourStats(newMap);
        // Provenance stamp: a masked map's recipe is unrecoverable from its
        // grid bytes (Moorhen persists only the bytes). Stamp the uniqueId so
        // the capturer can round-trip maps WE masked as maskMaps recipes.
        // Maps masked outside the scene context carry no such marker and are
        // dropped-and-logged on capture. Keep in sync with the capturer's
        // MASK_PROVENANCE_PREFIX reader.
        newMap.uniqueId = `${MASK_PROVENANCE_PREFIX}${name}`;
        dispatch(addMap(newMap));
        return newMap;
      } catch (err) {
        console.warn(`[scene] map masking failed for "${name}":`, err);
        return null;
      }
    },
    [commandCentre, moorhenInstanceRef, dispatch],
  );

  // Apply a scene YAML: validate, then hand to the resolver with the
  // refs/state it needs. The optional assets map carries bytes for any
  // `bundle:` refs the YAML uses; the Scenes panel passes it when a
  // .scene.zip is loaded.
  const handleApplyScene = useCallback(
    async (
      yamlText: string,
      assets: SceneBundleAssets = new Map(),
    ): Promise<SceneResolveResult> => {
      // Refresh the bundle-assets ref so the (stable-identity) fetchers
      // see this apply's assets. Cleared back to empty by the next
      // non-bundled apply.
      bundleAssetsRef.current = assets;
      const scene = parseScene(yamlText);
      lastAppliedSuperposeRef.current = scene.superpose;
      lastAppliedMaskMapsRef.current = scene.maskMaps;
      lastAppliedDomainsRef.current = scene.domains;
      // Live camera snapshot for view.clip: { front, back } (clip = zoom*depth,
      // fog offset by fogClipOffset).
      const gl = readCameraState(store.getState() as moorhen.State);
      const result = await applyScene({
        scene,
        molecules,
        maps,
        dispatch,
        glRef: { zoom: gl.zoom, fogClipOffset: gl.fogClipOffset },
        fetcher: handleFetchSceneFile,
        dictionaryFetcher: handleFetchSceneDictionary,
        dictionaryLoader: handleLoadSceneDictionary,
        mapFetcher: handleFetchSceneMap,
        mapMasker: handleMaskSceneMap,
      });
      dispatch(setRequestDrawScene(true));
      return result;
    },
    [
      store,
      molecules,
      maps,
      dispatch,
      handleFetchSceneFile,
      handleFetchSceneDictionary,
      handleLoadSceneDictionary,
      handleFetchSceneMap,
      handleMaskSceneMap,
    ],
  );
  useEffect(() => {
    applySceneRef.current = handleApplyScene;
  }, [handleApplyScene]);

  // -------------------------------------------------------------------------
  // Files and jobs, by database record
  // -------------------------------------------------------------------------

  const fetchFile = useCallback(async (
    fileId: number,
    opts: {
      /** The dictionaries to attach to a coordinate file. Omitted: ask the
       *  server for the ones that belong with the file (its job's). */
      dictionaryFiles?: DictionaryFile[];
    } = {},
  ) => {
    const fileInfo: FileRecord | null = await apiGet(`files/${fileId}`);
    if (!fileInfo) {
      console.warn(`File with ID ${fileId} not found.`);
      return;
    }
    const url = fileUrl(fileId);
    if (COORDINATE_TYPES.has(fileInfo.type)) {
      const molName = fileInfo.annotation || fileInfo.job_param_name || fileInfo.name || `file_${fileId}`;
      // A coordinate file brings the dictionaries of the job it belongs to.
      const dictionaryFiles = opts.dictionaryFiles ?? (await fetchCompanionDictionaryFiles(fileId));
      await loadStructure(url, molName, { dictionaries: await fetchDictionaryTexts(dictionaryFiles) });
    } else if (fileInfo.type === "application/CCP4-mtz-map") {
      const molName = fileInfo.name || fileInfo.job_param_name || `file_${fileId}`;
      await fetchMap(url, molName, fileInfo.sub_type || 1, fileInfo.annotation || molName);
    } else if (fileInfo.type === SCENE_FILE_TYPE) {
      // A scene the job authored: apply it as written.
      await applySceneRef.current?.(await apiText(url));
    } else if (fileInfo.type === "application/CCP4-map") {
      const molName = fileInfo.annotation || fileInfo.name || fileInfo.job_param_name || `file_${fileId}`;
      await fetchMapFile(url, molName, {
        isMask: isMaskSubType(fileInfo.sub_type),
        subType: fileInfo.sub_type ?? undefined,
        description: molName,
      });
    } else if (fileInfo.type === DICTIONARY_TYPE) {
      await loadMoleculesFromDictionary({ text: await apiText(url), fileId }, url);
    }
  }, [loadStructure, fetchMap, fetchMapFile, loadMoleculesFromDictionary]);

  /**
   * Load a job's outputs into the viewer: its best coordinate file (mmCIF
   * over PDB) with the job's own dictionaries attached to that molecule
   * alone, then its maps.
   *
   * The job's dictionaries are its own and its inputs (a refinement's ligand
   * usually comes from an earlier acedrg job, which a filter over this job's
   * own files never saw). They go on this job's molecule and nothing else:
   * each coordinate set carries its own chemistry, and a molecule from a job
   * without dictionaries inherits none ("everything is called LIG").
   *
   * A job that authored a scene is the authority on how its outputs should
   * look: the scene is applied instead, and nothing is inferred.
   */
  const fetchJobFiles = useCallback(async (jobId: number, opts: JobLoadOptions = {}) => {
    if (!commandCentre.current) return;

    const files: FileRecord[] | null = await apiGet(`files/?job=${jobId}`);
    if (!files || !Array.isArray(files)) return;

    // Only job output files (directory=1), not imported files (directory=2)
    const jobOutputFiles = files.filter((f) => f.directory === 1);

    const sceneFiles = jobOutputFiles.filter((f) => f.type === SCENE_FILE_TYPE);
    if (sceneFiles.length > 0 && applySceneRef.current) {
      // The job-level scene (param SCENE) wins over per-item ones.
      const chosen = sceneFiles.find((f) => f.job_param_name === "SCENE") ?? sceneFiles[0];
      try {
        await applySceneRef.current(await apiText(fileUrl(chosen.id)));
        return;
      } catch (err) {
        console.warn("[fetchJobFiles] the job's scene could not be applied; inferring instead:", err);
      }
    }

    // STEP 1: the job's dictionaries, from the database.
    const dictionaryFiles = await fetchJobDictionaryFiles(jobId);
    const dictionaries = await fetchDictionaryTexts(dictionaryFiles, currentContext().uuid ?? undefined);

    // STEP 2: the best coordinate file (prefer mmCIF over PDB)
    const coordFiles = jobOutputFiles.filter((f) => COORDINATE_TYPES.has(f.type));
    const mmcifFile = coordFiles.find((f) => (f.name ?? "").toLowerCase().endsWith(".cif"));
    const coordFile = mmcifFile || coordFiles[0] || null;

    let loadedMolecule: moorhen.Molecule | null = null;
    if (coordFile) {
      const molName = coordFile.annotation || coordFile.job_param_name || coordFile.name || `job_${jobId}`;
      loadedMolecule = await loadStructure(fileUrl(coordFile.id), molName, {
        dictionaries,
        centre: opts.centre === true,
      });
    }
    opts.onLoaded?.({ dictionaryFiles, dictionaries, coordFile, molecule: loadedMolecule });

    // STEP 2b: No coordinates: draw the monomers the job's dictionaries
    // describe instead (ImportDictionary, MakeLink run without a model).
    // Imported files count here: an imported dictionary is ImportDictionary's result.
    if (!coordFile) {
      const ownDictionaries = files.filter((f) => f.type === DICTIONARY_TYPE);
      if (ownDictionaries.length > 0) {
        await loadDictionaryOnlyJob(jobId, ownDictionaries as DictionaryFile[], { centre: opts.centre === true });
      }
    }

    // STEP 3: every map (MTZ coefficients and real-space CCP4 maps / masks)
    for (const file of jobOutputFiles) {
      if (file.type === "application/CCP4-mtz-map") {
        const mapName = file.name || file.job_param_name || `job_${jobId}`;
        await fetchMap(fileUrl(file.id), mapName, file.sub_type || 1, file.annotation || mapName);
      } else if (file.type === "application/CCP4-map") {
        const mapName = file.annotation || file.name || file.job_param_name || `job_${jobId}`;
        await fetchMapFile(fileUrl(file.id), mapName, {
          isMask: isMaskSubType(file.sub_type),
          subType: file.sub_type ?? undefined,
          description: mapName,
        });
      }
    }

    // STEP 4: a molecule loaded with dictionaries is redrawn, so bond orders
    // and ligand geometry come from the per-molecule dictionary.
    if (loadedMolecule && dictionaries.length > 0) {
      try {
        await loadedMolecule.redraw();
      } catch (err) {
        console.warn("[fetchJobFiles] Redraw after dict association failed:", err);
      }
      dispatch(setRequestDrawScene(true));
    }
  }, [commandCentre, dispatch, fetchMap, fetchMapFile, loadStructure, loadDictionaryOnlyJob, currentContext]);

  // Handle map contour level changes from the control panel slider.
  // Synchronous, so slider dragging stays smooth: the redraw is not awaited.
  const handleMapContourLevelChange = useCallback(
    (molNo: number, level: number) => {
      // Despite the TypeScript def saying 'level', Moorhen reads 'contourLevel'.
      dispatch(setContourLevel({ molNo, contourLevel: level } as any));
      const map = maps.find((m) => m.molNo === molNo);
      if (map) {
        map.drawMapContour().catch((err) => {
          console.error("Failed to redraw map contour:", err);
        });
      }
      dispatch(setRequestDrawScene(true));
    },
    [dispatch, maps],
  );

  // -------------------------------------------------------------------------
  // Scenes: capture and export
  // -------------------------------------------------------------------------

  // Promote an editor YAML to a fully self-contained scene + asset bundle:
  // every URL-resolvable ref is fetched into assets/ and rewritten to
  // bundle: form; ligand dicts the captured YAML omits (library monomers and
  // custom in-memory dicts) are re-collected from the monomer library or the
  // molecule's own stored dict.
  const handlePromoteSceneToPortable = useCallback(
    async (
      yamlText: string,
      currentAssets: SceneBundleAssets,
    ): Promise<{ yamlText: string; assets: SceneBundleAssets; warnings: string[] }> => {
      const scene = parseScene(yamlText);
      const { scene: out, assets, warnings } = await promoteSceneToPortable({
        scene,
        existingAssets: currentAssets,
        resolveUrl: resolveSceneRefUrl,
        molecules,
        monomerLibraryPath: molecules[0]?.monomerLibraryPath,
      });
      return { yamlText: serialiseScene(out), assets, warnings };
    },
    [molecules, resolveSceneRefUrl],
  );

  // STRAIGHT lift: keep each file's real provenance (fileId+projectId / url /
  // path) so the captured YAML re-applies against the SAME source. Bundling
  // belongs only to the explicit "Save self-contained" action above. The
  // capture never inlines dict cifText: project ligand dicts travel as terse
  // fileId refs, and a dict without project provenance is left out.
  const handleCaptureScene = useCallback(async (): Promise<{
    scene: MoorhenScene;
    hints: SceneLiftHints;
    assets: SceneBundleAssets;
  }> => {
    const state = store.getState() as moorhen.State;
    // sceneSettings carries the effect toggles (SSAO / edge-detect / shadows /
    // depth-blur / perspective) the lifter folds into hints.effects.
    const sceneSettingsState = (state as unknown as {
      sceneSettings: {
        doSSAO?: boolean | null;
        ssaoRadius?: number | null;
        ssaoBias?: number | null;
        doEdgeDetect?: boolean | null;
        edgeDetectDepthThreshold?: number | null;
        edgeDetectNormalThreshold?: number | null;
        edgeDetectDepthScale?: number | null;
        edgeDetectNormalScale?: number | null;
        doShadow?: boolean | null;
        useOffScreenBuffers?: number | boolean | null;
        depthBlurRadius?: number | null;
        depthBlurDepth?: number | null;
        doPerspectiveProjection?: boolean | null;
      };
    }).sceneSettings;
    const activeMapMolNo =
      (state as unknown as { generalStates?: { activeMap?: moorhen.Map | null } })
        .generalStates?.activeMap?.molNo;
    const { scene, hints } = await liftSceneStraight({
      molecules,
      glRef: readCameraState(state),
      sceneSettings: sceneSettingsState,
      superpose: lastAppliedSuperposeRef.current,
      maskMaps: lastAppliedMaskMapsRef.current,
      domains: lastAppliedDomainsRef.current,
      projectId: captureProjectRef.current?.projectId,
      projectName: captureProjectRef.current?.projectName,
      // First molecule's monomerLibraryPath is the canonical Moorhen root
      // (all molecules share it). Without it the lifter falls back to
      // STANDARD_MONOMERS only.
      monomerLibraryPath: molecules[0]?.monomerLibraryPath,
      maps,
      // Moorhen keeps contour settings in parallel per-attribute lists keyed
      // by molNo; flattened to one record per map for the lifter.
      mapState: collectMapRenderState(state),
      activeMapMolNo: typeof activeMapMolNo === "number" ? activeMapMolNo : undefined,
      dictSources: dictSourcesRef.current,
    });
    return { scene, hints, assets: new Map() as SceneBundleAssets };
  }, [store, molecules, maps]);

  return {
    loadStructureFromText,
    loadStructure,
    fetchMolecule,
    fetchMap,
    fetchMapFile,
    loadedCrystals,
    loadMoleculesFromDictionary,
    loadDictionaryOnlyJob,
    fetchFile,
    fetchJobFiles,
    handleMapContourLevelChange,
    handleFetchSceneFile,
    handleFetchSceneDictionary,
    handleLoadSceneDictionary,
    handleFetchSceneMap,
    handleMaskSceneMap,
    resolveSceneRefUrl,
    handleApplyScene,
    applySceneRef,
    handleCaptureScene,
    handlePromoteSceneToPortable,
    dictSourcesRef,
    bundleAssetsRef,
  };
}
