"use client";

import {
  addMolecule,
  addMap,
  setActiveMap,
  setContourLevel,
  setTheme,
  setBackgroundColor,
  setRequestDrawScene,
  showMolecule,
  MoorhenContainer,
  MoorhenMolecule,
  MoorhenMap,
} from "moorhen/react-lib";
import { MoorhenInstanceProvider, MoorhenMenuSystem, setShownSidePanel } from "moorhen/react-lib";
// @ts-ignore - moorhen 0.23 type may lack .d.ts depending on build
import type { MoorhenPanel } from "moorhen/react-lib";

import {
  RefObject,
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { moorhen } from "moorhen/types/moorhen";
import { useDispatch, useSelector, useStore } from "react-redux";
import { webGL } from "moorhen/types/mgWebGL";
import { MoorhenCcp4i2TabbedPanel } from "./moorhen-ccp4i2-tabbed-panel";
import { apiGet, apiText, apiArrayBuffer, apiPost, apiUpload } from "../../api-fetch";
import { useTheme } from "../../theme/theme-provider";
import { useMoorhenViewState } from "../../hooks/use-moorhen-view-state";
import { parseScene, serialiseScene } from "../../lib/moorhen-scene";
import {
  applyScene,
  SceneFileFetcher,
  SceneMapFetcher,
  SceneResolveResult,
} from "../../lib/moorhen-scene-resolver";
import type { SceneFileRef } from "../../types/moorhen-scene";
import type { SceneBundleAssets } from "./moorhen-scenes-panel";
import {
  liftSceneToBundle,
  MapRenderState,
  promoteSceneToPortable,
  SceneLiftHints,
  SceneRefUrlResolver,
} from "../../lib/moorhen-scene-lifter";
import type { MoorhenScene } from "../../types/moorhen-scene";
import {
  MoorhenFallback,
  MoorhenErrorBoundary,
  useMoorhenCapabilities,
  isSafariBrowser,
} from "./moorhen-capability-check";
import { usePopcorn } from "../../providers/popcorn-provider";

export interface MoorhenWrapperProps {
  fileIds?: number[];
  viewParam?: string | null;
  jobId?: number | null;
}

const MoorhenWrapper: React.FC<MoorhenWrapperProps> = ({ fileIds, viewParam, jobId }) => {
  const capabilities = useMoorhenCapabilities();
  const [isSafari] = useState(() => isSafariBrowser());
  const { setMessage } = usePopcorn();
  const dispatch = useDispatch();
  const theme = useTheme();

  // View state hook for URL parameter support
  const { getViewUrl } = useMoorhenViewState({
    viewParam: viewParam ?? null,
    onViewRestored: () => {},
  });
  // Container ref for measuring available height below the AppBar
  const moorhenContainerRef = useRef<HTMLDivElement>(null);
  const setMoorhenDimensions = useCallback((): [number, number] => {
    if (moorhenContainerRef.current) {
      const top = moorhenContainerRef.current.getBoundingClientRect().top;
      return [window.innerWidth, window.innerHeight - top];
    }
    return [window.innerWidth, window.innerHeight];
  }, []);

  const glRef: RefObject<webGL.MGWebGL | null> = useRef(null);
  const commandCentre = useRef<null | moorhen.CommandCentre>(null);
  const moleculesRef = useRef<null | moorhen.Molecule[]>(null);
  const mapsRef = useRef<null | moorhen.Map[]>(null);
  const activeMapRef = useRef<moorhen.Map>(null);
  const lastHoveredAtom = useRef<null | moorhen.HoveredAtom>(null);
  const prevActiveMoleculeRef = useRef<null | moorhen.Molecule>(null);
  const timeCapsuleRef = useRef(null);
  const cootInitialized = useSelector(
    (state: moorhen.State) => state.generalStates.cootInitialized
  );
  const molecules = useSelector(
    (state: moorhen.State) => state.molecules.moleculeList
  );
  const maps = useSelector(
    (state: moorhen.State) => (state as unknown as { maps: moorhen.Map[] }).maps || []
  );
  const store = useStore();

  useEffect(() => {
    dispatch(
      setBackgroundColor(theme.mode === "light" ? [1, 1, 1, 1] : [0, 0, 0, 1])
    );
    dispatch(setTheme(theme.mode === "light" ? "flatly" : "darkly"));
  }, [theme.mode]);

  // Auto-open the (single) CCP4i2 side panel on mount. Multi-panel
  // registration in Moorhen's extension API has a bug where only one
  // tab renders even with both panels in extraSidePanels; sub-tabs for
  // scenes vs controls are handled inside the panel content via MUI Tabs
  // instead.
  useEffect(() => {
    dispatch(setShownSidePanel("ccp4i2Controls"));
  }, [dispatch]);

  const monomerLibraryPath =
    "https://raw.githubusercontent.com/MonomerLibrary/monomers/master/";

  const backgroundColor = useSelector(
    (state: moorhen.State) => state.sceneSettings.backgroundColor
  );
  const defaultBondSmoothness = useSelector(
    (state: moorhen.State) => state.sceneSettings.defaultBondSmoothness
  );

  // URL prefix for Moorhen to load its resources (CSS, pixmaps, monomers, etc.)
  // In web browsers, use API route for CORP headers (COEP compatibility)
  // In Electron, serve directly from public/MoorhenAssets
  const isElectron = typeof window !== "undefined" && !!(window as any).electronAPI;
  const urlPrefix = isElectron ? "/MoorhenAssets" : "/api/moorhen/MoorhenAssets";

  // Note: Don't subscribe to glRef state here - it changes every frame during rotation
  // and would cause constant re-renders. Access origin directly from store when needed.
  const getOrigin = useCallback(() => {
    const state = store.getState() as moorhen.State;
    return state.glRef.origin;
  }, [store]);

  /**
   * Core "fetch text, build MoorhenMolecule, register in store" path.
   * Returns the molecule so callers (notably the scene fetcher) can use
   * the result; fetchMolecule below wraps this and returns void to keep
   * existing callers' contract unchanged.
   */
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
  ): Promise<moorhen.Molecule | null> => {
    if (!commandCentre.current) return null;
    const newMolecule = new MoorhenMolecule(
      commandCentre as RefObject<moorhen.CommandCentre>,
      store as any,
      monomerLibraryPath
    );
    newMolecule.setBackgroundColour(backgroundColor);
    newMolecule.defaultBondOptions.smoothness = defaultBondSmoothness;
    try {
      await newMolecule.loadToCootFromString(coordText, molName);
      if (newMolecule.molNo === -1) {
        throw new Error("Cannot read the fetched molecule...");
      }
      newMolecule.uniqueId = uniqueIdMarker;

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

      await newMolecule.centreOn("/*/*/*/*", false, true);
      dispatch(addMolecule(newMolecule));
      dispatch(showMolecule({ molNo: newMolecule.molNo } as any));
      return newMolecule;
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot parse coords for ${molName} (${uniqueIdMarker})`);
      return null;
    }
  }, [commandCentre, store, monomerLibraryPath, backgroundColor, defaultBondSmoothness, dispatch]);

  const loadStructure = useCallback(async (
    url: string,
    molName: string,
  ): Promise<moorhen.Molecule | null> => {
    try {
      const pdbData = await apiText(url);
      return loadStructureFromText(pdbData, molName, url);
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch PDB entry from ${url}, doing nothing...`);
      return null;
    }
  }, [loadStructureFromText]);

  const fetchMolecule = useCallback(async (url: string, molName: string) => {
    await loadStructure(url, molName);
  }, [loadStructure]);

  const fetchMap = useCallback(async (
    url: string,
    mapName: string,
    mapSubType: number = 1
  ) => {
    if (!commandCentre.current) return;
    const newMap = new MoorhenMap(
      commandCentre as RefObject<moorhen.CommandCentre>,
      store as any
    );
    const isDiffMap = mapSubType === 2 || mapSubType === 3;
    try {
      const mtzData = await apiArrayBuffer(url);
      await newMap.loadToCootFromMtzData(new Uint8Array(mtzData), mapName, {
        F: "F",
        PHI: "PHI",
        useWeight: false,
        isDifference: isDiffMap,
      });
      newMap.uniqueId = url;
      (newMap as any).mapSubType = mapSubType;
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
    return newMap;
  }, [commandCentre, store, dispatch]);

  const fetchDict = useCallback(async (url: string) => {
    if (!commandCentre.current) return;
    const fileContent = await apiText(url);
    // Load all monomers in the dictionary into coot's global store
    await commandCentre.current.cootCommand(
      {
        returnType: "status",
        command: "read_dictionary_string",
        commandArgs: [fileContent, -999999],
        changesMolecules: [],
      },
      false
    );
    // Extract monomer codes from the CIF content (data_comp_XXX blocks, excluding comp_list)
    const monomerCodes: string[] = [];
    const blockPattern = /^data_comp_(\S+)/gm;
    let match;
    while ((match = blockPattern.exec(fileContent)) !== null) {
      const code = match[1];
      if (code !== "list") {
        monomerCodes.push(code);
      }
    }
    // Fallback: if no comp_ blocks found, try the legacy block name pattern
    if (monomerCodes.length === 0) {
      const legacyPattern = /^data_(\S+)/gm;
      while ((match = legacyPattern.exec(fileContent)) !== null) {
        monomerCodes.push(match[1]);
      }
    }
    const originCoords = getOrigin().map((coord: number) => -coord);
    let centredFirst = false;
    for (const code of monomerCodes) {
      const result = (await commandCentre.current.cootCommand(
        {
          returnType: "status",
          command: "get_monomer_and_position_at",
          commandArgs: [code, -999999, ...originCoords],
        },
        true
      )) as moorhen.WorkerResponse<number>;
      if (result.data.result.status === "Completed") {
        const newMolecule = new MoorhenMolecule(
          commandCentre as RefObject<moorhen.CommandCentre>,
          store as any,
          monomerLibraryPath
        );
        newMolecule.uniqueId = `${url}#${code}`;
        newMolecule.molNo = result.data.result.result;
        newMolecule.name = code;
        newMolecule.setBackgroundColour(backgroundColor);
        newMolecule.defaultBondOptions.smoothness = defaultBondSmoothness;
        newMolecule.coordsFormat = "mmcif";
        await Promise.all([
          newMolecule.fetchDefaultColourRules(),
          newMolecule.addDict(fileContent),
        ]);
        if (!centredFirst) {
          // Pass "" (not "/*/*/*/*"): centreAndAlignViewOn appends "*" and
          // "CA" to the CID internally, so passing "/*/*/*/*" produces the
          // malformed "/*/*/*/**" selector and a WebAssembly.Exception from
          // Coot. The empty-string branch in moorhen's implementation falls
          // back to the correct "/*/*/*/*" wildcard internally.
          newMolecule.centreAndAlignViewOn("", false, 100);
          centredFirst = true;
        }
        await newMolecule.fetchIfDirtyAndDraw("ligands");
        dispatch(addMolecule(newMolecule));
        dispatch(showMolecule({ molNo: newMolecule.molNo } as any));
      }
    }
  }, [commandCentre, store, monomerLibraryPath, backgroundColor, defaultBondSmoothness, getOrigin, dispatch]);

  const fetchFile = useCallback(async (fileId: number) => {
    const fileInfo = await apiGet(`files/${fileId}`);
    if (!fileInfo) {
      console.warn(`File with ID ${fileId} not found.`);
      return;
    }
    if (fileInfo.type === "chemical/x-pdb") {
      const url = `/api/proxy/ccp4i2/files/${fileId}/download/`;
      const molName = fileInfo.annotation || fileInfo.job_param_name || fileInfo.name || `file_${fileId}`;
      await fetchMolecule(url, molName);
    } else if (fileInfo.type === "application/CCP4-mtz-map") {
      const url = `/api/proxy/ccp4i2/files/${fileId}/download/`;
      const molName = fileInfo.name || fileInfo.job_param_name;
      const mapSubType = fileInfo.sub_type || 1;
      await fetchMap(url, molName, mapSubType);
    } else if (fileInfo.type === "application/refmac-dictionary") {
      const url = `/api/proxy/ccp4i2/files/${fileId}/download/`;
      await fetchDict(url);
    }
  }, [fetchMolecule, fetchMap, fetchDict]);

  /**
   * Load all output files from a job into Moorhen.
   *
   * Dictionaries are loaded globally first so coot can parse ligand geometry
   * in coordinates. After loading coordinates, dictionaries are re-associated
   * with the specific molecule molNo (not the global -999999). This handles
   * the "everything is called LIG/DRG" problem in fragment campaigns: each
   * coordinate set gets its own dictionary association, so different ligand
   * geometries with the same residue name coexist correctly.
   */
  const fetchJobFiles = useCallback(async (jobId: number) => {
    if (!commandCentre.current) return;

    const files = await apiGet(`files/?job=${jobId}`);
    if (!files || !Array.isArray(files)) return;

    // Only job output files (directory=1), not imported files (directory=2)
    const jobOutputFiles = files.filter((f: { directory: number }) => f.directory === 1);

    // STEP 1: Load all ligand dictionaries into coot's global store FIRST.
    // This ensures coot understands ligand geometry when parsing coordinates.
    const dictFiles = files.filter(
      (f: { type: string }) => f.type === "application/refmac-dictionary"
    );
    const dictContents: string[] = [];
    for (const dictFile of dictFiles) {
      const dictUrl = `/api/proxy/ccp4i2/files/${dictFile.id}/download/`;
      try {
        const content = await apiText(dictUrl);
        // Load globally so coordinate parsing works
        await commandCentre.current.cootCommand(
          {
            returnType: "status",
            command: "read_dictionary_string",
            commandArgs: [content, -999999],
            changesMolecules: [],
          },
          false
        );
        dictContents.push(content);
      } catch (err) {
        console.warn("[fetchJobFiles] Failed to load dictionary:", err);
      }
    }

    // STEP 2: Load the best coordinate file (prefer mmCIF over PDB)
    const coordFiles = jobOutputFiles.filter(
      (f: { type: string }) =>
        f.type === "chemical/x-pdb" ||
        f.type === "chemical/x-cif" ||
        f.type === "chemical/x-mmcif"
    );
    const mmcifFile = coordFiles.find((f: { name: string }) =>
      f.name.toLowerCase().endsWith(".cif")
    );
    const coordFile = mmcifFile || coordFiles[0];

    let loadedMolecule: moorhen.Molecule | null = null;
    if (coordFile) {
      const url = `/api/proxy/ccp4i2/files/${coordFile.id}/download/`;
      const molName = coordFile.annotation || coordFile.job_param_name || coordFile.name || `job_${jobId}`;
      const newMolecule = new MoorhenMolecule(
        commandCentre as RefObject<moorhen.CommandCentre>,
        store as any,
        monomerLibraryPath
      );
      newMolecule.setBackgroundColour(backgroundColor);
      newMolecule.defaultBondOptions.smoothness = defaultBondSmoothness;
      try {
        const pdbData = await apiText(url);
        await newMolecule.loadToCootFromString(pdbData, molName);
        if (newMolecule.molNo === -1) throw new Error("Cannot read coordinates");
        newMolecule.uniqueId = url;

        // STEP 2b: Re-associate dictionaries with THIS molecule's molNo.
        // This overrides the global (-999999) association so that this
        // coordinate set's ligands use the correct geometry even when
        // other coordinate sets have ligands with the same residue name.
        for (const dictContent of dictContents) {
          await commandCentre.current!.cootCommand(
            {
              returnType: "status",
              command: "read_dictionary_string",
              commandArgs: [dictContent, newMolecule.molNo],
              changesMolecules: [newMolecule.molNo!],
            },
            false
          );
          await newMolecule.addDict(dictContent);
        }

        // Add representations — ribbon + ligands
        try {
          await newMolecule.addRepresentation("CRs", "/*/*/*/*");
        } catch {
          await newMolecule.addRepresentation("CBs", "/*/*/*/*");
        }
        try {
          await newMolecule.addRepresentation("ligands", "/*/*/*/*");
        } catch {
          console.warn("[fetchJobFiles] Ligands representation failed");
        }

        dispatch(addMolecule(newMolecule));
        dispatch(showMolecule({ molNo: newMolecule.molNo } as any));
        loadedMolecule = newMolecule;
      } catch (err) {
        console.warn("[fetchJobFiles] Failed to load coordinates:", err);
      }
    }

    // STEP 3: Load all maps
    for (const file of jobOutputFiles) {
      if (file.type === "application/CCP4-mtz-map") {
        const url = `/api/proxy/ccp4i2/files/${file.id}/download/`;
        const mapName = file.name || file.job_param_name;
        const mapSubType = file.sub_type || 1;
        await fetchMap(url, mapName, mapSubType);
      }
    }

    // STEP 4: If we loaded a molecule with dictionaries, redraw to pick up
    // correct bond orders / ligand geometry from the per-molecule dict.
    if (loadedMolecule && dictContents.length > 0) {
      try {
        await loadedMolecule.redraw();
      } catch (err) {
        console.warn("[fetchJobFiles] Redraw after dict association failed:", err);
      }
      dispatch(setRequestDrawScene(true));
    }
  }, [commandCentre, store, monomerLibraryPath, backgroundColor, defaultBondSmoothness, dispatch, fetchMap, getOrigin]);

  // Handle map contour level changes from the control panel slider
  const handleMapContourLevelChange = useCallback(
    (molNo: number, level: number) => {
      dispatch(setContourLevel({ molNo, contourLevel: level } as any));
      const map = maps.find((m) => m.molNo === molNo);
      if (map) {
        map.drawMapContour().catch((err) => {
          console.error("Failed to redraw map contour:", err);
        });
      }
      dispatch(setRequestDrawScene(true));
    },
    [dispatch, maps]
  );

  // Run servalcat refinement using the current job as context
  const [servalcatStatus, setServalcatStatus] = useState<string | null>(null);
  const handleRunServalcat = useCallback(
    async (mol: moorhen.Molecule) => {
      if (!jobId) {
        setServalcatStatus("No source job — open from a job to use refinement");
        return;
      }
      try {
        setServalcatStatus("Creating servalcat job...");

        // Look up the job's UUID and project
        const jobInfo = await apiGet(`jobs/${jobId}`);
        if (!jobInfo?.uuid || !jobInfo?.project) {
          setServalcatStatus("Could not retrieve job information");
          return;
        }

        // Create servalcat_pipe job with context_job_uuid for auto-population
        const createResult = await apiPost(
          `projects/${jobInfo.project}/create_task/`,
          {
            task_name: "servalcat_pipe",
            title: `Servalcat refinement of ${mol.name}`,
            context_job_uuid: jobInfo.uuid,
          }
        );
        const newJobId = createResult?.data?.new_job?.id;
        if (!newJobId) {
          setServalcatStatus("Failed to create servalcat job");
          return;
        }

        // Upload Moorhen coordinates as XYZIN (preserve original format)
        setServalcatStatus("Uploading coordinates...");
        const coordText = await mol.getAtoms();
        const ismmCIF = (mol as any).coordsFormat === "mmcif";
        const mimeType = ismmCIF ? "chemical/x-cif" : "chemical/x-pdb";
        const ext = ismmCIF ? ".cif" : ".pdb";
        const coordBlob = new Blob([coordText], { type: mimeType });
        const coordFile = new File([coordBlob], `${mol.name || "coords"}${ext}`);
        const formData = new FormData();
        formData.append("file", coordFile);
        formData.append("object_path", "servalcat_pipe.inputData.XYZIN");
        await apiUpload(`jobs/${newJobId}/upload_file_param/`, formData);

        // Run the job
        setServalcatStatus("Running servalcat refinement...");
        await apiPost(`jobs/${newJobId}/run/`, {});

        setServalcatStatus("Servalcat refinement submitted");
        // Clear status after a few seconds
        setTimeout(() => setServalcatStatus(null), 4000);
      } catch (err) {
        console.error("Failed to run servalcat:", err);
        setServalcatStatus(
          `Failed: ${err instanceof Error ? err.message : "Unknown error"}`
        );
      }
    },
    [jobId]
  );

  // Scene fetcher: invoked by the resolver for file refs that aren't
  // already loaded. Routes by ref shape:
  //   pdb:                       → PDBe proxy (mmCIF)
  //   fileId + projectId         → ccp4i2 proxy (per-project file id)
  //   url:                       → direct URL fetch
  // Returns the molecule that ends up in the store, or null if the fetch
  // path fails (the resolver then leaves it in unresolvedFiles).
  //
  // Bundle assets are kept in a ref that handleApplyScene refreshes
  // before each apply — keeps the fetcher callbacks stable while
  // letting per-apply asset maps reach them.
  const bundleAssetsRef = useRef<SceneBundleAssets>(new Map());

  const handleFetchSceneFile: SceneFileFetcher = useCallback(
    async (ref: SceneFileRef) => {
      // Bundle: decode bytes from the in-memory asset map and hand the
      // text straight to loadStructureFromText — no network round-trip,
      // no Blob URL (which would get prefixed by the ccp4i2 api-fetch
      // helper and 404).
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
          );
        } catch (err) {
          console.warn(`[scene] failed to load bundle coord ${ref.bundle}:`, err);
          return null;
        }
      }
      if (ref.pdb) {
        const pdbId = ref.pdb.toLowerCase();
        const url = `/api/proxy/pdbe/entry-files/download/${pdbId}.cif`;
        return loadStructure(url, ref.name || pdbId);
      }
      if (ref.fileId !== undefined && ref.projectId) {
        const url = `/api/proxy/ccp4i2/files/${ref.fileId}/download/`;
        return loadStructure(url, ref.name || `file_${ref.fileId}`);
      }
      if (ref.url) {
        return loadStructure(ref.url, ref.name || ref.url);
      }
      return null;
    },
    [loadStructure],
  );

  // Dictionary fetcher: returns raw CIF text for a `kind: dictionary`
  // ref. PDB id form is not allowed (validator should reject) — only
  // url, path, fileId+projectId, cifText, and bundle are valid for dicts.
  // Dicts may contain multiple `data_comp_*` blocks; the resolver hands
  // the whole text to Coot in one call so all blocks get parsed in one shot.
  const handleFetchSceneDictionary = useCallback(
    async (ref: SceneFileRef): Promise<string | null> => {
      // Inline text: no network needed. The lifter produces these for
      // dicts loaded from job outputs where there's no stable URL.
      if (ref.cifText) return ref.cifText;
      // Bundle: decode the asset bytes as UTF-8 text.
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
      let url: string | null = null;
      if (ref.fileId !== undefined && ref.projectId) {
        url = `/api/proxy/ccp4i2/files/${ref.fileId}/download/`;
      } else if (ref.url) {
        url = ref.url;
      }
      if (!url) return null;
      try {
        return await apiText(url);
      } catch (err) {
        console.warn(`[scene] failed to fetch dictionary ${ref.name}:`, err);
        return null;
      }
    },
    [],
  );

  // Dictionary loader: tells Coot to parse the given CIF and associate
  // it with the given molNo (use -999999 for global). Matches the
  // existing fragment-campaign pattern in fetchJobFiles.
  const handleLoadSceneDictionary = useCallback(
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
    [],
  );

  // Centralised URL-derivation for any URL-shaped ref kind. Shared by
  // the fetcher (apply path), the map fetcher, and the promoter (export
  // path) so they agree on what counts as "a resolvable ref". Declared
  // here (above the consumers) so the closures pick it up cleanly.
  const resolveSceneRefUrl: SceneRefUrlResolver = useCallback((ref) => {
    if (ref.pdb) return `/api/proxy/pdbe/entry-files/download/${ref.pdb.toLowerCase()}.cif`;
    if (ref.fileId !== undefined && ref.projectId) {
      return `/api/proxy/ccp4i2/files/${ref.fileId}/download/`;
    }
    if (ref.url) return ref.url;
    return null;
  }, []);

  // Map fetcher: produces a loaded MoorhenMap for a scene `kind: "mtz"`
  // ref + SceneMap column spec. Bundle-asset bytes short-circuit the
  // network; URL-resolvable refs go via the existing scene-file URL
  // logic. Always uses loadToCootFromMtzData (taking raw bytes) so the
  // bundle path doesn't need a Blob round-trip.
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
        const url = resolveSceneRefUrl(ref);
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
        const newMap = new MoorhenMap(
          commandCentre as RefObject<moorhen.CommandCentre>,
          store as any,
        );
        await newMap.loadToCootFromMtzData(
          new Uint8Array(bytes as ArrayBuffer),
          sceneMap.name,
          {
            F: sceneMap.columns.F,
            PHI: sceneMap.columns.PHI,
            Fobs: sceneMap.columns.Fobs,
            SigFobs: sceneMap.columns.SigFobs,
            FreeR: sceneMap.columns.FreeR,
            useWeight: !!sceneMap.columns.useWeight,
            calcStructFact: !!sceneMap.columns.calcStructFact,
            isDifference: !!sceneMap.isDifference,
          } as moorhen.selectedMtzColumns,
        );
        if (newMap.molNo === -1) return null;
        if (uniqueId) newMap.uniqueId = uniqueId;
        dispatch(addMap(newMap));
        return newMap;
      } catch (err) {
        console.warn(`[scene] failed to load MTZ for ${ref.name}:`, err);
        return null;
      }
    },
    [commandCentre, store, dispatch, resolveSceneRefUrl],
  );

  // Apply a scene YAML: validate, then hand to the resolver with the
  // refs/state it needs. Defined here because the wrapper owns dispatch
  // and the command-centre ref. The optional assets map carries bytes
  // for any `bundle:` refs the YAML uses; the Scenes panel passes it
  // when a .scene.zip is loaded.
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
      return applyScene({
        scene,
        molecules,
        maps,
        dispatch,
        fetcher: handleFetchSceneFile,
        dictionaryFetcher: handleFetchSceneDictionary,
        dictionaryLoader: handleLoadSceneDictionary,
        mapFetcher: handleFetchSceneMap,
      });
    },
    [
      molecules,
      maps,
      dispatch,
      handleFetchSceneFile,
      handleFetchSceneDictionary,
      handleLoadSceneDictionary,
      handleFetchSceneMap,
    ],
  );

  // Resolve project context from jobId so downloaded scenes can include
  // a stable projectId on their file references. Looked up once on mount;
  // null when the wrapper isn't tied to a job.
  const [projectInfo, setProjectInfo] = useState<{ id: string; name?: string } | null>(null);
  useEffect(() => {
    if (!jobId) return;
    let cancelled = false;
    (async () => {
      try {
        const jobInfo = await apiGet(`jobs/${jobId}`);
        if (cancelled || !jobInfo?.project) return;
        // jobInfo.project is the project pk; fetch the uuid + name.
        const proj = await apiGet(`projects/${jobInfo.project}`);
        if (cancelled || !proj?.uuid) return;
        setProjectInfo({ id: proj.uuid, name: proj.name });
      } catch (err) {
        console.warn("[wrapper] failed to resolve project for scene authoring:", err);
      }
    })();
    return () => { cancelled = true; };
  }, [jobId]);

  // Promote an editor YAML to a fully self-contained scene + asset
  // bundle: every URL-resolvable ref is fetched into assets/ and
  // rewritten to bundle: form; ligand dicts that live in Moorhen's
  // monomer library are re-collected (the lifter omits them so the
  // captured YAML stays small, but self-contained mode wants them).
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

  const handleCaptureScene = useCallback(async (): Promise<{
    scene: MoorhenScene;
    hints: SceneLiftHints;
    assets: SceneBundleAssets;
  }> => {
    const state = store.getState() as moorhen.State;
    const glRefState = (state as unknown as { glRef: {
      origin: number[] | Float32Array;
      quat: number[] | Float32Array;
      zoom: number;
      clipStart?: number;
      clipEnd?: number;
      fogStart?: number;
      fogEnd?: number;
    } }).glRef;
    // Flatten Moorhen's per-attribute contour slices into one
    // MapRenderState per molNo. mapContourSettings is keyed by
    // molNo across several parallel lists in the redux slice; we
    // gather them so the lifter can read off a single object.
    const mapState = collectMapRenderState(state);
    const activeMapMolNo =
      (state as unknown as { generalStates?: { activeMap?: moorhen.Map | null } })
        .generalStates?.activeMap?.molNo;
    return liftSceneToBundle({
      molecules,
      glRef: glRefState,
      projectId: projectInfo?.id,
      projectName: projectInfo?.name,
      // First molecule's monomerLibraryPath is the canonical Moorhen
      // root (all molecules share it via the wrapper's construction).
      // Without it the lifter falls back to STANDARD_MONOMERS only.
      monomerLibraryPath: molecules[0]?.monomerLibraryPath,
      maps,
      mapState,
      activeMapMolNo: typeof activeMapMolNo === "number" ? activeMapMolNo : undefined,
    });
  }, [store, molecules, maps, projectInfo]);

  // Custom side panels: CCP4i2 controls + Scenes (YAML editor + lifter +
  // apply). Scene operations live in the Scenes panel exclusively; the
  // CCP4i2 control row no longer carries the inline Apply/Download buttons.
  // Single Moorhen side-panel registration whose content is our own MUI
  // Tabs container hosting Controls and Scenes sub-panels. This sidesteps
  // an upstream Moorhen quirk where registering two extraSidePanels only
  // surfaces one tab in the panel switcher.
  const extraSidePanels: Record<string, MoorhenPanel> = useMemo(() => ({
    ccp4i2Controls: {
      icon: "MatSymSettings",
      label: "CCP4i2",
      panelContent: (
        <MoorhenCcp4i2TabbedPanel
          onFileSelect={fetchFile}
          onJobLoad={fetchJobFiles}
          getViewUrl={getViewUrl}
          molecules={molecules}
          maps={maps}
          onMapContourLevelChange={handleMapContourLevelChange}
          onRunServalcat={jobId ? handleRunServalcat : undefined}
          servalcatStatus={servalcatStatus}
          onApplyScene={handleApplyScene}
          onCaptureScene={handleCaptureScene}
          onPromoteSceneToPortable={handlePromoteSceneToPortable}
          cootInitialized={cootInitialized}
        />
      ),
    },
  }), [fetchFile, fetchJobFiles, getViewUrl, molecules, maps, handleMapContourLevelChange, jobId, handleRunServalcat, servalcatStatus, handleApplyScene, handleCaptureScene, handlePromoteSceneToPortable, cootInitialized]);

  // Moorhen 1.0 requires the InstanceProvider to be seeded with a menu system
  // (it builds the per-instance MoorhenInstance from it). One per wrapper.
  const menuSystem = useMemo(() => new MoorhenMenuSystem(), []);

  const collectedProps = useMemo(() => ({
    glRef,
    timeCapsuleRef,
    commandCentre,
    moleculesRef,
    mapsRef,
    activeMapRef: activeMapRef as React.RefObject<moorhen.Map>,
    lastHoveredAtom,
    prevActiveMoleculeRef,
    monomerLibraryPath,
    urlPrefix,
    store,
    viewOnly: false,
    extraSidePanels,
    setMoorhenDimensions,
  }), [urlPrefix, store, extraSidePanels, setMoorhenDimensions]);

  useEffect(() => {
    if (fileIds && cootInitialized) {
      fileIds.forEach((fileId) => {
        fetchFile(fileId);
      });
    }
  }, [fileIds, cootInitialized]);

  // Show Safari advisory as a non-blocking snackbar
  const isElectronEnv = typeof window !== "undefined" && !!(window as any).electronAPI;
  useEffect(() => {
    if (isSafari && !isElectronEnv) {
      setMessage(
        "Safari may have compatibility issues with the Moorhen viewer. For the best experience, consider using Chrome, Edge, or Firefox.",
        "warning"
      );
    }
  }, [isSafari, isElectronEnv, setMessage]);

  return (
    <div ref={moorhenContainerRef}>
      <MoorhenErrorBoundary
        fallback={
          <MoorhenFallback
            reason="runtime_error"
            capabilities={capabilities}
          />
        }
      >
        {store && (
          <MoorhenInstanceProvider menuSystem={menuSystem}>
            <MoorhenContainer {...collectedProps} />
          </MoorhenInstanceProvider>
        )}
      </MoorhenErrorBoundary>
    </div>
  );
};

export default MoorhenWrapper;

// --------------------------------------------------------------------------
// Helpers — capture-path
// --------------------------------------------------------------------------

interface ContourEntry {
  molNo: number;
  contourLevel?: number;
  radius?: number;
  alpha?: number;
  style?: "lines" | "solid" | "lit-lines";
  rgb?: { r: number; g: number; b: number };
}

/**
 * Flatten Moorhen's mapContourSettings slice — which stores per-attribute
 * arrays keyed by molNo — into a single MapRenderState per molNo for the
 * lifter to consume. Defensive: every field is optional, and we tolerate
 * the slice being absent (older Moorhen versions, or no maps loaded).
 */
function collectMapRenderState(
  state: unknown,
): Record<number, MapRenderState> {
  const out: Record<number, MapRenderState> = {};
  const slice = (state as { mapContourSettings?: Record<string, unknown> })
    .mapContourSettings;
  if (!slice) return out;
  const ensure = (molNo: number): MapRenderState => {
    if (!out[molNo]) out[molNo] = {};
    return out[molNo];
  };
  const merge = (
    rows: unknown,
    setter: (s: MapRenderState, v: ContourEntry) => void,
  ) => {
    if (!Array.isArray(rows)) return;
    for (const row of rows as ContourEntry[]) {
      if (typeof row?.molNo === "number") setter(ensure(row.molNo), row);
    }
  };
  merge(slice.contourLevels, (s, v) => { if (v.contourLevel !== undefined) s.contourLevel = v.contourLevel; });
  merge(slice.mapRadii, (s, v) => { if (v.radius !== undefined) s.radius = v.radius; });
  merge(slice.mapAlpha, (s, v) => { if (v.alpha !== undefined) s.alpha = v.alpha; });
  merge(slice.mapStyles, (s, v) => { if (v.style) s.style = v.style; });
  merge(slice.mapColours, (s, v) => { if (v.rgb) s.colour = rgb01ToHex(v.rgb); });
  merge(slice.positiveMapColours, (s, v) => { if (v.rgb) s.positiveColour = rgb01ToHex(v.rgb); });
  merge(slice.negativeMapColours, (s, v) => { if (v.rgb) s.negativeColour = rgb01ToHex(v.rgb); });
  const visible = (slice as { visibleMaps?: unknown }).visibleMaps;
  if (Array.isArray(visible)) {
    const visibleSet = new Set<number>(
      (visible as unknown[]).filter((n): n is number => typeof n === "number"),
    );
    // visibleMaps lists the *visible* molNos; map state already includes
    // every molNo we've seen, but make sure we don't drop a map that has
    // no other contour settings touched.
    for (const molNo of visibleSet) ensure(molNo);
    for (const molNo of Object.keys(out).map(Number)) {
      out[molNo].visible = visibleSet.has(molNo);
    }
  }
  return out;
}

/** Convert Moorhen's {r,g,b} 0-1 shape to a 6-hex string. */
function rgb01ToHex(rgb: { r: number; g: number; b: number }): string {
  const to8 = (v: number) =>
    Math.max(0, Math.min(255, Math.round(v * 255)))
      .toString(16)
      .padStart(2, "0");
  return `#${to8(rgb.r)}${to8(rgb.g)}${to8(rgb.b)}`;
}
