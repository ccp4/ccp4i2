"use client";

/**
 * Campaign-specific Moorhen wrapper with site navigation and member project switching.
 *
 * This is a specialized version of MoorhenWrapper for fragment screening campaigns.
 * It provides:
 * - Binding site navigation
 * - Member project switching via dropdown
 * - Standard view state URL support
 */

import {
  addMolecule,
  showMolecule,
  addMap,
  removeMolecule,
  removeMap,
  setActiveMap,
  setContourLevel,
  setTheme,
  setBackgroundColor,
  setOrigin,
  setQuat,
  setZoom,
  setRequestDrawScene,
  MoorhenContainer,
  MoorhenMolecule,
  MoorhenMap,
} from "moorhen/react-lib";
import { setShownSidePanel, MoorhenInstanceProvider, MoorhenMenuSystem } from "moorhen/react-lib";
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
import { apiText, apiArrayBuffer, apiGet, apiPost, apiUpload } from "../../api-fetch";
import { useTheme } from "../../theme/theme-provider";
import { useMoorhenViewState } from "../../hooks/use-moorhen-view-state";
import { useCampaignsApi } from "../../lib/campaigns-api";
import { usePopcorn } from "../../providers/popcorn-provider";
import {
  ProjectGroup,
  CampaignSite,
  MemberProjectWithSummary,
} from "../../types/campaigns";
import { Project } from "../../types/models";
import {
  MoorhenFallback,
  MoorhenErrorBoundary,
  useMoorhenCapabilities,
  isSafariBrowser,
} from "./moorhen-capability-check";
import {
  applyScene,
  SceneFileFetcher,
  SceneDictionaryFetcher,
  SceneDictionaryLoader,
  SceneMapFetcher,
  SceneResolveResult,
} from "../../lib/moorhen-scene-resolver";
import { parseScene, serialiseScene } from "../../lib/scene";
import { applyMaskDefaults, isMaskSubType, markMaskMap, ccp4Mode0ToFloat, ccp4DodgeEmClamp } from "../../lib/moorhen-map-file";
import type { MoorhenScene, SceneFileRef } from "../../types/moorhen-scene";
import { CampaignMoorhenTabbedPanel } from "./campaign-moorhen-tabbed-panel";
import type { SceneBundleAssets } from "./moorhen-scenes-panel";

type FileSource =
  | { type: "none" }
  | { type: "files"; fileIds: number[] }
  | { type: "job"; jobId: number };

export interface CampaignMoorhenWrapperProps {
  campaign: ProjectGroup;
  fileSource: FileSource;
  /** Campaign summary scene (built server-side). When set, it's seeded into
   *  the Scenes panel and auto-applied through that panel's own parse/apply
   *  path — the single rendering pathway shared with hand-edited scenes. */
  summaryScene?: MoorhenScene | null;
  viewParam?: string | null;
  sites: CampaignSite[];
  onUpdateSites: (sites: CampaignSite[]) => Promise<void>;
  memberProjects: MemberProjectWithSummary[];
  selectedMemberProjectId: number | null;
  onSelectMemberProject: (projectId: number | null) => void;
  parentProject: Project | null | undefined;
}

const CampaignMoorhenWrapper: React.FC<CampaignMoorhenWrapperProps> = ({
  campaign,
  fileSource,
  summaryScene,
  viewParam,
  sites,
  onUpdateSites,
  memberProjects,
  selectedMemberProjectId,
  onSelectMemberProject,
  parentProject,
}) => {
  const capabilities = useMoorhenCapabilities();
  const [isSafari] = useState(() => isSafariBrowser());
  const dispatch = useDispatch();
  const theme = useTheme();
  const campaignsApi = useCampaignsApi();
  const { setMessage } = usePopcorn();

  // Representation visibility state (lifted from control panel for URL capture)
  const [visibleRepresentations, setVisibleRepresentations] = useState<string[]>(["CRs"]);

  // Apply representations to all molecules (used when restoring from URL)
  const applyRepresentationsToMolecules = useCallback(
    async (newReps: string[], currentMolecules: moorhen.Molecule[]) => {
      // Default is ["CRs"], which is what molecules load with
      const defaultReps = ["CRs"];

      // Determine what to add/remove compared to default
      const toAdd = newReps.filter((r) => !defaultReps.includes(r));
      const toRemove = defaultReps.filter((r) => !newReps.includes(r));

      for (const mol of currentMolecules) {
        for (const rep of toAdd) {
          try {
            await mol.addRepresentation(rep, "/*/*/*/*");
          } catch (err) {
            console.error(`Failed to add ${rep} to ${mol.name}:`, err);
          }
        }
        for (const rep of toRemove) {
          try {
            mol.clearBuffersOfStyle(rep);
          } catch (err) {
            console.error(`Failed to remove ${rep} from ${mol.name}:`, err);
          }
        }
      }
    },
    []
  );

  // View state hook for URL parameter support
  const { getViewUrl, initialRepresentations } = useMoorhenViewState({
    viewParam: viewParam ?? null,
    onViewRestored: () => {},
    representations: visibleRepresentations,
  });

  const hasInitializedReps = useRef(false);

  // Container ref for measuring available height below the AppBar.
  // Moorhen 1.0 no longer accepts a setMoorhenDimensions callback; it takes a
  // static `size` prop and otherwise defaults to the full window.innerHeight,
  // which overflows the viewport because the viewer renders under a toolbar.
  // So we measure the container's offset and feed Moorhen an explicit size,
  // keeping it in sync on window resize.
  const moorhenContainerRef = useRef<HTMLDivElement>(null);
  const [size, setSize] = useState<[number, number]>(() =>
    typeof window === "undefined"
      ? [0, 0]
      : [window.innerWidth, window.innerHeight]
  );
  useEffect(() => {
    const measure = () => {
      const top = moorhenContainerRef.current?.getBoundingClientRect().top ?? 0;
      setSize([window.innerWidth, window.innerHeight - top]);
    };
    measure();
    window.addEventListener("resize", measure);
    return () => window.removeEventListener("resize", measure);
  }, []);

  const glRef: RefObject<webGL.MGWebGL | null> = useRef(null);
  const commandCentre = useRef<null | moorhen.CommandCentre>(null);
  const moleculesRef = useRef<null | moorhen.Molecule[]>(null);
  const mapsRef = useRef<null | moorhen.Map[]>(null);
  const activeMapRef = useRef<moorhen.Map>(null);
  const lastHoveredAtom = useRef<null | moorhen.HoveredAtom>(null);
  const prevActiveMoleculeRef = useRef<null | moorhen.Molecule>(null);
  const timeCapsuleRef = useRef(null);
  const loadedFileSource = useRef<FileSource | null>(null);

  // Ligand dictionary file ID for 2D structure display (first dict file found)
  const [ligandDictFileId, setLigandDictFileId] = useState<number | null>(null);
  const [ligandName, setLigandName] = useState<string | null>(null);
  // Store ALL loaded dictionary contents so we can add them to molecules
  const loadedDictContents = useRef<string[]>([]);

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

  // Initialize representations from URL if available (after molecules load)
  useEffect(() => {
    if (!hasInitializedReps.current && initialRepresentations && molecules.length > 0) {
      setVisibleRepresentations(initialRepresentations);
      applyRepresentationsToMolecules(initialRepresentations, molecules);
      hasInitializedReps.current = true;
    }
  }, [initialRepresentations, molecules, applyRepresentationsToMolecules]);

  useEffect(() => {
    dispatch(
      setBackgroundColor(theme.mode === "light" ? [1, 1, 1, 1] : [0, 0, 0, 1])
    );
    dispatch(setTheme(theme.mode === "light" ? "flatly" : "darkly"));
  }, [theme.mode, dispatch]);

  // Auto-open the campaign controls side panel
  useEffect(() => {
    dispatch(setShownSidePanel("campaignControls"));
  }, [dispatch]);

  const monomerLibraryPath =
    "https://raw.githubusercontent.com/MonomerLibrary/monomers/master/";

  const backgroundColor = useSelector(
    (state: moorhen.State) => state.sceneSettings.backgroundColor
  );
  const defaultBondSmoothness = useSelector(
    (state: moorhen.State) => state.sceneSettings.defaultBondSmoothness
  );

  const isElectron =
    typeof window !== "undefined" && !!(window as any).electronAPI;
  // In web browsers, use API route for CORP headers (COEP compatibility)
  // In Electron, serve directly from public/MoorhenAssets
  const urlPrefix = isElectron ? "/MoorhenAssets" : "/api/moorhen/MoorhenAssets";

  const getOrigin = useCallback(() => {
    const state = store.getState() as moorhen.State;
    return (state as unknown as { glRef: { origin: number[] } }).glRef.origin;
  }, [store]);

  // Cleanup all loaded molecules and maps
  const cleanupLoadedContent = useCallback(() => {
    // Remove all molecules
    for (const mol of molecules) {
      dispatch(removeMolecule(mol));
      mol.delete();
    }
    // Remove all maps
    for (const map of maps) {
      dispatch(removeMap(map));
      map.delete();
    }
    // Trigger redraw
    dispatch(setRequestDrawScene(true));
    // Reset ligand info
    setLigandDictFileId(null);
    setLigandName(null);
    loadedDictContents.current = [];
    // Reset representation state to default
    setVisibleRepresentations(["CRs"]);
    hasInitializedReps.current = false;
  }, [molecules, maps, dispatch]);

  // ----------------------------------------------------------------------
  // Scene apply path (used by the campaign "Summary View")
  //
  // A summary scene (built server-side by lib/campaign_scene.py) overlays
  // every discovered fragment on the parent ribbon, each fragment scoped to
  // its own restraint dictionary. We drive it through the shared scene
  // resolver (applyScene) rather than the ad-hoc job loader above, because
  // the resolver already does the load-global-then-scope-per-molecule
  // dictionary dance the fragment-campaign case needs.
  // ----------------------------------------------------------------------

  // Bundle assets (from an opened .scene.zip in the Scenes panel) live in a
  // ref so the stable-identity fetchers below can reach this apply's assets;
  // refreshed at the top of each apply.
  const bundleAssetsRef = useRef<SceneBundleAssets>(new Map());

  // Structure load for the scene fetcher. Mirrors the generic wrapper's
  // proven loadStructureFromText: load coords, add default reps, centre, and
  // crucially dispatch(showMolecule) so the molecule is actually visible. The
  // resolver then hides these loader-default reps and applies the scene's own
  // (so "scene owns the look" still holds) — but without showMolecule the
  // molecule loads invisibly and the viewer stays blank.
  const loadSceneStructureFromText = useCallback(
    async (
      coordText: string,
      molName: string,
      uniqueId: string,
    ): Promise<moorhen.Molecule | null> => {
      if (!commandCentre.current) return null;
      const newMolecule = new MoorhenMolecule(
        commandCentre as RefObject<moorhen.CommandCentre>,
        store as any,
        monomerLibraryPath,
      );
      newMolecule.setBackgroundColour(backgroundColor);
      newMolecule.defaultBondOptions.smoothness = defaultBondSmoothness;
      try {
        await newMolecule.loadToCootFromString(coordText, molName);
        if (newMolecule.molNo === -1) throw new Error("Cannot read coordinates");
        newMolecule.uniqueId = uniqueId;
        // Ribbon first (protein overview), fall back to sticks. These get
        // hidden by the resolver and replaced with the scene's reps.
        try {
          await newMolecule.addRepresentation("CRs", "/*/*/*/*");
        } catch {
          await newMolecule.addRepresentation("CBs", "/*/*/*/*");
        }
        try {
          await newMolecule.addRepresentation("ligands", "/*/*/*/*");
        } catch {
          /* no ligands present — fine */
        }
        await newMolecule.centreOn("/*/*/*/*", false, true);
        dispatch(addMolecule(newMolecule));
        dispatch(showMolecule({ molNo: newMolecule.molNo } as never));
        return newMolecule;
      } catch (err) {
        console.warn(`[scene] failed to load ${molName} (${uniqueId}):`, err);
        return null;
      }
    },
    [commandCentre, store, monomerLibraryPath, backgroundColor, defaultBondSmoothness, dispatch],
  );

  const loadSceneStructure = useCallback(
    async (url: string, molName: string): Promise<moorhen.Molecule | null> => {
      try {
        const pdbData = await apiText(url);
        return loadSceneStructureFromText(pdbData, molName, url);
      } catch (err) {
        console.warn(`[scene] failed to fetch ${molName} from ${url}:`, err);
        return null;
      }
    },
    [loadSceneStructureFromText],
  );

  const handleFetchSceneFile: SceneFileFetcher = useCallback(
    async (ref: SceneFileRef) => {
      // Bundle: decode bytes from the in-memory asset map (no network).
      if (ref.bundle) {
        const buf = bundleAssetsRef.current.get(ref.bundle);
        if (!buf) {
          console.warn(`[scene] bundle lookup miss: ${ref.bundle}`);
          return null;
        }
        const coordText = new TextDecoder("utf-8").decode(buf);
        return loadSceneStructureFromText(
          coordText,
          ref.name || ref.bundle,
          `bundle:${ref.bundle}`,
        );
      }
      if (ref.fileId !== undefined && ref.projectId) {
        return loadSceneStructure(
          `/api/proxy/ccp4i2/files/${ref.fileId}/download/`,
          ref.name || `file_${ref.fileId}`,
        );
      }
      if (ref.pdb) {
        const pdbId = ref.pdb.toLowerCase();
        return loadSceneStructure(
          `/api/proxy/pdbe/entry-files/download/${pdbId}.cif`,
          ref.name || pdbId,
        );
      }
      if (ref.url) return loadSceneStructure(ref.url, ref.name || ref.url);
      return null;
    },
    [loadSceneStructure, loadSceneStructureFromText],
  );

  const handleFetchSceneDictionary: SceneDictionaryFetcher = useCallback(
    async (ref: SceneFileRef): Promise<string | null> => {
      if (ref.cifText) return ref.cifText;
      if (ref.bundle) {
        const buf = bundleAssetsRef.current.get(ref.bundle);
        if (!buf) return null;
        return new TextDecoder("utf-8").decode(buf);
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
    [],
  );

  // Map fetcher: loads a scene maps[] entry. Real-space CCP4 map files
  // (kind: "map", incl. masks) load via loadToCootFromMapData; MTZ refs via
  // loadToCootFromMtzData with the column spec. applyMapState (in the resolver)
  // then applies contour/style/colour, incl. the mask defaults.
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
        let url: string | null = null;
        if (ref.fileId !== undefined && ref.projectId) {
          url = `/api/proxy/ccp4i2/files/${ref.fileId}/download/`;
        } else if (ref.url) {
          url = ref.url;
        }
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
        if (ref.kind === "map") {
          // mode-0 -> float (sane stats); masks also dodge coot's EM cell-clamp.
          const mapBytes = ccp4Mode0ToFloat(bytes as ArrayBuffer);
          await newMap.loadToCootFromMapData(
            new Uint8Array(sceneMap.isMask ? ccp4DodgeEmClamp(mapBytes) : mapBytes),
            sceneMap.name,
            !!sceneMap.isDifference,
          );
          (newMap as any).isCcp4MapFile = true;
          if (sceneMap.isMask) markMaskMap(newMap);
        } else {
          const cols = sceneMap.columns ?? {};
          await newMap.loadToCootFromMtzData(
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
          );
        }
        if (newMap.molNo === -1) return null;
        if (uniqueId) newMap.uniqueId = uniqueId;
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
    [commandCentre, store, dispatch],
  );

  // Core apply: hand a parsed scene + bundle assets to the resolver.
  const runScene = useCallback(
    async (
      scene: MoorhenScene,
      assets: SceneBundleAssets = new Map(),
    ): Promise<SceneResolveResult> => {
      bundleAssetsRef.current = assets;
      const result = await applyScene({
        scene,
        molecules,
        maps,
        dispatch,
        fetcher: handleFetchSceneFile,
        dictionaryFetcher: handleFetchSceneDictionary,
        dictionaryLoader: handleLoadSceneDictionary,
        mapFetcher: handleFetchSceneMap,
      });
      dispatch(setRequestDrawScene(true));
      return result;
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

  // Scenes-panel path: parse the editor YAML and apply (with bundle assets
  // when opened from a .scene.zip). This is the single apply entry — the
  // summary view drives it too, via the panel's auto-apply.
  const handleApplyScene = useCallback(
    (yamlText: string, assets: SceneBundleAssets): Promise<SceneResolveResult> => {
      return runScene(parseScene(yamlText), assets);
    },
    [runScene],
  );

  // Load files when fileSource changes
  useEffect(() => {
    if (!cootInitialized) return;

    // Check if file source changed
    const sourceChanged =
      JSON.stringify(loadedFileSource.current) !== JSON.stringify(fileSource);
    if (!sourceChanged) return;

    // Cleanup existing content before loading new files
    if (loadedFileSource.current !== null) {
      cleanupLoadedContent();
    }

    loadedFileSource.current = fileSource;

    if (fileSource.type === "files" && fileSource.fileIds.length > 0) {
      fileSource.fileIds.forEach((fileId) => {
        fetchFile(fileId);
      });
    } else if (fileSource.type === "job") {
      fetchJobFiles(fileSource.jobId);
    }
    // The summary scene is applied by the Scenes panel (auto-apply), not here,
    // so it shares one rendering pathway with hand-edited scenes.
  }, [fileSource, cootInitialized, cleanupLoadedContent]);

  // Dimension updates are handled by Moorhen's MainContainer automatically

  const fetchFile = async (fileId: number) => {
    const fileInfo = await apiGet(`files/${fileId}`);
    if (!fileInfo) {
      console.warn(`File with ID ${fileId} not found.`);
      return;
    }
    if (fileInfo.type === "chemical/x-pdb") {
      const url = `/api/proxy/ccp4i2/files/${fileId}/download/`;
      const molName = fileInfo.annotation || fileInfo.job_param_name;
      await fetchMolecule(url, molName);
    } else if (fileInfo.type === "application/CCP4-mtz-map") {
      const url = `/api/proxy/ccp4i2/files/${fileId}/download/`;
      const molName = fileInfo.name || fileInfo.job_param_name;
      // subType: 1=normal, 2=difference, 3=anomalous difference
      const mapSubType = fileInfo.sub_type || 1;
      await fetchMap(url, molName, mapSubType);
    } else if (fileInfo.type === "application/CCP4-map") {
      const url = `/api/proxy/ccp4i2/files/${fileId}/download/`;
      const molName = fileInfo.annotation || fileInfo.name || fileInfo.job_param_name;
      await fetchMapFile(url, molName, { isMask: isMaskSubType(fileInfo.sub_type) });
    }
  };

  const fetchJobFiles = async (jobId: number) => {
    const files = await apiGet(`files/?job=${jobId}`);
    if (!files || !Array.isArray(files)) return;

    // Filter to only JOB_DIR files (directory=1), exclude imported files (directory=2)
    const jobOutputFiles = files.filter((f: { directory: number }) => f.directory === 1);

    // STEP 1: Load ALL ligand dictionaries FIRST (before coordinates)
    // This ensures coot understands ligand geometry when parsing coordinates.
    // A dictionary file may contain multiple monomers — read_dictionary_string
    // loads all of them into coot's global store.
    const ligandDictFiles = files.filter(
      (f: { type: string }) => f.type === "application/refmac-dictionary"
    );
    if (ligandDictFiles.length > 0) {
      loadedDictContents.current = [];
      for (const dictFile of ligandDictFiles) {
        const dictUrl = `/api/proxy/ccp4i2/files/${dictFile.id}/download/`;
        await fetchDict(dictUrl);
      }
      // Use the first dictionary file for 2D display in the control panel
      const firstDict = ligandDictFiles[0];
      setLigandDictFileId(firstDict.id);
      const name = firstDict.name?.replace(/\.cif$/i, "") ||
                   firstDict.annotation ||
                   "Ligand";
      setLigandName(name);
    } else {
      loadedDictContents.current = [];
      setLigandDictFileId(null);
      setLigandName(null);
    }

    // STEP 2: Find and load coordinate files
    // Check for both PDB and mmCIF types
    const coordFiles = jobOutputFiles.filter(
      (f: { type: string }) =>
        f.type === "chemical/x-pdb" ||
        f.type === "chemical/x-cif" ||
        f.type === "chemical/x-mmcif"
    );
    // Prefer mmCIF (.cif) over PDB (.pdb) for coordinates
    const mmcifFile = coordFiles.find((f: { name: string }) =>
      f.name.toLowerCase().endsWith(".cif")
    );
    const coordFile = mmcifFile || coordFiles[0];

    // Load the single best coordinate file
    if (coordFile) {
      const url = `/api/proxy/ccp4i2/files/${coordFile.id}/download/`;
      const molName = coordFile.annotation || coordFile.job_param_name;
      await fetchMolecule(url, molName);
    }

    // STEP 3: Load map files (MTZ coefficients and real-space CCP4 maps / masks)
    for (const file of jobOutputFiles) {
      if (file.type === "application/CCP4-mtz-map") {
        const url = `/api/proxy/ccp4i2/files/${file.id}/download/`;
        const molName = file.name || file.job_param_name;
        // subType: 1=normal, 2=difference, 3=anomalous difference
        const mapSubType = file.sub_type || 1;
        await fetchMap(url, molName, mapSubType);
      } else if (file.type === "application/CCP4-map") {
        const url = `/api/proxy/ccp4i2/files/${file.id}/download/`;
        const molName = file.annotation || file.name || file.job_param_name;
        await fetchMapFile(url, molName, { isMask: isMaskSubType(file.sub_type) });
      }
    }
  };

  /**
   * Load a ligand dictionary into coot's global dictionary store.
   * This should be called BEFORE loading coordinates so coot understands ligand geometry.
   */
  const fetchDict = async (url: string): Promise<string | null> => {
    if (!commandCentre.current) return null;
    try {
      const fileContent = await apiText(url);
      // Load dictionary globally into coot (molNo=-999999 means global)
      await commandCentre.current.cootCommand(
        {
          returnType: "status",
          command: "read_dictionary_string",
          commandArgs: [fileContent, -999999],
          changesMolecules: [],
        },
        false
      );
      // Store content so we can add it to molecules later
      loadedDictContents.current.push(fileContent);
      return fileContent;
    } catch (err) {
      console.error("[fetchDict] Failed to load dictionary:", err);
      return null;
    }
  };

  const fetchMolecule = async (url: string, molName: string) => {
    if (!commandCentre.current) return;
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
      if (newMolecule.molNo === -1) {
        throw new Error("Cannot read the fetched molecule...");
      }
      newMolecule.uniqueId = url;

      // Add all loaded dictionaries to molecule
      // This ensures the molecule understands geometry for all monomers
      for (const dictContent of loadedDictContents.current) {
        try {
          await newMolecule.addDict(dictContent);
        } catch (err) {
          console.warn("[fetchMolecule] Failed to add dictionary:", err);
        }
      }
      // Try ribbon representation first (better for protein overview)
      // Fall back to CBs if ribbons fail (e.g., no protein backbone)
      try {
        await newMolecule.addRepresentation("CRs", "/*/*/*/*");
      } catch {
        // Ribbons failed, fall back to CBs
        await newMolecule.addRepresentation("CBs", "/*/*/*/*");
      }

      // Always try to add ligand representation
      try {
        await newMolecule.addRepresentation("ligands", "/*/*/*/*");
      } catch {
        console.warn("[fetchMolecule] Ligands representation failed");
      }

      await newMolecule.centreOn("/*/*/*/*", false, true);
      dispatch(addMolecule(newMolecule));
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch PDB entry from ${url}`);
    }
  };

  const fetchMap = async (
    url: string,
    mapName: string,
    mapSubType: number = 1
  ) => {
    if (!commandCentre.current) return;
    const newMap = new MoorhenMap(
      commandCentre as RefObject<moorhen.CommandCentre>,
      store as any
    );
    // subType: 1=normal, 2=difference, 3=anomalous difference
    // Both difference and anomalous maps use isDifference=true for contouring
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
      // Store the original sub_type for proper labeling and coloring
      (newMap as any).mapSubType = mapSubType;
      // Set custom colors for anomalous maps (orange/purple instead of green/red)
      if (mapSubType === 3) {
        newMap.defaultPositiveMapColour = { r: 1.0, g: 0.65, b: 0.0 }; // Orange
        newMap.defaultNegativeMapColour = { r: 0.6, g: 0.3, b: 0.8 }; // Purple
      }
      if (newMap.molNo === -1) throw new Error("Cannot read the fetched map...");
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
  };

  // Load a real-space CCP4 map file (application/CCP4-map), e.g. a mask. Uses
  // loadToCootFromMapData (not the MTZ path); masks get the shared translucent
  // solid defaults and are never the active map.
  const fetchMapFile = async (
    url: string,
    mapName: string,
    opts: { isMask?: boolean } = {}
  ) => {
    if (!commandCentre.current) return;
    const newMap = new MoorhenMap(
      commandCentre as RefObject<moorhen.CommandCentre>,
      store as any
    );
    try {
      // Convert mode-0 (int8) CCP4 maps to float so Moorhen reads sane stats
      // (no-op if already float). For masks, also nudge the P1/orthogonal cell
      // off 90° so coot contours periodically instead of clamping to the cell box.
      let mapData = ccp4Mode0ToFloat(await apiArrayBuffer(url));
      if (opts.isMask) mapData = ccp4DodgeEmClamp(mapData);
      await newMap.loadToCootFromMapData(new Uint8Array(mapData), mapName, false);
      if (newMap.molNo === -1) throw new Error("Cannot read the fetched map file...");
      newMap.uniqueId = url;
      // Tag so the lifter captures it as a kind: "map" ref (not MTZ).
      (newMap as any).isCcp4MapFile = true;
      if (opts.isMask) {
        markMaskMap(newMap);
      }
      dispatch(addMap(newMap));
      if (opts.isMask) {
        await applyMaskDefaults(dispatch, newMap as any);
      }
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch map file from ${url}`);
    }
  };

  // Navigate to a site
  const handleGoToSite = useCallback(
    (site: CampaignSite) => {
      dispatch(setOrigin(site.origin));
      if (site.quat) {
        dispatch(setQuat(site.quat));
      }
      if (site.zoom) {
        dispatch(setZoom(site.zoom));
      }
      dispatch(setRequestDrawScene(true));
    },
    [dispatch]
  );

  // Save current view as a site
  const handleSaveCurrentAsSite = useCallback(
    async (name: string) => {
      const state = store.getState() as unknown as {
        glRef: { origin: number[]; quat: number[]; zoom: number };
      };
      const newSite: CampaignSite = {
        name,
        origin: Array.from(state.glRef.origin).slice(0, 3) as [
          number,
          number,
          number
        ],
        quat: Array.from(state.glRef.quat).slice(0, 4) as [
          number,
          number,
          number,
          number
        ],
        zoom: state.glRef.zoom,
      };
      await onUpdateSites([...sites, newSite]);
    },
    [store, sites, onUpdateSites]
  );

  // Delete a site
  const handleDeleteSite = useCallback(
    async (index: number) => {
      const newSites = sites.filter((_, i) => i !== index);
      await onUpdateSites(newSites);
    },
    [sites, onUpdateSites]
  );

  // Update a site (rename and optionally update position)
  const handleUpdateSite = useCallback(
    async (index: number, name: string, updatePosition: boolean) => {
      const existingSite = sites[index];
      let updatedSite: CampaignSite;

      if (updatePosition) {
        // Capture current view position
        const state = store.getState() as unknown as {
          glRef: { origin: number[]; quat: number[]; zoom: number };
        };
        updatedSite = {
          name,
          origin: Array.from(state.glRef.origin).slice(0, 3) as [
            number,
            number,
            number
          ],
          quat: Array.from(state.glRef.quat).slice(0, 4) as [
            number,
            number,
            number,
            number
          ],
          zoom: state.glRef.zoom,
        };
      } else {
        // Keep existing position, just update name
        updatedSite = {
          ...existingSite,
          name,
        };
      }

      const newSites = [...sites];
      newSites[index] = updatedSite;
      await onUpdateSites(newSites);
    },
    [store, sites, onUpdateSites]
  );

  // Handle map contour level changes
  // NOTE: This must be synchronous for slider dragging to work smoothly.
  // The map redraw is fired off without awaiting to prevent blocking.
  const handleMapContourLevelChange = useCallback(
    (molNo: number, level: number) => {
      // Update Redux state - Moorhen's getMapContourParams() reads contourLevel from here
      // Note: Despite TypeScript def saying 'level', Moorhen internally expects 'contourLevel'
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      dispatch(setContourLevel({ molNo, contourLevel: level } as any));

      // Find the map and redraw its contour (fire-and-forget, don't block)
      const map = maps.find((m) => m.molNo === molNo);
      if (map) {
        map.drawMapContour().catch((err) => {
          console.error("Failed to redraw map contour:", err);
        });
      }

      // Trigger scene redraw
      dispatch(setRequestDrawScene(true));
    },
    [dispatch, maps]
  );

  // Handle tagging the currently selected project with a site name
  const handleTagProjectWithSite = useCallback(
    async (siteName: string) => {
      if (!selectedMemberProjectId) {
        setMessage("Please select a member project first");
        return;
      }
      try {
        await campaignsApi.addTagByText(selectedMemberProjectId, siteName);
        const memberProject = memberProjects.find((p) => p.id === selectedMemberProjectId);
        const projectName = memberProject?.name || `Project ${selectedMemberProjectId}`;
        setMessage(`Tagged "${projectName}" with "${siteName}"`);
      } catch (err) {
        console.error("Failed to tag project:", err);
        setMessage("Failed to tag project");
      }
    },
    [selectedMemberProjectId, campaignsApi, setMessage, memberProjects]
  );

  // Run servalcat_pipe refinement on a molecule
  const handleRunServalcat = useCallback(
    async (mol: moorhen.Molecule) => {
      // Refinement always operates on a member project (child) - the reflections
      // and coordinates belong to the specific crystal, not the campaign parent
      if (!selectedMemberProjectId) {
        setMessage("No member project selected");
        return;
      }
      const projectId = selectedMemberProjectId;

      const mp = memberProjects.find((p) => p.id === selectedMemberProjectId);
      if (!mp?.uuid) {
        setMessage("Cannot determine project UUID");
        return;
      }
      const projectUuid = mp.uuid;
      const projectDbId = projectUuid.replace(/-/g, "");

      setMessage("Creating servalcat refinement job...");

      try {
        // Step 1: Find observation reflections from this project's top-level jobs
        type ObsFile = { id: number; job: number; uuid: string; name: string;
                         content: number | null; sub_type: number | null };
        const obsFiles = await apiGet(
          `files/?type=application/CCP4-mtz-observed&directory=1` +
          `&job__project=${projectId}&job__parent__isnull=true`
        ) as ObsFile[] | null;

        if (!obsFiles || obsFiles.length === 0) {
          setMessage("No reflection data found in project");
          return;
        }

        // Sort by job ID descending (most recent first), then prefer
        // anomalous data (IPAIR content & 1, FPAIR content & 2) over IMEAN/FMEAN
        const sorted = [...obsFiles].sort((a, b) => b.job - a.job);
        const hasAnomalous = (f: ObsFile) => f.content !== null && (f.content & 3) !== 0;
        const reflectionFile = sorted.find(hasAnomalous) || sorted[0];

        // Step 2: Create servalcat_pipe job
        const jobResponse = await apiPost<{
          status: string;
          data: { new_job: { id: number; uuid: string } };
        }>(`projects/${projectId}/create_task/`, {
          task_name: "servalcat_pipe",
          title: `Servalcat refinement of ${mol.name}`,
        });
        const newJobId = jobResponse.data.new_job.id;

        // Step 3: Upload coordinates from Moorhen (preserve original format)
        setMessage("Uploading coordinates...");
        const coordText = await mol.getAtoms();
        const ismmCIF = (mol as any).coordsFormat === "mmcif";
        const mimeType = ismmCIF ? "chemical/x-cif" : "chemical/x-pdb";
        const ext = ismmCIF ? ".cif" : ".pdb";
        const coordBlob = new Blob([coordText], { type: mimeType });
        const coordFile = new File([coordBlob], `${mol.name || "coords"}${ext}`);
        const coordFormData = new FormData();
        coordFormData.append("file", coordFile);
        coordFormData.append("object_path", "servalcat_pipe.inputData.XYZIN");
        await apiUpload(`jobs/${newJobId}/upload_file_param/`, coordFormData);

        // Step 4: Link reflections via database file reference
        setMessage("Linking reflection data...");
        const reflDbFileId = reflectionFile.uuid.replace(/-/g, "");
        const hklinValue: Record<string, string | number> = {
          project: projectDbId,
          dbFileId: reflDbFileId,
        };
        if (reflectionFile.sub_type !== null) {
          hklinValue.subType = reflectionFile.sub_type;
        }
        await apiPost(`jobs/${newJobId}/set_parameter/`, {
          object_path: "servalcat_pipe.inputData.HKLIN",
          value: hklinValue,
        });

        // Step 5: Upload dictionary if available
        // Uses upload_file_param (not set_parameter) because DICT_LIST starts empty
        // and upload_file_param handles list expansion automatically
        if (ligandDictFileId) {
          setMessage("Uploading ligand dictionary...");
          const dictUrl = `/api/proxy/ccp4i2/files/${ligandDictFileId}/download/`;
          const dictText = await apiText(dictUrl);
          const dictBlob = new Blob([dictText], { type: "application/refmac-dictionary" });
          const dictFile = new File([dictBlob], "ligand.cif");
          const dictFormData = new FormData();
          dictFormData.append("file", dictFile);
          dictFormData.append("object_path", "servalcat_pipe.inputData.DICT_LIST[0]");
          await apiUpload(`jobs/${newJobId}/upload_file_param/`, dictFormData);
        }

        // Step 6: Run the job
        setMessage("Running servalcat refinement...");
        await apiPost(`jobs/${newJobId}/run/`, {});

        setMessage("Servalcat refinement job submitted successfully");
      } catch (error) {
        console.error("Failed to create servalcat job:", error);
        setMessage(
          `Failed to create servalcat job: ${error instanceof Error ? error.message : "Unknown error"}`
        );
      }
    },
    [selectedMemberProjectId, memberProjects, ligandDictFileId, setMessage]
  );

  // Moorhen 1.0 requires the InstanceProvider to be seeded with a menu system
  // (it builds the per-instance MoorhenInstance from it).
  const menuSystem = useMemo(() => new MoorhenMenuSystem(), []);

  // When viewing the campaign summary, serialise the scene to YAML so the
  // Scenes panel can show it in the editor and auto-apply it through its own
  // parse/apply path (one shared rendering pathway).
  const summarySceneYaml = useMemo(() => {
    if (!summaryScene) return undefined;
    try {
      return serialiseScene(summaryScene);
    } catch (err) {
      console.warn("[scene] failed to serialise summary scene for editor:", err);
      return undefined;
    }
  }, [summaryScene]);

  // Custom side panel: campaign Controls + Scenes editor under one Moorhen
  // side-panel registration (matching the job/file viewers' tabbed panel).
  const extraSidePanels: Record<string, MoorhenPanel> = {
    campaignControls: {
      icon: "MatSymSettings",
      label: "Campaign",
      panelContent: (
        <CampaignMoorhenTabbedPanel
          controlPanelProps={{
            campaign,
            sites,
            onGoToSite: handleGoToSite,
            onSaveCurrentAsSite: handleSaveCurrentAsSite,
            onUpdateSite: handleUpdateSite,
            onDeleteSite: handleDeleteSite,
            memberProjects,
            selectedMemberProjectId,
            onSelectMemberProject,
            parentProject,
            getViewUrl,
            molecules,
            visibleRepresentations,
            onRepresentationsChange: setVisibleRepresentations,
            ligandDictFileId,
            ligandName,
            maps,
            onMapContourLevelChange: handleMapContourLevelChange,
            onTagProjectWithSite: handleTagProjectWithSite,
            onFileSelect: fetchFile,
            onRunServalcat: handleRunServalcat,
          }}
          onApplyScene={handleApplyScene}
          cootInitialized={cootInitialized}
          initialSceneYaml={summarySceneYaml}
          autoApplyInitialScene={!!summarySceneYaml}
        />
      ),
    },
  };

  const collectedProps = {
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
    size,
  };

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
      <MoorhenErrorBoundary fallback={<MoorhenFallback reason="runtime_error" capabilities={capabilities} />}>
        {store && (
          <MoorhenInstanceProvider menuSystem={menuSystem}>
            <MoorhenContainer {...collectedProps} />
          </MoorhenInstanceProvider>
        )}
      </MoorhenErrorBoundary>
    </div>
  );
};

export default CampaignMoorhenWrapper;
