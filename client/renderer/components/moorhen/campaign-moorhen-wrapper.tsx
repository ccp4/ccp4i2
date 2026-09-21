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
import type { MoorhenInstance, MoorhenPanel } from "moorhen/react-lib";

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
import { readCameraState } from "../../lib/moorhen-view-state";
import { useCampaignsApi } from "../../lib/campaigns-api";
import { usePopcorn } from "../../providers/popcorn-provider";
import {
  ProjectGroup,
  CampaignSite,
  NewCampaignSite,
  SiteEvaluation,
  SiteVerdict,
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
import { applyMaskDefaults, isMaskSubType, markMaskMap, ccp4Mode0ToFloat, ccp4DodgeEmClamp, requireMoorhenInstance, primeEmMapHeaderInfo } from "../../lib/moorhen-map-file";
import {
  COORDINATE_TYPES,
  fetchCompanionDictionaryFiles,
  fetchDictionaryTexts,
  fetchJobDictionaryFiles,
  loadWithDictionaries,
  type DictionaryToAttach,
} from "../../lib/moorhen-dictionaries";
import { candidateLigandCodes, placeLigand } from "../../lib/ligand-codes";
import type { MoorhenScene, SceneFileRef } from "../../types/moorhen-scene";
import { isElectronWindow, moorhenUrlPrefix } from "../../lib/moorhen-asset-path";
import { prefetchMoorhenWasm } from "../../lib/moorhen-wasm-prefetch";
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
  /** A site to move to once the scene is up — the `site` URL parameter, which
   *  is how a verdict chip in the campaign overview opens its site. */
  initialSiteId?: number | null;
  sites: CampaignSite[];
  onAddSite: (site: NewCampaignSite) => Promise<void>;
  onUpdateSite: (
    siteId: number,
    changes: Partial<NewCampaignSite>
  ) => Promise<void>;
  onDeleteSite: (siteId: number) => Promise<void>;
  /** Verdicts recorded for the selected dataset, empties included. */
  evaluations: SiteEvaluation[];
  onSetVerdict: (siteId: number, verdict: SiteVerdict) => Promise<void>;
  onClearVerdict: (siteId: number) => Promise<void>;
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
  initialSiteId,
  sites,
  onAddSite,
  onUpdateSite,
  onDeleteSite,
  evaluations,
  onSetVerdict,
  onClearVerdict,
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
            // reps arrive as string[] from URL-restore state; they are the
            // ToggleButton representation vocabulary, so narrow to the union.
            await mol.addRepresentation(rep as moorhen.RepresentationStyles, "/*/*/*/*");
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
  // Filled by MoorhenContainer on mount. Molecules and maps are built from it.
  const moorhenInstanceRef = useRef<null | MoorhenInstance>(null);
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
  // What "Add ligand here" needs: the codes the member's dictionaries define
  // and the file its coordinates came from (the molecule is found by that,
  // not by which molecule is active; see lib/ligand-codes).
  const [ligandCodes, setLigandCodes] = useState<string[]>([]);
  const [memberCoordFileId, setMemberCoordFileId] = useState<number | null>(null);
  // Store ALL loaded dictionary contents so we can add them to molecules

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

  // In web browsers, use API route for CORP headers (COEP compatibility), with
  // the Moorhen version in the path so the route's immutable cache is honest
  // across upgrades. In Electron, serve directly from public/MoorhenAssets.
  const isElectron = isElectronWindow();
  const urlPrefix = moorhenUrlPrefix(isElectron);

  // Start the WASM download now, rather than after the data archives.
  //
  // MoorhenCommandCentre.init() awaits three .tar.gz fetches before posting
  // CootInitialize, and only that message makes the worker fetch the WASM --
  // the biggest file on the page, last in the queue. Warming the cache here
  // lets the two run together; the worker's own fetch then joins it. No effect
  // in Electron, which reads these files from disk.
  useEffect(() => {
    if (isElectron) return;
    return prefetchMoorhenWasm(urlPrefix);
  }, [isElectron, urlPrefix]);

  const getOrigin = useCallback(() => {
    return readCameraState(store.getState() as moorhen.State).origin;
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
      opts: { dictionaries?: DictionaryToAttach[] } = {},
    ): Promise<moorhen.Molecule | null> => {
      if (!commandCentre.current) return null;
      const newMolecule = new MoorhenMolecule(requireMoorhenInstance(moorhenInstanceRef));
      newMolecule.setBackgroundColour(backgroundColor);
      newMolecule.defaultBondOptions.smoothness = defaultBondSmoothness;
      try {
        // The element's own dictionaries go on BEFORE Moorhen goes looking
        // for missing monomers. Without this the ligand is bonded from
        // whatever the monomer library happens to return for its code -- or
        // from nothing at all for a novel fragment -- and the scoped
        // read_dictionary_string that the resolver does afterwards arrives
        // too late to be what the molecule was built from. See
        // lib/moorhen-dictionaries: deferring that fetch is the whole point.
        await loadWithDictionaries(newMolecule as any, opts.dictionaries ?? [], () =>
          newMolecule.loadToCootFromString(coordText, molName),
        );
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
    async (
      url: string,
      molName: string,
      opts: { dictionaries?: DictionaryToAttach[] } = {},
    ): Promise<moorhen.Molecule | null> => {
      try {
        const pdbData = await apiText(url);
        return loadSceneStructureFromText(pdbData, molName, url, opts);
      } catch (err) {
        console.warn(`[scene] failed to fetch ${molName} from ${url}:`, err);
        return null;
      }
    },
    [loadSceneStructureFromText],
  );

  const handleFetchSceneFile: SceneFileFetcher = useCallback(
    async (ref: SceneFileRef, fetchOpts) => {
      // The element's scoped dictionaries, which the resolver hands us here.
      // Dropping this argument is what left every campaign summary ligand
      // bonded without its own chemistry.
      const loadOpts = { dictionaries: fetchOpts?.dictionaries ?? [] };
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
          loadOpts,
        );
      }
      // A ccp4i2 fileId is globally unique, so it alone builds the URL;
      // projectId is advisory and must not gate the fetch (the generic
      // wrapper has always treated it that way).
      if (ref.fileId !== undefined) {
        return loadSceneStructure(
          `/api/proxy/ccp4i2/files/${ref.fileId}/download/`,
          ref.name || `file_${ref.fileId}`,
          loadOpts,
        );
      }
      if (ref.pdb) {
        const pdbId = ref.pdb.toLowerCase();
        return loadSceneStructure(
          `/api/proxy/pdbe/entry-files/download/${pdbId}.cif`,
          ref.name || pdbId,
          loadOpts,
        );
      }
      if (ref.url) return loadSceneStructure(ref.url, ref.name || ref.url, loadOpts);
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
        const mapInstance = requireMoorhenInstance(moorhenInstanceRef);
        let newMap: moorhen.Map;
        if (ref.kind === "map") {
          // mode-0 -> float (sane stats); masks also dodge coot's EM cell-clamp.
          const mapBytes = ccp4Mode0ToFloat(bytes as ArrayBuffer);
          newMap = await MoorhenMap.loadToCootFromMapData(
            new Uint8Array(sceneMap.isMask ? ccp4DodgeEmClamp(mapBytes) : mapBytes),
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
    if (COORDINATE_TYPES.has(fileInfo.type)) {
      const url = `/api/proxy/ccp4i2/files/${fileId}/download/`;
      const molName = fileInfo.annotation || fileInfo.job_param_name;
      // A coordinate file brings the dictionaries of the job it belongs to.
      const dictionaries = await fetchDictionaryTexts(await fetchCompanionDictionaryFiles(fileId));
      await fetchMolecule(url, molName, dictionaries);
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

  /**
   * Load a job's best coordinate file, with the job's own dictionaries
   * attached to that molecule alone, and its maps.
   *
   * `asCurrentMember` (the default) is the campaign's own use: the job is the
   * selected member's, so its first dictionary drives the 2D ligand panel and
   * the view centres on it. With it false the job is an addition brought in
   * from the project browser ("load all job outputs"): nothing about the
   * current member is disturbed and the camera stays where it is.
   */
  const fetchJobFiles = async (
    jobId: number,
    opts: { asCurrentMember?: boolean } = {},
  ) => {
    const asCurrentMember = opts.asCurrentMember !== false;
    const files = await apiGet(`files/?job=${jobId}`);
    if (!files || !Array.isArray(files)) return;

    // Filter to only JOB_DIR files (directory=1), exclude imported files (directory=2)
    const jobOutputFiles = files.filter((f: { directory: number }) => f.directory === 1);

    // STEP 1: The job's dictionaries from the database: its own files AND its
    // inputs (a refinement's ligand usually comes from an earlier acedrg job,
    // which a filter over this job's own files never saw). They are attached
    // to this job's molecule only, never to Coot's global store, so a second
    // job brought into the view keeps its own chemistry for a ligand of the
    // same name.
    const ligandDictFiles = await fetchJobDictionaryFiles(jobId);
    const dictionaries = await fetchDictionaryTexts(ligandDictFiles);
    if (asCurrentMember) {
      if (ligandDictFiles.length > 0) {
        // Use the first dictionary file for 2D display in the control panel
        const firstDict = ligandDictFiles[0];
        setLigandDictFileId(firstDict.id);
        const name = firstDict.name?.replace(/\.cif$/i, "") ||
                     firstDict.annotation ||
                     "Ligand";
        setLigandName(name);
        setLigandCodes(candidateLigandCodes(dictionaries.map((d) => d.text)));
      } else {
        setLigandDictFileId(null);
        setLigandName(null);
        setLigandCodes([]);
      }
    }

    // STEP 2: Find and load coordinate files
    // Check for both PDB and mmCIF types
    const coordFiles = jobOutputFiles.filter((f: { type: string }) => COORDINATE_TYPES.has(f.type));
    // Prefer mmCIF (.cif) over PDB (.pdb) for coordinates
    const mmcifFile = coordFiles.find((f: { name: string }) =>
      f.name.toLowerCase().endsWith(".cif")
    );
    const coordFile = mmcifFile || coordFiles[0];
    if (asCurrentMember) setMemberCoordFileId(coordFile ? coordFile.id : null);

    // Load the single best coordinate file
    if (coordFile) {
      const url = `/api/proxy/ccp4i2/files/${coordFile.id}/download/`;
      const molName = coordFile.annotation || coordFile.job_param_name;
      await fetchMolecule(url, molName, dictionaries, { centre: asCurrentMember });
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

  /** Bring an extra job into the view from the project browser. */
  const importJobFiles = (jobId: number) => fetchJobFiles(jobId, { asCurrentMember: false });

  const fetchMolecule = async (
    url: string,
    molName: string,
    dictionaries: DictionaryToAttach[] = [],
    opts: { centre?: boolean } = {},
  ) => {
    if (!commandCentre.current) return;
    const newMolecule = new MoorhenMolecule(requireMoorhenInstance(moorhenInstanceRef));
    newMolecule.setBackgroundColour(backgroundColor);
    newMolecule.defaultBondOptions.smoothness = defaultBondSmoothness;
    try {
      const pdbData = await apiText(url);
      // The molecule's own dictionaries, attached to it alone, before Moorhen
      // goes looking for missing monomers (see lib/moorhen-dictionaries).
      await loadWithDictionaries(newMolecule as any, dictionaries, () =>
        newMolecule.loadToCootFromString(pdbData, molName),
      );
      if (newMolecule.molNo === -1) {
        throw new Error("Cannot read the fetched molecule...");
      }
      newMolecule.uniqueId = url;
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

      if (opts.centre !== false) await newMolecule.centreOn("/*/*/*/*", false, true);
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
      newMap.uniqueId = url;
      // Store the original sub_type for proper labeling and coloring
      (newMap as any).mapSubType = mapSubType;
      // Before addMap: an EM-flagged MTZ map crashes the viewer otherwise.
      primeEmMapHeaderInfo(newMap);
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
    let newMap: moorhen.Map | undefined;
    try {
      // Convert mode-0 (int8) CCP4 maps to float so Moorhen reads sane stats
      // (no-op if already float). For masks, also nudge the P1/orthogonal cell
      // off 90° so coot contours periodically instead of clamping to the cell box.
      let mapData = ccp4Mode0ToFloat(await apiArrayBuffer(url));
      if (opts.isMask) mapData = ccp4DodgeEmClamp(mapData);
      newMap = await MoorhenMap.loadToCootFromMapData(
        new Uint8Array(mapData),
        mapName,
        false,
        requireMoorhenInstance(moorhenInstanceRef),
      );
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

  // View a site: a fresh page on the site's scene, the way the campaign
  // overview opens the summary. In a new tab, because the current session
  // holds a dataset someone is looking at, and a site view replaces the
  // loaded molecules with the site's hits.
  const handleViewSite = useCallback(
    (site: CampaignSite) => {
      window.open(
        `/ccp4i2/moorhen-page/campaign/${campaign.id}?summary=1&site=${site.id}`,
        "_blank"
      );
    },
    [campaign.id]
  );

  // Move to the site named in the URL, once there is a scene to move around.
  //
  // Waits for coot: dispatching an origin before the molecules are drawn puts
  // the camera in the right place and then has it reset underneath us. Fires
  // once, so a user who navigates away from the site is not dragged back by a
  // later re-render.
  const appliedInitialSite = useRef(false);
  useEffect(() => {
    if (appliedInitialSite.current) return;
    if (!initialSiteId || !cootInitialized || sites.length === 0) return;
    const site = sites.find((s) => s.id === initialSiteId);
    if (!site) return;
    appliedInitialSite.current = true;
    handleGoToSite(site);
  }, [initialSiteId, cootInitialized, sites, handleGoToSite]);

  // Save current view as a site
  const handleSaveCurrentAsSite = useCallback(
    async (name: string) => {
      const camera = readCameraState(store.getState() as moorhen.State);
      const newSite: NewCampaignSite = {
        name,
        origin: Array.from(camera.origin).slice(0, 3) as [
          number,
          number,
          number
        ],
        quat: Array.from(camera.quat).slice(0, 4) as [
          number,
          number,
          number,
          number
        ],
        zoom: camera.zoom,
      };
      await onAddSite(newSite);
    },
    [store, onAddSite]
  );

  // Delete a site, by its id. Deleting takes that site's verdicts with it.
  const handleDeleteSite = useCallback(
    async (siteId: number) => {
      await onDeleteSite(siteId);
    },
    [onDeleteSite]
  );

  // Rename a site, and optionally move it to the current view.
  //
  // Addressed by id, not by position: a site's verdicts hang off its id, so a
  // rename has to reach the same row rather than replace a list entry.
  const handleUpdateSite = useCallback(
    async (siteId: number, name: string, updatePosition: boolean) => {
      const changes: Partial<NewCampaignSite> = { name };

      if (updatePosition) {
        const camera = readCameraState(store.getState() as moorhen.State);
        changes.origin = Array.from(camera.origin).slice(0, 3) as [
          number,
          number,
          number
        ];
        changes.quat = Array.from(camera.quat).slice(0, 4) as [
          number,
          number,
          number,
          number
        ];
        changes.zoom = camera.zoom;
      }

      await onUpdateSite(siteId, changes);
    },
    [store, onUpdateSite]
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

  // Record what was found at a site in the selected dataset.
  //
  // This replaces tagging the project with the site's name. A tag could not
  // say that somebody looked and found nothing -- an untagged project was
  // both "empty" and "not yet looked at" -- and it broke on a rename.
  const handleSetVerdict = useCallback(
    async (siteId: number, verdict: SiteVerdict) => {
      if (!selectedMemberProjectId) {
        setMessage("Please select a member project first");
        return;
      }
      const site = sites.find((s) => s.id === siteId);
      try {
        await onSetVerdict(siteId, verdict);
        setMessage(`Recorded ${verdict} at "${site?.name ?? "site"}"`);
      } catch (err) {
        console.error("Failed to record verdict:", err);
        setMessage("Failed to record verdict");
      }
    },
    [selectedMemberProjectId, sites, onSetVerdict, setMessage]
  );

  // Withdraw a verdict: back to nobody having looked, which is not the same
  // as recording "empty".
  const handleClearVerdict = useCallback(
    async (siteId: number) => {
      if (!selectedMemberProjectId) return;
      const site = sites.find((s) => s.id === siteId);
      try {
        await onClearVerdict(siteId);
        setMessage(`Withdrew the verdict at "${site?.name ?? "site"}"`);
      } catch (err) {
        console.error("Failed to withdraw verdict:", err);
        setMessage("Failed to withdraw verdict");
      }
    },
    [selectedMemberProjectId, sites, onClearVerdict, setMessage]
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
        // Step 1: Find observation reflections and the free-R set from this
        // project's top-level jobs.
        //
        // Job outputs (directory 1) are preferred over imports (directory 2),
        // but imports are NOT excluded. A member processed from unmerged data
        // always has an F_SIGF_OUT from aimless; one processed from merged
        // data has none unless the refinement step reindexed, because the
        // observations pass through unchanged -- its only observed file is
        // the imported F_SIGF_IN. Asking for directory=1 alone reported "no
        // reflection data" for every such member.
        type MtzFile = { id: number; job: number; uuid: string; name: string;
                         directory: number; content: number | null;
                         sub_type: number | null };
        const topLevelFiles = async (mimeType: string) =>
          ((await apiGet(
            `files/?type=${mimeType}` +
            `&job__project=${projectId}&job__parent__isnull=true`
          )) as MtzFile[] | null) || [];
        // Outputs before imports, then most recent job first
        const byPreference = (a: MtzFile, b: MtzFile) =>
          a.directory - b.directory || b.job - a.job;

        const obsFiles = await topLevelFiles("application/CCP4-mtz-observed");
        if (obsFiles.length === 0) {
          setMessage("No reflection data found in project");
          return;
        }

        // Among outputs if there are any (else among imports), prefer
        // anomalous data (IPAIR content & 1, FPAIR content & 2) over IMEAN/FMEAN
        const sorted = [...obsFiles].sort(byPreference);
        const candidates = sorted.filter((f) => f.directory === sorted[0].directory);
        const hasAnomalous = (f: MtzFile) => f.content !== null && (f.content & 3) !== 0;
        const reflectionFile = candidates.find(hasAnomalous) || candidates[0];

        // The free-R set. Preferring outputs matters more here than for the
        // observations: the imported set is the campaign's shared one, from
        // the reference crystal, while FREERFLAG_OUT is that set reconciled
        // with this dataset's cell and resolution.
        const freeRFile = (
          await topLevelFiles("application/CCP4-mtz-freerflag")
        ).sort(byPreference)[0];

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

        // Without this the refinement ran with no free set at all, so its
        // R-free meant nothing and could not be compared between siblings.
        if (freeRFile) {
          await apiPost(`jobs/${newJobId}/set_parameter/`, {
            object_path: "servalcat_pipe.inputData.FREERFLAG",
            value: {
              project: projectDbId,
              dbFileId: freeRFile.uuid.replace(/-/g, ""),
            },
          });
        }

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

  // Place the member's ligand at the view centre. Nothing is saved: the
  // result lives in the browser until it is pushed, and it is not pushed
  // automatically because a placed, unfitted ligand in an arbitrary
  // orientation is not something anyone wants silently written into their
  // project.
  const handleAddLigand = useCallback(
    async (code: string) => {
      try {
        await placeLigand(molecules, memberCoordFileId, code);
        dispatch(setRequestDrawScene(true));
        setMessage(`Ligand ${code} added. Push to CCP4i2 to keep it.`, "success");
      } catch (err) {
        setMessage(
          `Could not add ${code}: ${err instanceof Error ? err.message : String(err)}`,
          "error",
        );
      }
    },
    [molecules, memberCoordFileId, dispatch, setMessage]
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
            onViewSite: handleViewSite,
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
            ligandCodes,
            memberCoordFileId,
            onAddLigand: handleAddLigand,
            maps,
            onMapContourLevelChange: handleMapContourLevelChange,
            evaluations,
            onSetVerdict: handleSetVerdict,
            onClearVerdict: handleClearVerdict,
            onFileSelect: fetchFile,
            onJobLoad: importJobFiles,
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
    moorhenInstanceRef,
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
