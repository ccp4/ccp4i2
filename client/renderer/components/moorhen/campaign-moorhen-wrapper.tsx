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
  addMap,
  removeMolecule,
  removeMap,
  setActiveMap,
  setContourLevel,
  setWidth,
  setHeight,
  setTheme,
  setBackgroundColor,
  setOrigin,
  setQuat,
  setZoom,
  setRequestDrawScene,
  MoorhenContainer,
  MoorhenMolecule,
  MoorhenMap,
} from "moorhen";

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
import { useCCP4i2Window } from "../../app-context";
import { CampaignControlPanel } from "./campaign-control-panel";
import { apiText, apiArrayBuffer, apiGet } from "../../api-fetch";
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

// Detect Safari browser (has WASM threading issues with Moorhen)
const isSafariBrowser = (): boolean => {
  if (typeof navigator === "undefined") return false;
  const ua = navigator.userAgent;
  const isSafari = /^((?!chrome|android).)*safari/i.test(ua);
  const isIOS = /iPad|iPhone|iPod/.test(ua);
  return isSafari || isIOS;
};

// Safari warning component
const SafariWarning: React.FC = () => (
  <div
    style={{
      display: "flex",
      flexDirection: "column",
      alignItems: "center",
      justifyContent: "center",
      height: "100%",
      padding: "40px",
      textAlign: "center",
      backgroundColor: "#fff3cd",
      border: "1px solid #ffc107",
      borderRadius: "8px",
      margin: "20px",
    }}
  >
    <h2 style={{ color: "#856404", marginBottom: "16px" }}>
      Browser Not Supported
    </h2>
    <p style={{ color: "#856404", maxWidth: "600px", lineHeight: "1.6" }}>
      The Moorhen molecular viewer requires WebAssembly threading features that
      are not fully supported in Safari. Please use{" "}
      <strong>Google Chrome</strong>, <strong>Microsoft Edge</strong>, or{" "}
      <strong>Firefox</strong> to view molecular structures.
    </p>
  </div>
);

type FileSource =
  | { type: "none" }
  | { type: "files"; fileIds: number[] }
  | { type: "job"; jobId: number };

export interface CampaignMoorhenWrapperProps {
  campaign: ProjectGroup;
  fileSource: FileSource;
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
  viewParam,
  sites,
  onUpdateSites,
  memberProjects,
  selectedMemberProjectId,
  onSelectMemberProject,
  parentProject,
}) => {
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
    onViewRestored: () => console.log("View state restored from URL"),
    representations: visibleRepresentations,
  });

  const hasInitializedReps = useRef(false);

  const glRef: RefObject<webGL.MGWebGL | null> = useRef(null);
  const commandCentre = useRef<null | moorhen.CommandCentre>(null);
  const moleculesRef = useRef<null | moorhen.Molecule[]>(null);
  const mapsRef = useRef<null | moorhen.Map[]>(null);
  const activeMapRef = useRef<null | moorhen.Map>(null);
  const lastHoveredAtom = useRef<null | moorhen.HoveredAtom>(null);
  const prevActiveMoleculeRef = useRef<null | moorhen.Molecule>(null);
  const timeCapsuleRef = useRef(null);
  const loadedFileSource = useRef<FileSource | null>(null);

  // Ligand dictionary file ID for 2D structure display
  const [ligandDictFileId, setLigandDictFileId] = useState<number | null>(null);
  const [ligandName, setLigandName] = useState<string | null>(null);
  // Store loaded dictionary content so we can add it to molecules
  const loadedDictContent = useRef<string | null>(null);

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
  const { cootModule } = useCCP4i2Window();

  // Initialize representations from URL if available (after molecules load)
  useEffect(() => {
    if (!hasInitializedReps.current && initialRepresentations && molecules.length > 0) {
      setVisibleRepresentations(initialRepresentations);
      applyRepresentationsToMolecules(initialRepresentations, molecules);
      hasInitializedReps.current = true;
    }
  }, [initialRepresentations, molecules, applyRepresentationsToMolecules]);

  useEffect(() => {
    if (typeof window !== "undefined") {
      (window as any).CCP4Module = cootModule;
    }
  }, [cootModule]);

  useEffect(() => {
    dispatch(
      setBackgroundColor(theme.mode === "light" ? [1, 1, 1, 1] : [0, 0, 0, 1])
    );
    dispatch(setTheme(theme.mode === "light" ? "flatly" : "darkly"));
  }, [theme.mode, dispatch]);

  const monomerLibraryPath =
    "https://raw.githubusercontent.com/MonomerLibrary/monomers/master/";

  const backgroundColor = useSelector(
    (state: moorhen.State) => state.sceneSettings.backgroundColor
  );
  const defaultBondSmoothness = useSelector(
    (state: moorhen.State) => state.sceneSettings.defaultBondSmoothness
  );

  const [windowWidth, setWindowWidth] = useState<number>(
    typeof window !== "undefined" ? window.innerWidth : 1200
  );
  const [windowHeight, setWindowHeight] = useState<number>(
    typeof window !== "undefined" ? window.innerHeight : 800
  );

  const rightPanelWidth = 60 * 8;
  const leftPanelWidth = useMemo(() => {
    return windowWidth - rightPanelWidth;
  }, [windowWidth, rightPanelWidth]);

  const setMoorhenDimensions = useCallback(() => {
    const result = [leftPanelWidth, windowHeight];
    return result;
  }, [leftPanelWidth, windowHeight]);

  const isElectron =
    typeof window !== "undefined" && !!(window as any).electronAPI;
  const urlPrefix = isElectron ? "/baby-gru" : "/api/moorhen/baby-gru";

  const collectedProps = useMemo(
    () => ({
      glRef,
      timeCapsuleRef,
      commandCentre,
      moleculesRef,
      mapsRef,
      activeMapRef,
      lastHoveredAtom,
      prevActiveMoleculeRef,
      setMoorhenDimensions,
      monomerLibraryPath,
      urlPrefix,
    }),
    [setMoorhenDimensions, urlPrefix]
  );

  const getOrigin = useCallback(() => {
    const state = store.getState() as moorhen.State;
    return (state as unknown as { glRef: { origin: number[] } }).glRef.origin;
  }, [store]);

  const handleResize = () => {
    setWindowWidth(window.innerWidth);
    setWindowHeight(window.innerHeight - 30);
  };

  useEffect(() => {
    window.addEventListener("resize", handleResize);
    return () => {
      window.removeEventListener("resize", handleResize);
    };
  }, []);

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
    loadedDictContent.current = null;
    // Reset representation state to default
    setVisibleRepresentations(["CRs"]);
    hasInitializedReps.current = false;
  }, [molecules, maps, dispatch]);

  // Load files when fileSource changes
  useEffect(() => {
    if (!cootInitialized || !cootModule) return;

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
  }, [fileSource, cootInitialized, cootModule, cleanupLoadedContent]);

  useEffect(() => {
    if (cootInitialized) {
      dispatch(setWidth(leftPanelWidth));
      dispatch(setHeight(windowHeight - 75));
    }
  }, [cootInitialized, leftPanelWidth, windowHeight, dispatch]);

  useEffect(() => {
    if (cootInitialized) {
      dispatch(
        setBackgroundColor(theme.mode === "light" ? [1, 1, 1, 1] : [0, 0, 0, 1])
      );
    }
  }, [cootInitialized, theme.mode, dispatch]);

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
    }
  };

  const fetchJobFiles = async (jobId: number) => {
    const files = await apiGet(`files/?job=${jobId}`);
    if (!files || !Array.isArray(files)) return;

    console.log("[fetchJobFiles] All files with types:", files.map((f: any) =>
      `${f.name} (type=${f.type}, dir=${f.directory})`
    ));

    // Filter to only JOB_DIR files (directory=1), exclude imported files (directory=2)
    const jobOutputFiles = files.filter((f: { directory: number }) => f.directory === 1);
    console.log("[fetchJobFiles] JOB_DIR files:", jobOutputFiles.map((f: any) =>
      `${f.name} (type=${f.type})`
    ));

    // STEP 1: Load ligand dictionary FIRST (before coordinates)
    // This ensures coot understands ligand geometry when parsing coordinates
    // Check all files (including imported) for dictionary
    const ligandDictFile = files.find(
      (f: { type: string }) => f.type === "application/refmac-dictionary"
    );
    if (ligandDictFile) {
      console.log("[fetchJobFiles] Loading ligand dictionary FIRST:", ligandDictFile.name);
      const dictUrl = `/api/proxy/ccp4i2/files/${ligandDictFile.id}/download/`;
      await fetchDict(dictUrl);
      // Store file ID for 2D display
      setLigandDictFileId(ligandDictFile.id);
      // Extract ligand name from filename (e.g., "LIG.cif" -> "LIG")
      const name = ligandDictFile.name?.replace(/\.cif$/i, "") ||
                   ligandDictFile.annotation ||
                   "Ligand";
      setLigandName(name);
    } else {
      loadedDictContent.current = null;
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
    console.log("[fetchJobFiles] Coord files (PDB/CIF):", coordFiles.map((f: any) =>
      `${f.name} (type=${f.type})`
    ));

    // Prefer mmCIF (.cif) over PDB (.pdb) for coordinates
    const mmcifFile = coordFiles.find((f: { name: string }) =>
      f.name.toLowerCase().endsWith(".cif")
    );
    const coordFile = mmcifFile || coordFiles[0];
    console.log("[fetchJobFiles] Selected coord file:", coordFile?.name);

    // Load the single best coordinate file
    if (coordFile) {
      const url = `/api/proxy/ccp4i2/files/${coordFile.id}/download/`;
      const molName = coordFile.annotation || coordFile.job_param_name;
      await fetchMolecule(url, molName);
    }

    // STEP 3: Load map files
    for (const file of jobOutputFiles) {
      if (file.type === "application/CCP4-mtz-map") {
        const url = `/api/proxy/ccp4i2/files/${file.id}/download/`;
        const molName = file.name || file.job_param_name;
        // subType: 1=normal, 2=difference, 3=anomalous difference
        const mapSubType = file.sub_type || 1;
        await fetchMap(url, molName, mapSubType);
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
      console.log("[fetchDict] Loaded dictionary globally");
      // Store content so we can add it to molecules later
      loadedDictContent.current = fileContent;
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
      glRef as RefObject<webGL.MGWebGL>,
      store,
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

      // Add dictionary to molecule if we have one loaded
      // This ensures the molecule understands ligand geometry
      if (loadedDictContent.current) {
        try {
          await newMolecule.addDict(loadedDictContent.current);
          console.log("[fetchMolecule] Added dictionary to molecule");
        } catch (err) {
          console.warn("[fetchMolecule] Failed to add dictionary:", err);
        }
      }

      // Try ribbon representation first (better for protein overview)
      // Fall back to CBs if ribbons fail (e.g., no protein backbone)
      try {
        await newMolecule.addRepresentation("CRs", "/*/*/*/*");
      } catch {
        console.log("[fetchMolecule] Ribbons failed, falling back to CBs");
        await newMolecule.addRepresentation("CBs", "/*/*/*/*");
      }

      // Always try to add ligand representation
      try {
        await newMolecule.addRepresentation("ligands", "/*/*/*/*");
      } catch {
        console.log("[fetchMolecule] Ligands representation failed");
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
      glRef as RefObject<webGL.MGWebGL>,
      store
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
      dispatch(setActiveMap(newMap));
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch map from ${url}`);
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

  if (isSafari) {
    return <SafariWarning />;
  }

  return (
    store &&
    cootModule && (
      <div
        style={{
          display: "flex",
          width: "100%",
          height: "100%",
          flexDirection: "row",
        }}
      >
        <div
          style={{
            flex: 1,
            minWidth: 0,
            height: "calc(100% - 120px)",
          }}
        >
          <MoorhenContainer {...collectedProps} />
        </div>

        <div
          style={{
            width: `${rightPanelWidth}px`,
            minWidth: `${rightPanelWidth}px`,
            maxWidth: `${rightPanelWidth}px`,
            height: "calc(100vh - 40px)",
            borderLeft: "1px solid #ddd",
            padding: "0px",
            fontSize: "14px",
            fontFamily: "monospace",
            overflow: "hidden",
            boxSizing: "border-box",
          }}
        >
          <CampaignControlPanel
            campaign={campaign}
            sites={sites}
            onGoToSite={handleGoToSite}
            onSaveCurrentAsSite={handleSaveCurrentAsSite}
            onUpdateSite={handleUpdateSite}
            onDeleteSite={handleDeleteSite}
            memberProjects={memberProjects}
            selectedMemberProjectId={selectedMemberProjectId}
            onSelectMemberProject={onSelectMemberProject}
            parentProject={parentProject}
            getViewUrl={getViewUrl}
            molecules={molecules}
            visibleRepresentations={visibleRepresentations}
            onRepresentationsChange={setVisibleRepresentations}
            ligandDictFileId={ligandDictFileId}
            ligandName={ligandName}
            maps={maps}
            onMapContourLevelChange={handleMapContourLevelChange}
            onTagProjectWithSite={handleTagProjectWithSite}
            onFileSelect={fetchFile}
          />
        </div>
      </div>
    )
  );
};

export default CampaignMoorhenWrapper;
