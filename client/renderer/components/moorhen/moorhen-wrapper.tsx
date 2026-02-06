"use client";

import {
  addMolecule,
  addMap,
  setActiveMap,
  setWidth,
  setHeight,
  setTheme,
  setBackgroundColor,
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
import { MoorhenControlPanel } from "./moorhen-control-panel";
import { apiGet, apiText, apiArrayBuffer } from "../../api-fetch";
import { useTheme } from "../../theme/theme-provider";
import { useMoorhenViewState } from "../../hooks/use-moorhen-view-state";
import {
  MoorhenFallback,
  MoorhenErrorBoundary,
  useMoorhenCapabilities,
  isSafariBrowser,
  SafariExperimentalWarning,
} from "./moorhen-capability-check";

export interface MoorhenWrapperProps {
  fileIds?: number[];
  viewParam?: string | null;
}

const MoorhenWrapper: React.FC<MoorhenWrapperProps> = ({ fileIds, viewParam }) => {
  const capabilities = useMoorhenCapabilities();
  const [isSafari] = useState(() => isSafariBrowser());
  const [safariOverride, setSafariOverride] = useState(false);
  const dispatch = useDispatch();
  const theme = useTheme();

  // View state hook for URL parameter support
  const { getViewUrl } = useMoorhenViewState({
    viewParam: viewParam ?? null,
    onViewRestored: () => console.log("View state restored from URL"),
  });
  const glRef: RefObject<webGL.MGWebGL | null> = useRef(null);
  const commandCentre = useRef<null | moorhen.CommandCentre>(null);
  const moleculesRef = useRef<null | moorhen.Molecule[]>(null);
  const mapsRef = useRef<null | moorhen.Map[]>(null);
  const activeMapRef = useRef<null | moorhen.Map>(null);
  const lastHoveredAtom = useRef<null | moorhen.HoveredAtom>(null);
  const prevActiveMoleculeRef = useRef<null | moorhen.Molecule>(null);
  const timeCapsuleRef = useRef(null);
  const cootInitialized = useSelector(
    (state: moorhen.State) => state.generalStates.cootInitialized
  );
  const store = useStore();
  const { cootModule, cootModuleError } = useCCP4i2Window();

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
  }, [theme.mode]);

  const monomerLibraryPath =
    "https://raw.githubusercontent.com/MonomerLibrary/monomers/master/";

  const backgroundColor = useSelector(
    (state: moorhen.State) => state.sceneSettings.backgroundColor
  );
  const defaultBondSmoothness = useSelector(
    (state: moorhen.State) => state.sceneSettings.defaultBondSmoothness
  );

  const [windowWidth, setWindowWidth] = useState<number>(window.innerWidth);
  const [windowHeight, setWindowHeight] = useState<number>(window.innerHeight);

  // Calculate dimensions for the split screen layout
  const rightPanelWidth = 60 * 8; // Approximate character width in pixels (8px per character)
  const leftPanelWidth = useMemo(() => {
    return windowWidth - rightPanelWidth;
  }, [windowWidth, rightPanelWidth]);

  const setMoorhenDimensions = useCallback(() => {
    const result = [leftPanelWidth, windowHeight];
    return result;
  }, [leftPanelWidth, windowHeight]);

  // URL prefix for Moorhen to load its resources
  // In web browsers, use API route for CORP headers (COEP compatibility)
  // In Electron, serve directly from public/baby-gru
  const isElectron = typeof window !== "undefined" && !!(window as any).electronAPI;
  const urlPrefix = isElectron ? "/baby-gru" : "/api/moorhen/baby-gru";

  const collectedProps = useMemo(() => ({
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
  }), [setMoorhenDimensions, urlPrefix]);

  // Note: Don't subscribe to glRef state here - it changes every frame during rotation
  // and would cause constant re-renders. Access origin directly from store when needed.
  const getOrigin = useCallback(() => {
    const state = store.getState() as moorhen.State;
    return state.glRef.origin;
  }, [store]);

  const handleResize = () => {
    setWindowWidth(window.innerWidth);
    setWindowHeight(window.innerHeight - 30);
    console.log("Window resized");
  };

  useEffect(() => {
    //What to do when the component mounts
    console.log("MoorhenWrapper mounted");

    window.addEventListener("resize", handleResize);

    return () => {
      console.log("MoorhenWrapper unmounted");
      window.removeEventListener("resize", handleResize);
    };
  }, []);

  useEffect(() => {
    if (fileIds && cootInitialized && cootModule) {
      fileIds.forEach((fileId) => {
        fetchFile(fileId);
      });
    }
  }, [fileIds, cootInitialized, cootModule]);

  // Track if we've logged initialization to avoid repeat logs
  const hasLoggedInit = useRef(false);

  // One-time initialization message when Coot is ready
  useEffect(() => {
    if (cootInitialized && !hasLoggedInit.current) {
      console.log("Coot is initialized, you can now load molecules and maps.");
      hasLoggedInit.current = true;
    }
  }, [cootInitialized]);

  // Handle dimension updates separately - these can change frequently
  useEffect(() => {
    if (cootInitialized) {
      dispatch(setWidth(leftPanelWidth));
      dispatch(setHeight(windowHeight - 75));
    }
  }, [cootInitialized, leftPanelWidth, windowHeight, dispatch]);

  // Handle theme changes separately
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
      // Anomalous maps (3) behave like difference maps for contouring
      const mapSubType = fileInfo.sub_type || 1;
      await fetchMap(url, molName, mapSubType);
    } else if (fileInfo.type === "application/refmac-dictionary") {
      const url = `/api/proxy/ccp4i2/files/${fileId}/download/`;
      await fetchDict(url);
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
      // Fetch PDB data using authenticated api-fetch, then load into Moorhen
      const pdbData = await apiText(url);
      await newMolecule.loadToCootFromString(pdbData, molName);
      if (newMolecule.molNo === -1) {
        throw new Error("Cannot read the fetched molecule...");
      }
      newMolecule.uniqueId = url; // Use URL as unique identifier
      await newMolecule.addRepresentation("CBs", "/*/*/*/*");
      await newMolecule.addRepresentation("ligands", "/*/*/*/*");
      await newMolecule.centreOn("/*/*/*/*", false, true);

      dispatch(addMolecule(newMolecule));
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch PDB entry from ${url}, doing nothing...`);
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
      // Fetch MTZ data using authenticated api-fetch, then load into Moorhen
      const mtzData = await apiArrayBuffer(url);
      await newMap.loadToCootFromMtzData(new Uint8Array(mtzData), mapName, {
        F: "F",
        PHI: "PHI",
        useWeight: false,
        isDifference: isDiffMap,
      });
      newMap.uniqueId = url; // Use URL as unique identifier
      // Store the original sub_type for proper labeling and coloring
      (newMap as any).mapSubType = mapSubType;
      // Set custom colors for anomalous maps (orange/purple instead of green/red)
      if (mapSubType === 3) {
        newMap.defaultPositiveMapColour = { r: 1.0, g: 0.65, b: 0.0 }; // Orange
        newMap.defaultNegativeMapColour = { r: 0.6, g: 0.3, b: 0.8 }; // Purple
      }
      if (newMap.molNo === -1)
        throw new Error("Cannot read the fetched map...");
      dispatch(addMap(newMap));
      dispatch(setActiveMap(newMap));
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch map from ${url}`);
    }
    return newMap;
  };

  const fetchDict = async (url: string) => {
    if (!commandCentre.current) return;
    const fileContent = await apiText(url);
    await commandCentre.current.cootCommand(
      {
        returnType: "status",
        command: "read_dictionary_string",
        commandArgs: [fileContent, -999999],
        changesMolecules: [],
      },
      false
    );
    const instanceName = "LIG";
    const result = (await commandCentre.current.cootCommand(
      {
        returnType: "status",
        command: "get_monomer_and_position_at",
        commandArgs: [instanceName, -999999, ...getOrigin().map((coord: number) => -coord)],
      },
      true
    )) as moorhen.WorkerResponse<number>;
    if (result.data.result.status === "Completed") {
      const newMolecule = new MoorhenMolecule(
        commandCentre as RefObject<moorhen.CommandCentre>,
        glRef as RefObject<webGL.MGWebGL>,
        store,
        monomerLibraryPath
      );
      newMolecule.uniqueId = url; // Use URL as unique identifier
      newMolecule.molNo = result.data.result.result;
      newMolecule.name = instanceName;
      newMolecule.setBackgroundColour(backgroundColor);
      newMolecule.defaultBondOptions.smoothness = defaultBondSmoothness;
      newMolecule.coordsFormat = "mmcif";
      await Promise.all([
        newMolecule.fetchDefaultColourRules(),
        newMolecule.addDict(fileContent),
      ]);
      newMolecule.centreAndAlignViewOn("/*/*/*/*", false, 100);
      await newMolecule.fetchIfDirtyAndDraw("ligands");
      dispatch(addMolecule(newMolecule));
    }
  };

  // Show fallback if Coot module failed to load
  if (cootModuleError) {
    return (
      <MoorhenFallback
        reason="load_error"
        error={cootModuleError}
        capabilities={capabilities}
      />
    );
  }

  // Show Safari warning - capabilities may look OK but WASM threading crashes
  const isElectronEnv = typeof window !== "undefined" && !!(window as any).electronAPI;
  if (isSafari && !safariOverride && !isElectronEnv) {
    return <SafariExperimentalWarning onProceed={() => setSafariOverride(true)} />;
  }

  // Show fallback if browser capabilities are missing and we're not in Electron
  // But only after we've checked capabilities (undefined during SSR)
  if (capabilities && !capabilities.isSupported && !isElectronEnv) {
    // Still attempt to render - the MoorhenErrorBoundary will catch failures
    // This allows browsers that have improved their support to work
    console.log("[Moorhen] Browser capabilities limited, attempting to load anyway...");
  }

  return (
    <MoorhenErrorBoundary
      fallback={
        <MoorhenFallback
          reason="runtime_error"
          capabilities={capabilities}
        />
      }
    >
      {store && cootModule && (
        <div
          style={{
            display: "flex",
            width: "100%",
            height: "100%",
            flexDirection: "row",
          }}
        >
          {/* Left panel - MoorhenContainer taking remaining space */}
          <div
            style={{
              flex: 1,
              minWidth: 0, // Allows flex item to shrink below content size
              height: "calc(100% - 120px)",
            }}
          >
            <MoorhenContainer {...collectedProps} />
          </div>

          {/* Right panel - Fixed width 80 characters */}
          <div
            style={{
              width: `${rightPanelWidth}px`,
              minWidth: `${rightPanelWidth}px`,
              maxWidth: `${rightPanelWidth}px`,
              minHeight: "calc(100% - 150px)",
              borderLeft: "1px solid #ddd",
              padding: "0px",
              fontSize: "14px",
              fontFamily: "monospace",
              overflowY: "auto",
              overflowX: "hidden",
              boxSizing: "border-box",
            }}
          >
            <MoorhenControlPanel onFileSelect={fetchFile} getViewUrl={getViewUrl} />
          </div>
        </div>
      )}
    </MoorhenErrorBoundary>
  );
};

export default MoorhenWrapper;
