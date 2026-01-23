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
      are not fully supported in Safari. Please use <strong>Google Chrome</strong>,{" "}
      <strong>Microsoft Edge</strong>, or <strong>Firefox</strong> to view
      molecular structures.
    </p>
    <p style={{ color: "#856404", marginTop: "16px", fontSize: "14px" }}>
      Alternatively, use the CCP4i2 desktop application (Electron) which works
      on all platforms.
    </p>
  </div>
);

export interface MoorhenWrapperProps {
  fileIds?: number[];
}

const MoorhenWrapper: React.FC<MoorhenWrapperProps> = ({ fileIds }) => {
  const [isSafari] = useState(() => isSafariBrowser());
  const dispatch = useDispatch();
  const theme = useTheme();
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
  const { cootModule } = useCCP4i2Window();

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
  const baseUrl = "https://www.ebi.ac.uk/pdbe/entry-files";

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

  const collectedProps = {
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
  };

  const { origin } = useSelector((state: moorhen.State) => state.glRef);

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
      const isDiffMap = fileInfo.sub_type === 2;
      await fetchMap(url, molName, isDiffMap);
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
    isDiffMap: boolean = false
  ) => {
    if (!commandCentre.current) return;
    const newMap = new MoorhenMap(
      commandCentre as RefObject<moorhen.CommandCentre>,
      glRef as RefObject<webGL.MGWebGL>,
      store
    );
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

  const fetchDict = async (
    url: string,
    newMolecules: moorhen.Molecule[] = []
  ) => {
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
        commandArgs: [instanceName, -999999, ...origin.map((coord) => -coord)],
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

  // Show Safari warning instead of Moorhen viewer
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
          <MoorhenControlPanel onFileSelect={fetchFile} />
        </div>
      </div>
    )
  );
};

export default MoorhenWrapper;
