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
import { apiGet, apiText } from "../../api-fetch";
import { useTheme } from "../../theme/theme-provider";

export interface MoorhenWrapperProps {
  fileIds?: number[];
}

const MoorhenWrapper: React.FC<MoorhenWrapperProps> = ({ fileIds }) => {
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

  useEffect(() => {
    if (cootInitialized) {
      console.log("Coot is initialized, you can now load molecules and maps.");
      dispatch(setWidth(leftPanelWidth));
      dispatch(setHeight(windowHeight - 75));
      dispatch(
        setBackgroundColor(theme.mode === "light" ? [1, 1, 1, 1] : [0, 0, 0, 1])
      );
      handleResize();
    }
  }, [cootInitialized, leftPanelWidth, windowHeight, theme.mode]);

  const fetchFile = async (fileId: number) => {
    const fileInfo = await apiGet(`files/${fileId}`);
    console.log(fileInfo);
    if (!fileInfo) {
      console.warn(`File with ID ${fileId} not found.`);
      return;
    }
    if (fileInfo.type === "chemical/x-pdb") {
      const url = `/api/proxy/files/${fileId}/download/`;
      const molName = fileInfo.annotation || fileInfo.job_param_name;
      await fetchMolecule(url, molName);
    } else if (fileInfo.type === "application/CCP4-mtz-map") {
      const url = `/api/proxy/files/${fileId}/download/`;
      const molName = fileInfo.name || fileInfo.job_param_name;
      const isDiffMap = fileInfo.sub_type !== 1 || false;
      await fetchMap(url, molName, isDiffMap);
    } else if (fileInfo.type === "application/refmac-dictionary") {
      const url = `/api/proxy/files/${fileId}/download/`;
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
      await newMolecule.loadToCootFromURL(url, molName);
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
      await newMap.loadToCootFromMtzURL(url, mapName, {
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
