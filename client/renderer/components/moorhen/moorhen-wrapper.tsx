"use client";

import {
  addMolecule,
  addMap,
  setActiveMap,
  setContourLevel,
  setTheme,
  setBackgroundColor,
  setRequestDrawScene,
  MoorhenContainer,
  MoorhenMolecule,
  MoorhenMap,
} from "moorhen";
// @ts-ignore - moorhen 0.23 exports may lack .d.ts depending on build
import { MoorhenInstanceProvider, setShownSidePanel } from "moorhen";
// @ts-ignore - moorhen 0.23 type may lack .d.ts depending on build
import type { MoorhenPanel } from "moorhen";

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

  // Auto-open the CCP4i2 controls side panel
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

  const fetchMolecule = useCallback(async (url: string, molName: string) => {
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
        console.log("[fetchMolecule] Ligands representation failed");
      }

      await newMolecule.centreOn("/*/*/*/*", false, true);
      dispatch(addMolecule(newMolecule));
    } catch (err) {
      console.warn(err);
      console.warn(`Cannot fetch PDB entry from ${url}, doing nothing...`);
    }
  }, [commandCentre, store, monomerLibraryPath, backgroundColor, defaultBondSmoothness, dispatch]);

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
          newMolecule.centreAndAlignViewOn("/*/*/*/*", false, 100);
          centredFirst = true;
        }
        await newMolecule.fetchIfDirtyAndDraw("ligands");
        dispatch(addMolecule(newMolecule));
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

  // Custom side panel containing our control panel
  const extraSidePanels: Record<string, MoorhenPanel> = useMemo(() => ({
    ccp4i2Controls: {
      icon: "MatSymSettings",
      label: "CCP4i2",
      panelContent: (
        <MoorhenControlPanel
          onFileSelect={fetchFile}
          getViewUrl={getViewUrl}
          molecules={molecules}
          maps={maps}
          onMapContourLevelChange={handleMapContourLevelChange}
        />
      ),
    },
  }), [fetchFile, getViewUrl, molecules, maps, handleMapContourLevelChange]);

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

  // Show Safari warning - capabilities may look OK but WASM threading crashes
  const isElectronEnv = typeof window !== "undefined" && !!(window as any).electronAPI;
  if (isSafari && !safariOverride && !isElectronEnv) {
    return <SafariExperimentalWarning onProceed={() => setSafariOverride(true)} />;
  }

  if (capabilities && !capabilities.isSupported && !isElectronEnv) {
    console.log("[Moorhen] Browser capabilities limited, attempting to load anyway...");
  }

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
          <MoorhenInstanceProvider>
            <MoorhenContainer {...collectedProps} />
          </MoorhenInstanceProvider>
        )}
      </MoorhenErrorBoundary>
    </div>
  );
};

export default MoorhenWrapper;
