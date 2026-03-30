/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
import { apiGet, apiText, apiArrayBuffer, apiPost, apiUpload } from "../../api-fetch";
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
  jobId?: number | null;
}

const MoorhenWrapper: React.FC<MoorhenWrapperProps> = ({ fileIds, viewParam, jobId }) => {
  const capabilities = useMoorhenCapabilities();
  const [isSafari] = useState(() => isSafariBrowser());
  const [safariOverride, setSafariOverride] = useState(false);
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
        console.warn("[fetchMolecule] Ligands representation failed");
      }

      await newMolecule.centreOn("/*/*/*/*", false, true);
      dispatch(addMolecule(newMolecule));
      dispatch(showMolecule({ molNo: newMolecule.molNo } as any));
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
              changesMolecules: [newMolecule.molNo],
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
        formData.append("objectPath", "servalcat_pipe.inputData.XYZIN");
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

  // Custom side panel containing our control panel
  const extraSidePanels: Record<string, MoorhenPanel> = useMemo(() => ({
    ccp4i2Controls: {
      icon: "MatSymSettings",
      label: "CCP4i2",
      panelContent: (
        <MoorhenControlPanel
          onFileSelect={fetchFile}
          onJobLoad={fetchJobFiles}
          getViewUrl={getViewUrl}
          molecules={molecules}
          maps={maps}
          onMapContourLevelChange={handleMapContourLevelChange}
          onRunServalcat={jobId ? handleRunServalcat : undefined}
          servalcatStatus={servalcatStatus}
        />
      ),
    },
  }), [fetchFile, fetchJobFiles, getViewUrl, molecules, maps, handleMapContourLevelChange, jobId, handleRunServalcat, servalcatStatus]);

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
    // Browser capabilities limited - attempting to load anyway
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
