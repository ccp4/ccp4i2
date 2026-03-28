/*
 * Copyright (C) 2026 Newcastle University
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
import React, { useEffect, useMemo } from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { useJob } from "../../../../utils";
import { useBoolToggle } from "../../task-elements/shared-hooks";
import { InputDataTab } from "./InputDataTab";
import { ImportantOptionsTab } from "./ImportantOptionsTab";
import { AdditionalOptionsTab } from "./AdditionalOptionsTab";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // ─── Auto-sync: REFERENCE_FOR_AIMLESS tracks MODE ──────────────────
  const { syncTo: syncToAimlessRef } = useTaskItem("REFERENCE_FOR_AIMLESS");

  // ─── Enum values for conditional rendering ─────────────────────────
  const { value: mode } = useTaskItem("MODE");
  const { value: chooseMode } = useTaskItem("CHOOSE_MODE");
  const { value: referenceDataset } = useTaskItem("REFERENCE_DATASET");
  const { value: sdcorrOptions } = useTaskItem("SDCORRECTION_OPTIONS");
  const { value: scalingProtocol } = useTaskItem("SCALING_PROTOCOL");
  const { value: scalesRotType } = useTaskItem("SCALES_ROTATION_TYPE");
  const { value: scalesBrotType } = useTaskItem("SCALES_BROTATION_TYPE");
  const { value: scalesTileType } = useTaskItem("SCALES_TILETYPE");
  const { value: runMode } = useTaskItem("RUN_MODE");

  // ─── Boolean toggles ───────────────────────────────────────────────
  const analysisOverride = useBoolToggle(useTaskItem, "ANALYSIS_OVERRIDE");
  const intensitiesOverride = useBoolToggle(
    useTaskItem,
    "INTENSITIES_OVERRIDE"
  );
  const sdcorrOverride = useBoolToggle(useTaskItem, "SDCORRECTION_OVERRIDE");
  const sdcorrRefine = useBoolToggle(useTaskItem, "SDCORRECTION_REFINE");
  const scalingDetails = useBoolToggle(useTaskItem, "SCALING_DETAILS");
  const outputUnmerged = useBoolToggle(useTaskItem, "OUTPUT_UNMERGED");
  const removeLattice = useBoolToggle(
    useTaskItem,
    "REMOVE_LATTICE_CENTERING"
  );
  const keepLattice = useBoolToggle(useTaskItem, "KEEP_LATTICE_CENTERING");
  const outlierOverride = useBoolToggle(useTaskItem, "OUTLIER_OVERRIDE");
  const expertOptions = useBoolToggle(useTaskItem, "EXPERT_OPTIONS");

  // ─── Auto-sync REFERENCE_FOR_AIMLESS with MODE ────────────────────
  useEffect(() => {
    syncToAimlessRef(mode === "MATCH");
  }, [mode, syncToAimlessRef]);

  // ─── Visibility helpers ────────────────────────────────────────────
  const vis = useMemo(
    () => ({
      // Tab 1: Symmetry modes
      isChooseMode: () => mode === "CHOOSE",
      isChooseSolution: () =>
        mode === "CHOOSE" && chooseMode === "SOLUTION_NO",
      isChooseSpacegroup: () =>
        mode === "CHOOSE" &&
        ["SPACEGROUP", "REINDEX_SPACE"].includes(chooseMode),
      isReindexSpace: () =>
        mode === "CHOOSE" && chooseMode === "REINDEX_SPACE",
      isChooseLauegroup: () =>
        mode === "CHOOSE" && chooseMode === "LAUEGROUP",
      // Tab 2: SD correction
      showSimilarity: () => sdcorrOptions === "SIMILAR",
      // Tab 3: Scaling protocol
      hasRotationScales: () =>
        ["ROTATION", "SECONDARY"].includes(scalingProtocol),
      hasSecondaryBeam: () => scalingProtocol === "SECONDARY",
      showTileXY: () =>
        scalingProtocol === "SECONDARY" &&
        !["DEFAULT", "NONE"].includes(scalesTileType),
      rotSpacing: () => scalesRotType === "SPACING",
      rotNbins: () => scalesRotType === "NBINS",
      brotSpacing: () => scalesBrotType === "SPACING",
      brotNbins: () => scalesBrotType === "NBINS",
      // Tab 3: Lattice centering
      showLatticeThreshold: () =>
        !removeLattice.value && !keepLattice.value,
      // Tab 3: Run mode
      showBatchList: () => runMode === "BYRANGE",
    }),
    [
      mode,
      chooseMode,
      sdcorrOptions,
      scalingProtocol,
      scalesRotType,
      scalesBrotType,
      scalesTileType,
      removeLattice.value,
      keepLattice.value,
      runMode,
    ]
  );

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        <CCP4i2Tab label="Input Data">
          <InputDataTab
            {...props}
            mode={mode}
            chooseMode={chooseMode}
            referenceDataset={referenceDataset}
            vis={vis}
          />
        </CCP4i2Tab>
        <CCP4i2Tab label="Important Options">
          <ImportantOptionsTab
            {...props}
            sdcorrOptions={sdcorrOptions}
            analysisOverride={analysisOverride}
            intensitiesOverride={intensitiesOverride}
            sdcorrOverride={sdcorrOverride}
            sdcorrRefine={sdcorrRefine}
            scalingDetails={scalingDetails}
            outputUnmerged={outputUnmerged}
            vis={vis}
          />
        </CCP4i2Tab>
        <CCP4i2Tab label="Additional Options">
          <AdditionalOptionsTab
            {...props}
            scalingProtocol={scalingProtocol}
            scalesRotType={scalesRotType}
            scalesBrotType={scalesBrotType}
            removeLattice={removeLattice}
            keepLattice={keepLattice}
            outlierOverride={outlierOverride}
            expertOptions={expertOptions}
            vis={vis}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
