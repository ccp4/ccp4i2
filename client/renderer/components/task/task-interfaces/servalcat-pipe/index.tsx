import React, { useMemo } from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { useJob } from "../../../../utils";
import { useFreeRWarning } from "../../../../providers/run-check-provider";
import { isTruthy } from "../../task-elements/shared-hooks";
import { InputDataTab } from "./InputDataTab";
import { ParameterisationTab } from "./ParameterisationTab";
import { RestraintsTab } from "./RestraintsTab";
import { AdvancedTab } from "./AdvancedTab";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const {
    useTaskItem,
    createPeerTask,
    validation,
    useFileDigest,
  } = useJob(job.id);

  // --- Reflection data ---
  const { item: HKLINItem, value: HKLINValue } = useTaskItem("HKLIN");
  const { value: freeRFlag } = useTaskItem("FREERFLAG");

  useFreeRWarning({
    job,
    taskName: "servalcat_pipe",
    freeRFlag,
    validation,
    createPeerTask,
  });

  // --- Values for visibility conditions ---
  const { value: dataMethod } = useTaskItem("DATA_METHOD");
  const { value: mergedOrUnmerged } = useTaskItem("MERGED_OR_UNMERGED");
  const { value: hydrUse } = useTaskItem("HYDR_USE");
  const { value: addWaters } = useTaskItem("ADD_WATERS");
  const { value: useAnomalous } = useTaskItem("USEANOMALOUS");
  const { value: weightOpt } = useTaskItem("WEIGHT_OPT");
  const { value: bfacSetUse } = useTaskItem("BFACSETUSE");
  const { value: randomizeUse } = useTaskItem("RANDOMIZEUSE");
  const { value: occupancyRefinement } = useTaskItem("OCCUPANCY_REFINEMENT");
  const { value: occupancyGroups } = useTaskItem("OCCUPANCY_GROUPS");
  const { value: occupancyComplete } = useTaskItem("OCCUPANCY_COMPLETE");
  const { value: occupancyIncomplete } = useTaskItem("OCCUPANCY_INCOMPLETE");
  const { value: useJelly } = useTaskItem("USE_JELLY");
  const { value: prosmartProteinToggle } = useTaskItem(
    "prosmartProtein.TOGGLE"
  );
  const { value: prosmartProteinAdvanced } = useTaskItem(
    "prosmartProtein.ADVANCED"
  );
  const { value: prosmartNucleicAcidToggle } = useTaskItem(
    "prosmartNucleicAcid.TOGGLE"
  );
  const { value: prosmartNucleicAcidAdvanced } = useTaskItem(
    "prosmartNucleicAcid.ADVANCED"
  );
  const { value: libgToggle } = useTaskItem("libg.TOGGLE");
  const { value: libgOption } = useTaskItem("libg.OPTION");
  const { value: libgAdvanced } = useTaskItem("libg.ADVANCED");
  const { value: platonyzerToggle } = useTaskItem("platonyzer.TOGGLE");
  const { value: metalCoordRun } = useTaskItem(
    "metalCoordPipeline.RUN_METALCOORD"
  );
  const { value: metalCoordGenOrUse } = useTaskItem(
    "metalCoordPipeline.GENERATE_OR_USE"
  );
  const { value: metalCoordAdvanced } = useTaskItem(
    "metalCoordPipeline.TOGGLE_ADVANCED"
  );
  const { value: scatteringFactors } = useTaskItem("SCATTERING_FACTORS");
  const { value: resCustom } = useTaskItem("RES_CUSTOM");
  const { value: blurUse } = useTaskItem("BLURUSE");
  const { value: runAdpAnalysis } = useTaskItem("RUN_ADP_ANALYSIS");
  const { value: runCoordAdpDev } = useTaskItem(
    "monitor.RUN_COORDADPDEV_ANALYSIS"
  );

  // Chain composition from XYZIN digest (for ProSMART/libg section visibility)
  const { item: XYZINItem, value: XYZINValue } = useTaskItem("XYZIN");
  const xyzinDigestPath =
    XYZINValue?.dbFileId && XYZINItem?._objectPath
      ? XYZINItem._objectPath
      : "";
  const { data: xyzinDigest } = useFileDigest(xyzinDigestPath);
  const xyzinComposition = xyzinDigest?.composition;
  const hasProteinChains =
    Array.isArray(xyzinComposition?.peptides) &&
    xyzinComposition.peptides.length > 0;
  const hasNucleotideChains =
    Array.isArray(xyzinComposition?.nucleics) &&
    xyzinComposition.nucleics.length > 0;
  const hasMetalSites =
    Array.isArray(xyzinComposition?.ligands) &&
    xyzinComposition.ligands.some(
      (l: any) => l.isMetal || l.type === "metal"
    );

  // Derived state
  const isXtal = dataMethod === "xtal";
  const isSpa = dataMethod === "spa";
  const isMerged = mergedOrUnmerged === "merged";

  const intensitiesAvailable = useMemo(() => {
    return [1, 3].includes(HKLINValue?.contentFlag);
  }, [HKLINValue?.contentFlag]);

  return (
    <Paper>
      <CCP4i2Tabs>
        {/* TAB 1: INPUT DATA */}
        <CCP4i2Tab label="Input Data" key="input">
          <InputDataTab
            {...props}
            isXtal={isXtal}
            isSpa={isSpa}
            isMerged={isMerged}
            intensitiesAvailable={intensitiesAvailable}
            hydrUse={hydrUse}
            addWaters={addWaters}
            useAnomalous={useAnomalous}
            HKLINValue={HKLINValue}
          />
        </CCP4i2Tab>

        {/* TAB 2: PARAMETERISATION */}
        <CCP4i2Tab label="Parameterisation" key="parameterisation">
          <ParameterisationTab
            {...props}
            occupancyGroups={occupancyGroups}
            occupancyComplete={occupancyComplete}
            occupancyIncomplete={occupancyIncomplete}
            occupancyRefinement={occupancyRefinement}
          />
        </CCP4i2Tab>

        {/* TAB 3: RESTRAINTS */}
        <CCP4i2Tab label="Restraints" key="restraints">
          <RestraintsTab
            {...props}
            weightOpt={weightOpt}
            useJelly={useJelly}
            prosmartProteinToggle={prosmartProteinToggle}
            prosmartProteinAdvanced={prosmartProteinAdvanced}
            prosmartNucleicAcidToggle={prosmartNucleicAcidToggle}
            prosmartNucleicAcidAdvanced={prosmartNucleicAcidAdvanced}
            libgToggle={libgToggle}
            libgOption={libgOption}
            libgAdvanced={libgAdvanced}
            platonyzerToggle={platonyzerToggle}
            metalCoordRun={metalCoordRun}
            metalCoordGenOrUse={metalCoordGenOrUse}
            metalCoordAdvanced={metalCoordAdvanced}
            hasProteinChains={hasProteinChains}
            hasNucleotideChains={hasNucleotideChains}
            hasMetalSites={hasMetalSites}
            xyzinComposition={xyzinComposition}
          />
        </CCP4i2Tab>

        {/* TAB 4: ADVANCED */}
        <CCP4i2Tab label="Advanced" key="advanced">
          <AdvancedTab
            {...props}
            isXtal={isXtal}
            isSpa={isSpa}
            bfacSetUse={bfacSetUse}
            randomizeUse={randomizeUse}
            scatteringFactors={scatteringFactors}
            resCustom={resCustom}
            blurUse={blurUse}
            runAdpAnalysis={runAdpAnalysis}
            runCoordAdpDev={runCoordAdpDev}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
