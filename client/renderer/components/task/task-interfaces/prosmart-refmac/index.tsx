import React, { useCallback } from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { useJob } from "../../../../utils";
import { useFreeRWarning } from "../../../../providers/run-check-provider";
import { InputDataTab } from "./InputDataTab";
import { ParameterisationTab } from "./ParameterisationTab";
import { RestraintsTab } from "./RestraintsTab";
import { OutputTab } from "./OutputTab";
import { AdvancedTab } from "./AdvancedTab";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, fetchDigest, createPeerTask, validation, useFileDigest } =
    useJob(job.id);

  // --- F_SIGF onChange: extract wavelength, set anomalous/twinning flags ---
  const { item: F_SIGFItem, value: F_SIGFValue } = useTaskItem("F_SIGF");
  const { value: freeRFlag } = useTaskItem("FREERFLAG");

  useFreeRWarning({
    job,
    taskName: "prosmart_refmac",
    freeRFlag,
    validation,
    createPeerTask,
  });

  const { forceUpdate: forceUpdateWAVELENGTH } = useTaskItem("WAVELENGTH");
  const { forceUpdate: forceUpdateUSEANOMALOUS } = useTaskItem("USEANOMALOUS");
  const { forceUpdate: forceUpdateUSE_TWIN } = useTaskItem("USE_TWIN");

  // --- Values for visibility conditions ---
  const { value: refinementMode } = useTaskItem("REFINEMENT_MODE");
  const { value: hydrUse } = useTaskItem("HYDR_USE");
  const { value: addWaters } = useTaskItem("ADD_WATERS");
  const { value: useAnomalous } = useTaskItem("USEANOMALOUS");
  const { value: solventMaskType } = useTaskItem("SOLVENT_MASK_TYPE");
  const { value: solventAdvanced } = useTaskItem("SOLVENT_ADVANCED");
  const { value: tlsMode } = useTaskItem("TLSMODE");
  const { value: bfacSetUse } = useTaskItem("BFACSETUSE");
  const { value: weightOpt } = useTaskItem("WEIGHT_OPT");
  const { value: occupancyRefinement } = useTaskItem("OCCUPANCY_REFINEMENT");
  const { value: occupancyGroups } = useTaskItem("OCCUPANCY_GROUPS");
  const { value: occupancyComplete } = useTaskItem("OCCUPANCY_COMPLETE");
  const { value: occupancyIncomplete } = useTaskItem("OCCUPANCY_INCOMPLETE");
  const { value: useNcs } = useTaskItem("USE_NCS");
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
  const { value: platonyzerToggle } = useTaskItem("platonyzer.TOGGLE");
  const { value: mapSharp } = useTaskItem("MAP_SHARP");
  const { value: mapSharpCustom } = useTaskItem("MAP_SHARP_CUSTOM");
  const { value: scatteringFactors } = useTaskItem("SCATTERING_FACTORS");
  const { value: resCustom } = useTaskItem("RES_CUSTOM");
  const { value: hdInitToggle } = useTaskItem("HD_INIT_TOGGLE");

  // Chain composition from XYZIN digest (for ProSMART section visibility)
  const { item: XYZINItem, value: XYZINValue } = useTaskItem("container.inputData.XYZIN");
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

  // --- F_SIGF onChange handler ---
  const handleF_SIGFChange = useCallback(async () => {
    if (!job || job.status !== 1) return;

    if (F_SIGFItem?._objectPath) {
      const digestData = await fetchDigest(F_SIGFItem._objectPath);
      const wavelength = digestData?.wavelengths?.at(-1);
      if (
        wavelength &&
        wavelength > 0 &&
        wavelength < 9 &&
        forceUpdateWAVELENGTH
      ) {
        await forceUpdateWAVELENGTH(wavelength);
      }
    }

    if (F_SIGFValue) {
      try {
        const contentFlag = F_SIGFValue.contentFlag;
        if (![1, 2].includes(contentFlag) && forceUpdateUSEANOMALOUS) {
          await forceUpdateUSEANOMALOUS(false);
        }
        if (![3].includes(contentFlag) && forceUpdateUSE_TWIN) {
          await forceUpdateUSE_TWIN(false);
        }
      } catch (error) {
        console.error("Error processing F_SIGF change:", error);
      }
    }
  }, [
    F_SIGFItem?._objectPath,
    F_SIGFValue,
    job,
    fetchDigest,
    forceUpdateWAVELENGTH,
    forceUpdateUSEANOMALOUS,
    forceUpdateUSE_TWIN,
  ]);

  return (
    <Paper>
      <CCP4i2Tabs>
        {/* TAB 1: INPUT DATA */}
        <CCP4i2Tab label="Input data" key="input">
          <InputDataTab
            {...props}
            handleF_SIGFChange={handleF_SIGFChange}
            refinementMode={refinementMode}
            hydrUse={hydrUse}
            addWaters={addWaters}
            useAnomalous={useAnomalous}
            F_SIGFItem={F_SIGFItem}
          />
        </CCP4i2Tab>

        {/* TAB 2: PARAMETERISATION */}
        <CCP4i2Tab label="Parameterisation" key="parameterisation">
          <ParameterisationTab
            {...props}
            refinementMode={refinementMode}
            solventMaskType={solventMaskType}
            solventAdvanced={solventAdvanced}
            tlsMode={tlsMode}
            bfacSetUse={bfacSetUse}
            weightOpt={weightOpt}
            occupancyRefinement={occupancyRefinement}
            occupancyGroups={occupancyGroups}
            occupancyComplete={occupancyComplete}
            occupancyIncomplete={occupancyIncomplete}
          />
        </CCP4i2Tab>

        {/* TAB 3: RESTRAINTS */}
        <CCP4i2Tab label="Restraints" key="restraints">
          <RestraintsTab
            {...props}
            useNcs={useNcs}
            useJelly={useJelly}
            prosmartProteinToggle={prosmartProteinToggle}
            prosmartProteinAdvanced={prosmartProteinAdvanced}
            prosmartNucleicAcidToggle={prosmartNucleicAcidToggle}
            prosmartNucleicAcidAdvanced={prosmartNucleicAcidAdvanced}
            platonyzerToggle={platonyzerToggle}
            hasProteinChains={hasProteinChains}
            hasNucleotideChains={hasNucleotideChains}
            xyzinComposition={xyzinComposition}
          />
        </CCP4i2Tab>

        {/* TAB 4: OUTPUT */}
        <CCP4i2Tab label="Output" key="output">
          <OutputTab
            {...props}
            mapSharp={mapSharp}
            mapSharpCustom={mapSharpCustom}
          />
        </CCP4i2Tab>

        {/* TAB 5: ADVANCED */}
        <CCP4i2Tab label="Advanced" key="advanced">
          <AdvancedTab
            {...props}
            scatteringFactors={scatteringFactors}
            hydrUse={hydrUse}
            hdInitToggle={hdInitToggle}
            bfacSetUse={bfacSetUse}
            resCustom={resCustom}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
