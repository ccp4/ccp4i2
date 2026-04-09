import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: ARCIMBOLDO_OPTIONS } = useTaskItem("ARCIMBOLDO_OPTIONS");
  const { value: ARCIMBOLDO_RUN } = useTaskItem("ARCIMBOLDO_RUN");
  const { value: BORGES_LIBRARY } = useTaskItem("BORGES_LIBRARY");
  const { value: DEVELOPER_MODE } = useTaskItem("DEVELOPER_MODE");
  const { value: LITE_MODELS } = useTaskItem("LITE_MODELS");
  const { value: LITE_PARTIAL } = useTaskItem("LITE_PARTIAL");
  const { value: SHREDDER_OPTIONS } = useTaskItem("SHREDDER_OPTIONS");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <CCP4i2TaskElement itemName="ARCIMBOLDO_RUN" {...props} qualifiers={{ guiLabel: "on", toolTip: "Define where to run ARCIMBOLDO" }} />
          <CCP4i2TaskElement itemName="COIL_COILED" {...props} qualifiers={{ guiLabel: "Run in coil coiled mode", toolTip: "Run in coil coiled mode" }} />
          <CCP4i2TaskElement itemName="SHREDDER_PREDICTED" {...props} qualifiers={{ guiLabel: "Run in predicted model mode", toolTip: "Run in predicted model mode" }} />
          {(ARCIMBOLDO_RUN === "local_grid" || ARCIMBOLDO_RUN === "remote_grid") && (
            <CCP4i2TaskElement itemName="RUN_MODE" {...props} qualifiers={{ guiLabel: "Grid configuration" }} />
          )}
          <CCP4i2TaskElement itemName="CONFIG_FILE" {...props} qualifiers={{ guiLabel: "Local path of the configuration file", toolTip: "This file will be referred to in the job's bor-file in the instruction setup_bor_path = ..." }} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Input data
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="MOLECULAR_WEIGHT" {...props} qualifiers={{ guiLabel: "Daltons", toolTip: "Composition of the asymmetric unit" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Model
          </Typography>
          {(ARCIMBOLDO_OPTIONS === "LITE") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="LITE_RMSD" {...props} qualifiers={{ guiLabel: "A", toolTip: "Define type of search models and its expected r.m.s.d. from target" }} />
              {(LITE_MODELS === "HELIX") && (
                <CCP4i2TaskElement itemName="HELIX_LENGTH" {...props} qualifiers={{ guiLabel: "residues", toolTip: "Define the number of helical fragments to search for and the length of the fragment" }} />
              )}
              {(LITE_MODELS === "CUSTOM") && (
                <CCP4i2TaskElement itemName="N_FRAGMENTS" {...props} qualifiers={{ guiLabel: "copies of custom model", toolTip: "Define the number of custom fragments to search for" }} />
              )}
              {(LITE_MODELS === "CUSTOM") && (
                <CCP4i2TaskElement itemName="PDB_LITE" {...props} />
              )}
              {(LITE_MODELS === "HELICES") && (
                <CCP4i2TaskElement itemName="LITE_HELICES_LIST" {...props} qualifiers={{ toolTip: "Define length of an Ideal helical fragment" }} />
              )}
              {(LITE_MODELS === "CUSTOMS") && (
                <CCP4i2TaskElement itemName="LITE_CUSTOMS_LIST" {...props} />
              )}
              <CCP4i2TaskElement itemName="LITE_PARTIAL" {...props} qualifiers={{ guiLabel: "Start from known partial structure", toolTip: "Start from known partial structure" }} />
              {(LITE_PARTIAL === true) && (
                <CCP4i2TaskElement itemName="LITE_FIXED" {...props} qualifiers={{ guiLabel: "Fixed in" }} />
              )}
            </CCP4i2ContainerElement>
          )}
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Borges libraries and Model Handling Parameters
          </Typography>
          {(ARCIMBOLDO_OPTIONS === "BORGES") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="BORGES_LIBRARY" {...props} qualifiers={{ guiLabel: "(topology: u - up, d - down)", toolTip: "Orientations of secondary structure elements are shown as u (up) and d (down)" }} />
              {(BORGES_LIBRARY === "CUSTOM") && (
                <CCP4i2TaskElement itemName="BORGES_CUSTOM" {...props} qualifiers={{ guiLabel: "Library in" }} />
              )}
              <CCP4i2TaskElement itemName="BORGES_GYRE_T" {...props} qualifiers={{ guiLabel: "Switch Phaser GYRE option", toolTip: "Use GYRE option when running Phaser's Rotation Function" }} />
              <CCP4i2TaskElement itemName="BORGES_GIMBLE_T" {...props} qualifiers={{ guiLabel: "Switch Phaser GIMBLE option", toolTip: "Use GIMBLE option when running Phaser's Rigid Body Refinement" }} />
              <CCP4i2TaskElement itemName="BORGES_MULTICOPY_T" {...props} qualifiers={{ guiLabel: "Switch MULTICOPY option", toolTip: "Use MULTICOPY mode" }} />
            </CCP4i2ContainerElement>
          )}
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Shredder models
          </Typography>
          {(ARCIMBOLDO_OPTIONS === "SHREDDER") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="PDB_SHREDDER" {...props} />
              <CCP4i2TaskElement itemName="SHREDDER_RMSD_T" {...props} qualifiers={{ guiLabel: "A" }} />
              <CCP4i2TaskElement itemName="SHREDDER_CONVERT_T" {...props} qualifiers={{ guiLabel: "Convert to polyalanine", toolTip: "Convert input model to polyalanine" }} />
              <CCP4i2TaskElement itemName="SHREDDER_MAKE_T" {...props} qualifiers={{ guiLabel: "Make all B-factors equal", toolTip: "Make all B-factors equal" }} />
              <CCP4i2TaskElement itemName="SHREDDER_OPTIONS" {...props} qualifiers={{ guiLabel: "Shredder mode", toolTip: "Define Shredder Mode" }} />
              {(SHREDDER_OPTIONS === "spherical") && (
                <CCP4i2TaskElement itemName="SHREDDER_COIL_T" {...props} qualifiers={{ guiLabel: "Maintain coil in the model", toolTip: "Maintain coil in the model" }} />
              )}
              {(SHREDDER_OPTIONS === "spherical") && (
                <CCP4i2TaskElement itemName="SHREDDER_GYRE_T" {...props} qualifiers={{ guiLabel: "Perform gyre refinement", toolTip: "Use GYRE option when running Phaser's Rotation Function" }} />
              )}
              {(SHREDDER_OPTIONS === "spherical") && (
                <CCP4i2TaskElement itemName="SHREDDER_GIMBLE_T" {...props} qualifiers={{ guiLabel: "Perform gimble refinement", toolTip: "Use GIMBLE option when running Phaser's Rigid Body Refinement" }} />
              )}
              {(SHREDDER_OPTIONS === "spherical") && (
                <CCP4i2TaskElement itemName="SHREDDER_LLG_T" {...props} qualifiers={{ guiLabel: "Perform LLG-guided pruning", toolTip: "Perform Phaser's LLG-guided pruning" }} />
              )}
              {(SHREDDER_OPTIONS === "spherical") && (
                <CCP4i2TaskElement itemName="SHREDDER_COMBINE_T" {...props} qualifiers={{ guiLabel: "Combine phases with alixe", toolTip: "Combine partial solutions using ALIXE - better phases but runs longer" }} />
              )}
              {(SHREDDER_OPTIONS === "spherical") && (
                <CCP4i2TaskElement itemName="SHREDDER_MULTICOPY_T" {...props} qualifiers={{ guiLabel: "Switch MULTICOPY option", toolTip: "Use MULTICOPY mode" }} />
              )}
            </CCP4i2ContainerElement>
          )}
        </CCP4i2Tab>
        <CCP4i2Tab key="advancedData" label="Advanced data">
          <CCP4i2TaskElement itemName="FRAGMENT_SIZE_T" {...props} qualifiers={{ guiLabel: "Fragment size", toolTip: "Number of amino acids in each fragment" }} />
          <CCP4i2TaskElement itemName="TNCS_T" {...props} qualifiers={{ guiLabel: "Switch Phaser TNCS option", toolTip: "Use TNCS option when running Phaser" }} />
          <CCP4i2TaskElement itemName="SHELXE_LINE" {...props} qualifiers={{ guiLabel: "shelxe_line =", toolTip: "Replace default SHELXE options" }} />
          <CCP4i2TaskElement itemName="KEYWORDS" {...props} qualifiers={{ toolTip: "A line key = value to be added to bor-file" }} />
        </CCP4i2Tab>
        <CCP4i2Tab key="developerOptions" label="Developer options">
          <CCP4i2TaskElement itemName="DEVELOPER_MODE" {...props} qualifiers={{ guiLabel: "Select run mode", toolTip: "Run mode" }} />
          {(DEVELOPER_MODE === "EXISTING") && (
            <CCP4i2TaskElement itemName="EXISTING_FOLDER" {...props} qualifiers={{ guiLabel: "Existing run" }} />
          )}
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;