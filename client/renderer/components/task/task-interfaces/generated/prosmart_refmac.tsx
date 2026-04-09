import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: ADD_WATERS } = useTaskItem("ADD_WATERS");
  const { value: BFACSETUSE } = useTaskItem("BFACSETUSE");
  const { value: HD_INIT_TOGGLE } = useTaskItem("HD_INIT_TOGGLE");
  const { value: HYDR_USE } = useTaskItem("HYDR_USE");
  const { value: RES_CUSTOM } = useTaskItem("RES_CUSTOM");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Main inputs
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2TaskElement itemName="F_SIGF" {...props} />
          <CCP4i2TaskElement itemName="WAVELENGTH" {...props} qualifiers={{ guiLabel: "Wavelength" }} />
          <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Experimental phase information
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="ABCD" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Additional geometry dictionaries
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="DICT_LIST" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="REFINEMENT_MODE" {...props} qualifiers={{ guiLabel: "Refinement mode (restrained or rigid body):" }} />
            <CCP4i2TaskElement itemName="NCYCLES" {...props} qualifiers={{ guiLabel: "Number of restrained refinement cycles:" }} />
            <CCP4i2TaskElement itemName="NCYCRIGID" {...props} qualifiers={{ guiLabel: "Number of rigid body refinement cycles:" }} />
            {(HYDR_USE === false) && (
              <CCP4i2TaskElement itemName="HYDR_USE" {...props} qualifiers={{ guiLabel: "Use riding hydrogens during refinement" }} />
            )}
            {(HYDR_USE === true) && (
              <CCP4i2TaskElement itemName="HYDR_ALL" {...props} qualifiers={{ guiLabel: "Use riding hydrogens during refinement" }} />
            )}
            <CCP4i2TaskElement itemName="ADD_WATERS" {...props} qualifiers={{ guiLabel: "Add waters" }} />
            {(ADD_WATERS === true) && (
              <CCP4i2TaskElement itemName="REFPRO_RSR_RWORK_LIMIT" {...props} qualifiers={{ guiLabel: "or lower" }} />
            )}
            <CCP4i2TaskElement itemName="USE_TWIN" {...props} qualifiers={{ guiLabel: "Crystal is twinned" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Parameterisation">
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Restraints">
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Output">
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Advanced">
          <CCP4i2TaskElement itemName="SCATTERING_FACTORS" {...props} qualifiers={{ guiLabel: "Diffraction experiment type:" }} />
          <CCP4i2TaskElement itemName="SCATTERING_ELECTRON" {...props} qualifiers={{ guiLabel: "indentindentForm factor calculation:" }} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Neutron refinement options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            {(HYDR_USE === false) && (
              <CCP4i2TaskElement itemName="HYDR_USE" {...props} qualifiers={{ guiLabel: "Use hydrogens during refinement" }} />
            )}
            {(HYDR_USE === true) && (
              <CCP4i2TaskElement itemName="HYDR_ALL" {...props} qualifiers={{ guiLabel: "Use hydrogens during refinement" }} />
            )}
            {(HYDR_USE === true) && (
              <CCP4i2TaskElement itemName="HD_INIT_TOGGLE" {...props} qualifiers={{ guiLabel: "Initialise hydrogen/deuterium fractions" }} />
            )}
            {(HD_INIT_TOGGLE === true) && (
              <CCP4i2TaskElement itemName="HD_INIT" {...props} />
            )}
            {(HYDR_USE === true) && (
              <CCP4i2TaskElement itemName="HD_FRACTION" {...props} qualifiers={{ guiLabel: "Refine hydrogen/deuterium fractions" }} />
            )}
            <CCP4i2TaskElement itemName="HD_FRACTION_TYPE" {...props} qualifiers={{ guiLabel: "for" }} />
            {(HYDR_USE === true) && (
              <CCP4i2TaskElement itemName="H_REFINE" {...props} qualifiers={{ guiLabel: "Refine hydrogen positions" }} />
            )}
            <CCP4i2TaskElement itemName="H_REFINE_SELECT" {...props} qualifiers={{ guiLabel: "for" }} />
            {(HYDR_USE === true) && (
              <CCP4i2TaskElement itemName="H_TORSION" {...props} qualifiers={{ guiLabel: "Use hydrogen torsion angle restraints" }} />
            )}
          </CCP4i2ContainerElement>
          <CCP4i2TaskElement itemName="RES_CUSTOM" {...props} qualifiers={{ guiLabel: "Use custom resolution limits" }} />
          {(RES_CUSTOM === true) && (
            <CCP4i2TaskElement itemName="RES_MAX" {...props} qualifiers={{ guiLabel: "high (dmin):" }} />
          )}
          <CCP4i2TaskElement itemName="BFACSETUSE" {...props} qualifiers={{ guiLabel: "Reset all B-factors at start" }} />
          {(BFACSETUSE === true) && (
            <CCP4i2TaskElement itemName="BFACSET" {...props} qualifiers={{ guiLabel: "&nbsp;to fixed value:" }} />
          )}
          <CCP4i2TaskElement itemName="ASUIN" {...props} />
          <CCP4i2TaskElement itemName="REFMAC_CLEANUP" {...props} qualifiers={{ guiLabel: "Clean up intermediate files at end of job" }} />
          <CCP4i2TaskElement itemName="EXTRAREFMACKEYWORDS" {...props} />
          <CCP4i2TaskElement itemName="REFMAC_KEYWORD_FILE" {...props} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;