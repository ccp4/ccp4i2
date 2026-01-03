import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: TOGGLE_LINK } = useTaskItem("TOGGLE_LINK");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            First monomer to be linked
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="MON_1_TYPE" {...props} qualifiers={{ guiLabel: "Get ligand description from" }} />
            <CCP4i2TaskElement itemName="DICT_1" {...props} />
            <CCP4i2TaskElement itemName="ATOM_NAME_1_CIF" {...props} qualifiers={{ guiLabel: "linking atom" }} />
            <CCP4i2TaskElement itemName="ATOM_NAME_1_TLC" {...props} qualifiers={{ guiLabel: "linking atom" }} />
            <CCP4i2TaskElement itemName="TOGGLE_DELETE_1" {...props} qualifiers={{ guiLabel: "Delete atom" }} />
            <CCP4i2TaskElement itemName="DELETE_1_LIST" {...props} />
            <CCP4i2TaskElement itemName="TOGGLE_CHANGE_1" {...props} qualifiers={{ guiLabel: "Change order of bond" }} />
            <CCP4i2TaskElement itemName="CHANGE_BOND_1_TYPE" {...props} qualifiers={{ guiLabel: "to" }} />
            <CCP4i2TaskElement itemName="TOGGLE_CHARGE_1" {...props} qualifiers={{ guiLabel: "Change the formal charge of an atom" }} />
            <CCP4i2TaskElement itemName="CHARGE_1_VALUE" {...props} qualifiers={{ guiLabel: "to" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Second monomer to be linked
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="MON_2_TYPE" {...props} qualifiers={{ guiLabel: "Get ligand description from" }} />
            <CCP4i2TaskElement itemName="DICT_2" {...props} />
            <CCP4i2TaskElement itemName="ATOM_NAME_2_CIF" {...props} qualifiers={{ guiLabel: "linking atom" }} />
            <CCP4i2TaskElement itemName="ATOM_NAME_2_TLC" {...props} qualifiers={{ guiLabel: "linking atom" }} />
            <CCP4i2TaskElement itemName="TOGGLE_DELETE_2" {...props} qualifiers={{ guiLabel: "Delete atom" }} />
            <CCP4i2TaskElement itemName="DELETE_2_LIST" {...props} />
            <CCP4i2TaskElement itemName="TOGGLE_CHANGE_2" {...props} qualifiers={{ guiLabel: "Change order of bond" }} />
            <CCP4i2TaskElement itemName="CHANGE_BOND_2_TYPE" {...props} qualifiers={{ guiLabel: "to" }} />
            <CCP4i2TaskElement itemName="TOGGLE_CHARGE_2" {...props} qualifiers={{ guiLabel: "Change the formal charge of an atom" }} />
            <CCP4i2TaskElement itemName="CHARGE_2_VALUE" {...props} qualifiers={{ guiLabel: "to" }} />
          </CCP4i2ContainerElement>
          <CCP4i2TaskElement itemName="BOND_ORDER" {...props} qualifiers={{ guiLabel: "Order of the bond between linked atoms" }} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Apply links to model (optional)
          </Typography>
          <CCP4i2TaskElement itemName="TOGGLE_LINK" {...props} />
          {(TOGGLE_LINK === true) && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              {(TOGGLE_LINK === true) && (
                <CCP4i2TaskElement itemName="XYZIN" {...props} />
              )}
              <CCP4i2TaskElement itemName="LINK_MODE" {...props} qualifiers={{ guiLabel: "Apply links" }} />
              <CCP4i2TaskElement itemName="LINK_DISTANCE" {...props} qualifiers={{ guiLabel: "times the dictionary value for this bond" }} />
              <CCP4i2TaskElement itemName="MODEL_LINK_DISTANCE" {...props} qualifiers={{ guiLabel: "Filter list by atom proximity:" }} />
            </CCP4i2ContainerElement>
          )}
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Advanced">
          <CCP4i2TaskElement itemName="EXTRA_ACEDRG_INSTRUCTIONS" {...props} />
          <CCP4i2TaskElement itemName="EXTRA_ACEDRG_KEYWORDS" {...props} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;