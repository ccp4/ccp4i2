import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { InlineField } from "../../task-elements/inline-field";
import { useJob } from "../../../../utils";
import { useBoolToggle } from "../../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: MON_1_TYPE } = useTaskItem("MON_1_TYPE");
  const { value: MON_2_TYPE } = useTaskItem("MON_2_TYPE");
  const toggleDelete1 = useBoolToggle(useTaskItem, "TOGGLE_DELETE_1");
  const toggleDelete2 = useBoolToggle(useTaskItem, "TOGGLE_DELETE_2");
  const toggleChange1 = useBoolToggle(useTaskItem, "TOGGLE_CHANGE_1");
  const toggleChange2 = useBoolToggle(useTaskItem, "TOGGLE_CHANGE_2");
  const toggleCharge1 = useBoolToggle(useTaskItem, "TOGGLE_CHARGE_1");
  const toggleCharge2 = useBoolToggle(useTaskItem, "TOGGLE_CHARGE_2");
  const toggleLink = useBoolToggle(useTaskItem, "TOGGLE_LINK");
  const { value: LINK_MODE } = useTaskItem("LINK_MODE");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          {/* --- First monomer to be linked --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "First monomer to be linked" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <InlineField label="Get ligand description from">
                <CCP4i2TaskElement itemName="MON_1_TYPE" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
              {MON_1_TYPE === "CIF" && (
                <CCP4i2TaskElement itemName="DICT_1" {...props} />
              )}
              {MON_1_TYPE === "CIF" && (
                <InlineField label="Residue name">
                  <CCP4i2TaskElement itemName="RES_NAME_1_CIF" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {MON_1_TYPE === "CIF" && (
                <InlineField label="linking atom">
                  <CCP4i2TaskElement itemName="ATOM_NAME_1_CIF" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {MON_1_TYPE === "TLC" && (
                <InlineField label="Residue name">
                  <CCP4i2TaskElement itemName="RES_NAME_1_TLC" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {MON_1_TYPE === "TLC" && (
                <InlineField label="linking atom">
                  <CCP4i2TaskElement itemName="ATOM_NAME_1_TLC" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              <CCP4i2TaskElement itemName="TOGGLE_DELETE_1" {...props} qualifiers={{ guiLabel: "Delete atom" }} />
              {toggleDelete1.value && (
                <CCP4i2TaskElement itemName="DELETE_1_LIST" {...props} />
              )}
              <CCP4i2TaskElement itemName="TOGGLE_CHANGE_1" {...props} qualifiers={{ guiLabel: "Change order of bond" }} />
              {toggleChange1.value && (
                <InlineField label="">
                  <CCP4i2TaskElement itemName="CHANGE_BOND_1_LIST" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {toggleChange1.value && (
                <InlineField label="to">
                  <CCP4i2TaskElement itemName="CHANGE_BOND_1_TYPE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              <CCP4i2TaskElement itemName="TOGGLE_CHARGE_1" {...props} qualifiers={{ guiLabel: "Change the formal charge of an atom" }} />
              {toggleCharge1.value && (
                <InlineField label="">
                  <CCP4i2TaskElement itemName="CHARGE_1_LIST" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {toggleCharge1.value && (
                <InlineField label="to">
                  <CCP4i2TaskElement itemName="CHARGE_1_VALUE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* --- Second monomer to be linked --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Second monomer to be linked" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <InlineField label="Get ligand description from">
                <CCP4i2TaskElement itemName="MON_2_TYPE" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
              {MON_2_TYPE === "CIF" && (
                <CCP4i2TaskElement itemName="DICT_2" {...props} />
              )}
              {MON_2_TYPE === "CIF" && (
                <InlineField label="Residue name">
                  <CCP4i2TaskElement itemName="RES_NAME_2_CIF" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {MON_2_TYPE === "CIF" && (
                <InlineField label="linking atom">
                  <CCP4i2TaskElement itemName="ATOM_NAME_2_CIF" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {MON_2_TYPE === "TLC" && (
                <InlineField label="Residue name">
                  <CCP4i2TaskElement itemName="RES_NAME_2_TLC" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {MON_2_TYPE === "TLC" && (
                <InlineField label="linking atom">
                  <CCP4i2TaskElement itemName="ATOM_NAME_2_TLC" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              <CCP4i2TaskElement itemName="TOGGLE_DELETE_2" {...props} qualifiers={{ guiLabel: "Delete atom" }} />
              {toggleDelete2.value && (
                <CCP4i2TaskElement itemName="DELETE_2_LIST" {...props} />
              )}
              <CCP4i2TaskElement itemName="TOGGLE_CHANGE_2" {...props} qualifiers={{ guiLabel: "Change order of bond" }} />
              {toggleChange2.value && (
                <InlineField label="">
                  <CCP4i2TaskElement itemName="CHANGE_BOND_2_LIST" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {toggleChange2.value && (
                <InlineField label="to">
                  <CCP4i2TaskElement itemName="CHANGE_BOND_2_TYPE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              <CCP4i2TaskElement itemName="TOGGLE_CHARGE_2" {...props} qualifiers={{ guiLabel: "Change the formal charge of an atom" }} />
              {toggleCharge2.value && (
                <InlineField label="">
                  <CCP4i2TaskElement itemName="CHARGE_2_LIST" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              {toggleCharge2.value && (
                <InlineField label="to">
                  <CCP4i2TaskElement itemName="CHARGE_2_VALUE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* --- Bond order --- */}
          <InlineField label="Order of the bond between linked atoms">
            <CCP4i2TaskElement itemName="BOND_ORDER" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>

          {/* --- Apply links to model (optional) --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Apply links to model (optional)" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="TOGGLE_LINK" {...props} qualifiers={{ guiLabel: "Apply links to model" }} />
            {toggleLink.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ initiallyOpen: true }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement itemName="XYZIN" {...props} />
                <InlineField label="Apply links">
                  <CCP4i2TaskElement itemName="LINK_MODE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
                {LINK_MODE === "AUTO" && (
                  <InlineField label="within" hint="times the dictionary value for this bond">
                    <CCP4i2TaskElement itemName="LINK_DISTANCE" {...props} qualifiers={{ guiLabel: " " }} />
                  </InlineField>
                )}
                {LINK_MODE === "MANUAL" && (
                  <InlineField label="Create link between residues:">
                    <CCP4i2TaskElement itemName="MODEL_RES_LIST" {...props} qualifiers={{ guiLabel: " " }} />
                  </InlineField>
                )}
                {LINK_MODE === "MANUAL" && (
                  <InlineField label="Filter list by atom proximity:">
                    <CCP4i2TaskElement itemName="MODEL_LINK_DISTANCE" {...props} qualifiers={{ guiLabel: " " }} />
                  </InlineField>
                )}
              </CCP4i2ContainerElement>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="controlParameters" label="Advanced">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Advanced AceDRG options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="EXTRA_ACEDRG_INSTRUCTIONS" {...props} />
            <CCP4i2TaskElement itemName="EXTRA_ACEDRG_KEYWORDS" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
