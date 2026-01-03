import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: MAXIMUM_COORDINATION_NUMBER } = useTaskItem("MAXIMUM_COORDINATION_NUMBER");
  const { value: SAVE_PDBMMCIF } = useTaskItem("SAVE_PDBMMCIF");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <CCP4i2TaskElement itemName="XYZIN" {...props} qualifiers={{ toolTip: "Atomic model" }} />
          <CCP4i2TaskElement itemName="LIGAND_CODE" {...props} qualifiers={{ guiLabel: "Monomer code" }} />
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Advanced parameters", initiallyOpen: true }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="DISTANCE_THRESHOLD" {...props} qualifiers={{ guiLabel: "Distance threshold: (range 0-1) A threshold d to select atoms is (r1 + r2)*(1 + d) where r1 and r2 are covalent radii." }} />
            <CCP4i2TaskElement itemName="MAXIMUM_COORDINATION_NUMBER" {...props} qualifiers={{ guiLabel: "Maximum coordination number:" }} />
            {(MAXIMUM_COORDINATION_NUMBER === "2") && (
              <CCP4i2TaskElement itemName="COORD02" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "3") && (
              <CCP4i2TaskElement itemName="COORD03" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "4") && (
              <CCP4i2TaskElement itemName="COORD04" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "5") && (
              <CCP4i2TaskElement itemName="COORD05" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "6") && (
              <CCP4i2TaskElement itemName="COORD06" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "7") && (
              <CCP4i2TaskElement itemName="COORD07" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "8") && (
              <CCP4i2TaskElement itemName="COORD08" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "9") && (
              <CCP4i2TaskElement itemName="COORD09" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "10") && (
              <CCP4i2TaskElement itemName="COORD10" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "11") && (
              <CCP4i2TaskElement itemName="COORD11" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "12") && (
              <CCP4i2TaskElement itemName="COORD12" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "13") && (
              <CCP4i2TaskElement itemName="COORD13" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "14") && (
              <CCP4i2TaskElement itemName="COORD14" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "15") && (
              <CCP4i2TaskElement itemName="COORD15" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "16") && (
              <CCP4i2TaskElement itemName="COORD16" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "17") && (
              <CCP4i2TaskElement itemName="COORD17" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            {(MAXIMUM_COORDINATION_NUMBER === "24") && (
              <CCP4i2TaskElement itemName="COORD24" {...props} qualifiers={{ guiLabel: "Coordination class:" }} />
            )}
            <CCP4i2TaskElement itemName="PROCRUSTES_DISTANCE_THRESHOLD" {...props} qualifiers={{ guiLabel: "Procrustes distance threshold: (range 0-1)" }} />
            <CCP4i2TaskElement itemName="MINIMUM_SAMPLE_SIZE" {...props} qualifiers={{ guiLabel: "Minimum sample size for statistics:" }} />
            <CCP4i2TaskElement itemName="USE_PDB" {...props} qualifiers={{ guiLabel: "Use COD structures based on the input PDB/mmCIF coordinates" }} />
            <CCP4i2TaskElement itemName="IDEAL_ANGLES" {...props} qualifiers={{ guiLabel: "Provide only ideal bond angles" }} />
            <CCP4i2TaskElement itemName="SIMPLE" {...props} qualifiers={{ guiLabel: "Simple distance based filtering" }} />
            <CCP4i2TaskElement itemName="SAVE_PDBMMCIF" {...props} qualifiers={{ guiLabel: "Update link records to metal sites in the atomic model" }} />
            {(SAVE_PDBMMCIF === true) && (
              <CCP4i2TaskElement itemName="KEEP_LINKS" {...props} qualifiers={{ guiLabel: "Keep existing link records to metal sites" }} />
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;