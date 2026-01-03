import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: INPUT_REPRESENTATION } = useTaskItem("INPUT_REPRESENTATION");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Script control">
          <CCP4i2TaskElement itemName="INPUT_REPRESENTATION" {...props} qualifiers={{ guiLabel: "Using file with representation" }} />
          {(INPUT_REPRESENTATION === 1) && (
            <CCP4i2TaskElement itemName="F_SIGF_AS_IPAIR" {...props} qualifiers={{ guiLabel: "Reflection object containing I+/I-" }} />
          )}
          {(INPUT_REPRESENTATION === 2) && (
            <CCP4i2TaskElement itemName="F_SIGF_AS_FPAIR" {...props} qualifiers={{ guiLabel: "Reflection object containing F+/F-" }} />
          )}
          {(INPUT_REPRESENTATION === 3) && (
            <CCP4i2TaskElement itemName="F_SIGF_AS_IMEAN" {...props} qualifiers={{ guiLabel: "Reflection object containing Imean" }} />
          )}
          <CCP4i2TaskElement itemName="OUTPUT_REPRESENTATION" {...props} qualifiers={{ guiLabel: "Representation of file to generate" }} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;