import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: PASTEORREAD } = useTaskItem("PASTEORREAD");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Optional objects from which to start definition">
          <CCP4i2TaskElement itemName="PASTEORREAD" {...props} qualifiers={{ guiLabel: "Paste or read alignment, or extract from HHPred or Blast search" }} />
          {(PASTEORREAD === "ALIGNIN") && (
            <CCP4i2TaskElement itemName="ALIGNIN" {...props} qualifiers={{ toolTip: "Alignment object or file" }} />
          )}
          {(PASTEORREAD === "HHPREDIN") && (
            <CCP4i2TaskElement itemName="HHPREDIN" {...props} qualifiers={{ toolTip: "HHPred results" }} />
          )}
          {(PASTEORREAD === "BLASTIN") && (
            <CCP4i2TaskElement itemName="BLASTIN" {...props} qualifiers={{ toolTip: "Blast results" }} />
          )}
          {(PASTEORREAD === "PASTE") && (
            <CCP4i2TaskElement itemName="SEQUENCETEXT" {...props} />
          )}
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Choose one alignment from the search file
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Annotation for the alignment
          </Typography>
          <CCP4i2TaskElement itemName="ANNOTATION" {...props} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;