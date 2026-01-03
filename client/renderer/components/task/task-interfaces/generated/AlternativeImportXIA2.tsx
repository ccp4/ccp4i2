import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Xia2 runs">
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Xia2 directory to import
          </Typography>
          <CCP4i2TaskElement itemName="XIA2_DIRECTORY" {...props} qualifiers={{ toolTip: "Browse to the directory \"xia2\"" }} />
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            When a valid top-level XIA2 directory is selected above, the list of successful
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            data reduction protocols that XIA2 performed will be summarised below.
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Select and delete ("-") the protocols you do not wish to import
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Protocols to keep from this XIA2 directory
          </Typography>
          <CCP4i2TaskElement itemName="runSummaries" {...props} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;