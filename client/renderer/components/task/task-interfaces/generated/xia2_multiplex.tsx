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
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Search for DIALS integrated files
          </Typography>
          <CCP4i2TaskElement itemName="SEARCH_ROOT_DIR" {...props} />
          <CCP4i2TaskElement itemName="DIALS_INTEGRATED" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Basic parameters
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Clustering
          </Typography>
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Advanced parameters">
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;