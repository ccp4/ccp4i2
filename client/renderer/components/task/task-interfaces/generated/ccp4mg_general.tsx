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
          <CCP4i2TaskElement itemName="XYZIN_LIST" {...props} />
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Electron density maps:
          </Typography>
          <CCP4i2TaskElement itemName="FPHIIN_LIST" {...props} />
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Difference density maps:
          </Typography>
          <CCP4i2TaskElement itemName="DELFPHIIN_LIST" {...props} />
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Sequences:
          </Typography>
          <CCP4i2TaskElement itemName="SEQUENCE_LIST" {...props} />
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Ligand dictionary:
          </Typography>
          <CCP4i2TaskElement itemName="DICT" {...props} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;