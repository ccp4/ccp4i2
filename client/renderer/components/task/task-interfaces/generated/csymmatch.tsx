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
        <CCP4i2Tab key="inputData" label="Input Data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Atomic model to be moved
          </Typography>
          <CCP4i2TaskElement itemName="XYZIN_QUERY" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Fixed (reference) atomic model
          </Typography>
          <CCP4i2TaskElement itemName="XYZIN_TARGET" {...props} />
          <CCP4i2TaskElement itemName="ORIGIN_HAND" {...props} qualifiers={{ guiLabel: "Try all possible origins and hands" }} />
          <CCP4i2TaskElement itemName="CONNECTIVITY_RADIUS" {...props} qualifiers={{ guiLabel: "Radius to use in stiching floating fragments to chains" }} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;