import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Reflection data - mean intensities or anomalous intensities
          </Typography>
          <CCP4i2TaskElement itemName="OBSIN" {...props} />
        </CCP4i2Tab>
        <CCP4i2Tab key="folder_1" label="Options">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            AU content as sequence or number of residues
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="SEQIN" {...props} />
            <CCP4i2TaskElement itemName="NRES" {...props} qualifiers={{ guiLabel: "OR - number of residues in the asymmetric unit" }} />
          </CCP4i2ContainerElement>
          <CCP4i2TaskElement itemName="NO_ANISO" {...props} qualifiers={{ guiLabel: "Do not perform anisotropy correction" }} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;