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
        <CCP4i2Tab key="inputData" label="Enter sequences">
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Optionally load existing AU content file to edit
          </Typography>
          <CCP4i2TaskElement itemName="ASUCONTENTIN" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Specify the protein/nucleic acid sequences in the crystal
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Solvent content analysis will be done if you provide experimental data. Run the task to commit to database
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="ASU_CONTENT" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Solvent content analysis
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="HKLIN" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;