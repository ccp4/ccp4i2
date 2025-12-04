import { LinearProgress, Paper } from "@mui/material";
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
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Inputs", initiallyOpen: true }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="MAPIN" {...props} qualifiers={{ toolTip: "Input map coefficients" }} />
            <CCP4i2TaskElement itemName="XYZIN" {...props} qualifiers={{ toolTip: "Input model" }} />
            <CCP4i2TaskElement itemName="FRACTION" {...props} qualifiers={{ toolTip: "Fractionof FC to subtract" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;