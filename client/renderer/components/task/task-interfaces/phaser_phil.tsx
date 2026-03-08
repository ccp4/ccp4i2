import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          {/* --- Input coordinates --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input coordinates" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="XYZIN" {...props} />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* --- Input reflections --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input reflections" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* --- Basic parameters (nestedAutoGenerate at expertLevel 0) --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Basic parameters" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="phaser__mode" {...props} />
              <CCP4i2TaskElement itemName="phaser__hklin__labin" {...props} />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="controlParameters" label="Advanced parameters">
          {/* Advanced parameters (nestedAutoGenerate at expertLevel 0+1) */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Advanced parameters" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="phaser__mode" {...props} />
            <CCP4i2TaskElement itemName="phaser__hklin__labin" {...props} />
            <CCP4i2TaskElement itemName="phaser__keywords__composition__by" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
