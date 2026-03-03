import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Reflection data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="HKLIN" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Search models" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
            <CCP4i2TaskElement itemName="ASUIN" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Additional input", initiallyOpen: false }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="DICT" {...props} />
            <CCP4i2TaskElement itemName="DAGIN" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="keywords" label="Keywords">
          <CCP4i2TaskElement
            itemName="controlParameters"
            {...props}
            qualifiers={{ guiLabel: "PhaserTNG Parameters" }}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
