import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Input coordinates */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input coordinates" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="XYZIN"
          {...props}
          qualifiers={{ guiLabel: "Atomic model" }}
        />
      </CCP4i2ContainerElement>

      {/* Input reflections */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input reflections" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="F_SIGF"
          {...props}
          qualifiers={{ guiLabel: "Reflections" }}
        />
        <CCP4i2TaskElement
          itemName="FREERFLAG"
          {...props}
          qualifiers={{ guiLabel: "Free R set" }}
        />
      </CCP4i2ContainerElement>

      {/* Control parameters */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Control parameters" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="MR_WHEN_R"
          {...props}
          qualifiers={{ guiLabel: "R-factor threshold to force MR" }}
        />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
