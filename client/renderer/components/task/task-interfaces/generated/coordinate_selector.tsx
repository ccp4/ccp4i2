import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="XYZIN"
          {...props}
          qualifiers={{ toolTip: "Input model from which subset will be selected" }}
        />
        <CCP4i2TaskElement
          itemName="OVERRIDE_SUBTYPE"
          {...props}
          qualifiers={{
            guiLabel: "These coordinates should be used as",
            toolTip: "Specify how coordinates should be used in other tasks",
          }}
        />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
