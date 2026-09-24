import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Map file to import" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="MAPIN" {...props} />
        {/* Let the user say what kind of map this is. Without this, every map
            imported through the desktop UI was silently subType 1 (normal), so
            half maps and masks were mis-typed (#524). The wrapper honours
            MAP_SUBTYPE; it just was never shown. */}
        <CCP4i2TaskElement itemName="MAP_SUBTYPE" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
