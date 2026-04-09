import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Coordinates */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Coordinates" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="XYZIN_LIST" {...props} />
      </CCP4i2ContainerElement>

      {/* Electron density maps */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Electron density maps" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="FPHIIN_LIST" {...props} />
      </CCP4i2ContainerElement>

      {/* Difference density maps */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Difference density maps" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="DELFPHIIN_LIST" {...props} />
      </CCP4i2ContainerElement>

      {/* Anomalous difference maps */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Anomalous difference maps" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="DELFPHIINANOM_LIST" {...props} />
      </CCP4i2ContainerElement>

      {/* Additional data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Additional data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="DICT" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
