import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

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
        <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
          Optionally load electron density map coefficients for display
        </Typography>
        <CCP4i2TaskElement itemName="FPHIIN_LIST" {...props} />
      </CCP4i2ContainerElement>

      {/* Difference density maps */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Difference density maps" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
          Optionally load difference density map coefficients for display
        </Typography>
        <CCP4i2TaskElement itemName="DELFPHIIN_LIST" {...props} />
      </CCP4i2ContainerElement>

      {/* Sequences */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Sequences" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
          Optionally load sequence files
        </Typography>
        <CCP4i2TaskElement itemName="SEQUENCE_LIST" {...props} />
      </CCP4i2ContainerElement>

      {/* Ligand dictionary */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Ligand dictionary" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
          Optionally load a custom ligand dictionary
        </Typography>
        <CCP4i2TaskElement itemName="DICT" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
