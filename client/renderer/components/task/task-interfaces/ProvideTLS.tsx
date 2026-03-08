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
      {/* Optional objects from which to start definition */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Optional objects from which to start definition" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2">Input objects from which to infer (coordinates) or copy (TLS) sets</Typography>
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <CCP4i2TaskElement itemName="TLSIN" {...props} qualifiers={{ guiLabel: "Starting TLS set" }} />
        <CCP4i2TaskElement itemName="TLSTEXT" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
