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
      {/* Input Data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input Data" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2">
          Select either the &apos;top&apos; XIA2 directory to import all processing runs or select one processing run
        </Typography>
        <Typography variant="body2">
          Each processing run in the XIA2 directory will be automatically imported as a separate job
        </Typography>
        <CCP4i2TaskElement itemName="XIA2_DIRECTORY" {...props} qualifiers={{ guiLabel: "XIA2 Directory" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
