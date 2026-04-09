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
          Assign two or more datasets in the list: the first one will be treated as &quot;native&quot;
        </Typography>
        <CCP4i2TaskElement itemName="MERGEDFILES" {...props} />
        <CCP4i2TaskElement itemName="RESOLUTION_MAX" {...props} qualifiers={{ guiLabel: "Maximum resolution used (\u00C5)" }} />
        <Typography variant="caption" color="text.secondary">
          SCALEIT will compare datasets to the first one (treated as &quot;native&quot;).
          Other &quot;derivative&quot; datasets will be scaled with a scale and anisotropic B-factor.
        </Typography>
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
