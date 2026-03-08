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
        <Typography variant="subtitle2">Enter structure to be edited</Typography>
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <Typography variant="subtitle2">Enter sequence alignment of edited structure and MR target</Typography>
        <CCP4i2TaskElement itemName="ALIGNIN" {...props} />
        <CCP4i2TaskElement itemName="TARGETINDEX" {...props} qualifiers={{ guiLabel: "Identifier of the target sequence in this alignemnt" }} />
      </CCP4i2ContainerElement>

      {/* Simple Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Simple Options" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2">How severe is side-chain truncation for non-conserved residues:</Typography>
        <CCP4i2TaskElement itemName="MODE" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
