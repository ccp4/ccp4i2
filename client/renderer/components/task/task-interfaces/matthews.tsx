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
        <Typography variant="subtitle2">Cell parameters taken from reflection data</Typography>
        <CCP4i2TaskElement itemName="HKLIN" {...props} />
        <Typography variant="subtitle2">Calculate molecular weight from...</Typography>
        <CCP4i2TaskElement itemName="ASUIN" {...props} qualifiers={{ guiLabel: "Contents of biological unit" }} />
        <CCP4i2TaskElement itemName="NRES" {...props} qualifiers={{ guiLabel: "Number of residues" }} />
        <CCP4i2TaskElement itemName="MOLWT" {...props} qualifiers={{ guiLabel: "Molecular weight" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
