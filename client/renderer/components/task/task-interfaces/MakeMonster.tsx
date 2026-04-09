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
        <Typography variant="subtitle2">Input reflection data objects</Typography>
        <CCP4i2TaskElement itemName="NObsObjects" {...props} qualifiers={{ guiLabel: "Number to include" }} />
        <CCP4i2TaskElement itemName="Obs_1" {...props} />

        <Typography variant="subtitle2">Input phase data objects</Typography>
        <CCP4i2TaskElement itemName="NPhiObjects" {...props} qualifiers={{ guiLabel: "Number to include" }} />
        <CCP4i2TaskElement itemName="Phi_1" {...props} />

        <Typography variant="subtitle2">Input map coefficients data objects</Typography>
        <CCP4i2TaskElement itemName="NFWTObjects" {...props} qualifiers={{ guiLabel: "Number to include" }} />
        <CCP4i2TaskElement itemName="FWT_1" {...props} />

        <Typography variant="subtitle2">Input difmap coefficients data objects</Typography>
        <CCP4i2TaskElement itemName="NDELFWTObjects" {...props} qualifiers={{ guiLabel: "Number to include" }} />
        <CCP4i2TaskElement itemName="DELFWT_1" {...props} />

        <Typography variant="subtitle2">Input FreeR flags data objects</Typography>
        <CCP4i2TaskElement itemName="NFREERObjects" {...props} qualifiers={{ guiLabel: "Number to include" }} />
        <CCP4i2TaskElement itemName="FREER_1" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
