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
        <Typography variant="subtitle2">Model to analyse</Typography>
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <Typography variant="subtitle2">Map coefficients</Typography>
        <Typography variant="body2">Electron density map</Typography>
        <CCP4i2TaskElement itemName="FPHIIN1" {...props} />
        <Typography variant="body2">Difference density map</Typography>
        <CCP4i2TaskElement itemName="FPHIIN2" {...props} />
      </CCP4i2ContainerElement>

      {/* Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2">Specify the resolution range you would like to use</Typography>
        <CCP4i2TaskElement itemName="RES_LOW" {...props} qualifiers={{ guiLabel: "Low resolution limit" }} />
        <CCP4i2TaskElement itemName="RES_HIGH" {...props} qualifiers={{ guiLabel: "High resolution limit" }} />
        <Typography variant="body2">Map values are averaged separately for main- and side-chains</Typography>
        <CCP4i2TaskElement itemName="MAIN_AVERAGING" {...props} qualifiers={{ guiLabel: "Main-chain averaging is to be performed across" }} />
        <CCP4i2TaskElement itemName="SIDE_AVERAGING" {...props} qualifiers={{ guiLabel: "Side-chain averaging is to be performed across" }} />
        <CCP4i2TaskElement itemName="SCALING" {...props} />
        <CCP4i2TaskElement itemName="SCALING_TYPE" {...props} qualifiers={{ guiLabel: "Rescale Q-Q plot using" }} />
        <Typography variant="body2">Adjust the rejection criteria for flagging up the outliers in Coot</Typography>
        <CCP4i2TaskElement itemName="SIGMA_RZ_MINUS" {...props} qualifiers={{ guiLabel: "Accuracy metrics: RZ- <" }} />
        <CCP4i2TaskElement itemName="SIGMA_RZ_PLUS" {...props} qualifiers={{ guiLabel: "RZ+ >" }} />
        <CCP4i2TaskElement itemName="SIGMA_RO" {...props} qualifiers={{ guiLabel: "Precision metric: RO <" }} />
        <CCP4i2TaskElement itemName="OUTPUT_PDB_FILE" {...props} qualifiers={{ guiLabel: "Output PDB file with per-atom metrics" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
