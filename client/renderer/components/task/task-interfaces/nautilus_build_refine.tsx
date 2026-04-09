import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const xyzinMode = useBoolToggle(useTaskItem, "XYZIN_MODE");

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
        <Typography variant="subtitle2">Select experimental data</Typography>
        <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        <CCP4i2TaskElement itemName="ABCD" {...props} />
        <CCP4i2TaskElement itemName="FREERFLAG" {...props} />

        <Typography variant="subtitle2">Enter the AU content containing the structure sequence(s)</Typography>
        <CCP4i2TaskElement itemName="ASUIN" {...props} />

        <CCP4i2TaskElement itemName="XYZIN_MODE" {...props} qualifiers={{ guiLabel: "Start from a partially built model" }} />
        {xyzinMode.value && (
          <CCP4i2TaskElement itemName="XYZIN" {...props} />
        )}
      </CCP4i2ContainerElement>

      {/* Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <Typography variant="subtitle2">Pipeline control</Typography>
        <Typography variant="caption" color="text.secondary">Warning: default values work best for most situations.</Typography>
        <CCP4i2TaskElement itemName="ITERATIONS" {...props} qualifiers={{ guiLabel: "Iterate" }} />
        <CCP4i2TaskElement itemName="NAUTILUS_CYCLES" {...props} qualifiers={{ guiLabel: "cycles of build" }} />
        <CCP4i2TaskElement itemName="REFMAC_CYCLES" {...props} qualifiers={{ guiLabel: "cycles of refinement" }} />
        <CCP4i2TaskElement itemName="NAUTILUS_ANISOTROPY_CORRECTION" {...props} qualifiers={{ guiLabel: "Apply anisotropy correction" }} />
      </CCP4i2ContainerElement>

      {/* Advanced Nautilus Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Advanced Nautilus Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="NAUTILUS_RESOLUTION" {...props} qualifiers={{ guiLabel: "Leave out reflection data beyond a resolution limit (\u00C5)" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
