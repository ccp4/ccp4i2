import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const buccaneerPhsinType = useBoolToggle(useTaskItem, "BUCCANEER_PHSIN_TYPE");

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
        <CCP4i2TaskElement itemName="BUCCANEER_PHSIN_TYPE" {...props} qualifiers={{ guiLabel: "Build model with phases coming from" }} />
        {buccaneerPhsinType.value && (
          <CCP4i2TaskElement itemName="BUCCANEER_MR_MODE_XYZIN" {...props} />
        )}
        {buccaneerPhsinType.value && (
          <CCP4i2TaskElement itemName="BUCCANEER_MR_MODE" {...props} qualifiers={{ guiLabel: "This model will be used to place and name chains, and" }} />
        )}
        <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        {buccaneerPhsinType.value && (
          <CCP4i2TaskElement itemName="ABCD_EXP" {...props} />
        )}
        <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
        <CCP4i2TaskElement itemName="ASUIN" {...props} />
        <CCP4i2TaskElement itemName="XYZIN_MODE" {...props} />
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <CCP4i2TaskElement itemName="KNOWN_STRUCTURE" {...props} qualifiers={{ guiLabel: "Optionally keep part of the known structure fixed" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
