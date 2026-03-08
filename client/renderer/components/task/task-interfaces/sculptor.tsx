import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: alignMode } = useTaskItem("ALIGNMENTORSEQUENCEIN");

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
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        {alignMode === "SEQUENCE" && (
          <CCP4i2TaskElement itemName="CHAINIDS" {...props} qualifiers={{ guiLabel: "Chain ID in template to manipulate" }} />
        )}
        <CCP4i2TaskElement itemName="ALIGNMENTORSEQUENCEIN" {...props} qualifiers={{ guiLabel: "Provide target sequence as :" }} />
        {alignMode === "ALIGNMENT" && (
          <CCP4i2TaskElement itemName="ALIGNIN" {...props} />
        )}
        {alignMode === "ALIGNMENT" && (
          <CCP4i2TaskElement itemName="TARGETINDEX" {...props} qualifiers={{ guiLabel: "Identifier of the target sequence in this alignemnt" }} />
        )}
        {alignMode === "SEQUENCE" && (
          <CCP4i2TaskElement itemName="SEQUENCEIN" {...props} />
        )}
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
