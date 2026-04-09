import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Input data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="XYZIN_LIST" {...props} />
        <CCP4i2TaskElement
          itemName="OVERRIDEID"
          {...props}
          qualifiers={{ guiLabel: "Sequence identity to place in header" }}
        />
        <CCP4i2TaskElement
          itemName="ALIGNIN"
          {...props}
          qualifiers={{ guiLabel: "Optional multiple sequence alignment" }}
        />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
