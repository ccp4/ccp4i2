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
      {/* Protocol */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Protocol" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="MODE" {...props} qualifiers={{ guiLabel: "from an MTZ file" }} />
      </CCP4i2ContainerElement>

      {/* Input Data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input Data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="HKLIN1" {...props} />
        <CCP4i2TaskElement itemName="FSIGF" {...props} qualifiers={{ guiLabel: "Columns to be included/excluded" }} />
        <CCP4i2TaskElement itemName="HKLIN2" {...props} />
      </CCP4i2ContainerElement>

      {/* Output Data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Output Data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="HKLOUT" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
