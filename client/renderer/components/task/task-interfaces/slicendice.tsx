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
      {/* Input Data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input Data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
        <CCP4i2TaskElement itemName="ASUIN" {...props} />
        <CCP4i2TaskElement itemName="NO_MOLS" {...props} qualifiers={{ guiLabel: "The number of monomers to search for" }} />
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <CCP4i2TaskElement itemName="BFACTOR_TREATMENT" {...props} />
      </CCP4i2ContainerElement>

      {/* Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="MIN_SPLITS" {...props} />
        <CCP4i2TaskElement itemName="MAX_SPLITS" {...props} />
        <CCP4i2TaskElement itemName="NPROC" {...props} />
        <CCP4i2TaskElement itemName="NCYC" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
