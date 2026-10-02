import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  // A menu, not a yes/no: each file field belongs to one of its choices.
  // (Read through useBoolToggle, "reference" was never true and no file
  // field ever showed.)
  const { value: symmetrySource } = useTaskItem("SYMMETRY_SOURCE");

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
        <CCP4i2TaskElement itemName="HKLIN" {...props} />
        <CCP4i2TaskElement itemName="HKLIN1" {...props} />
        <CCP4i2TaskElement itemName="HKLIN2" {...props} />
        <CCP4i2TaskElement itemName="N_BINS" {...props} qualifiers={{ guiLabel: "Number of resolution bins" }} />
        <CCP4i2TaskElement itemName="SYMMETRY_SOURCE" {...props} qualifiers={{ guiLabel: "Load symmetry from" }} />
        {symmetrySource === "reference" && (
          <CCP4i2TaskElement itemName="REFERENCEFILE" {...props} />
        )}
        {symmetrySource === "cellfile" && (
          <CCP4i2TaskElement itemName="CELLFILE" {...props} />
        )}
        {symmetrySource === "streamfile" && (
          <CCP4i2TaskElement itemName="STREAMFILE" {...props} />
        )}
        <CCP4i2TaskElement itemName="SPACEGROUP" {...props} qualifiers={{ guiLabel: "Space group" }} />
        <CCP4i2TaskElement itemName="CELL" {...props} qualifiers={{ guiLabel: "Unit cell" }} />
        <CCP4i2TaskElement itemName="WAVELENGTH" {...props} qualifiers={{ guiLabel: "Wavelength (A)" }} />
        <CCP4i2TaskElement itemName="D_MAX" {...props} qualifiers={{ guiLabel: "Low resolution cutoff" }} />
        <CCP4i2TaskElement itemName="D_MIN" {...props} qualifiers={{ guiLabel: "High resolution cutoff" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
