import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: compBy } = useTaskItem("COMP_BY");
  const { value: sgaltSelect } = useTaskItem("SGALT_SELECT");

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
        <CCP4i2TaskElement itemName="COMP_BY" {...props} qualifiers={{ guiLabel: "For estimating asymmetric unit contents:" }} />
        {compBy === "ASU" && (
          <CCP4i2TaskElement itemName="ASU_COMPONENTS" {...props} qualifiers={{ guiLabel: "Contents of asymmetric unit" }} />
        )}
        {compBy === "MW" && (
          <>
            <CCP4i2TaskElement itemName="ASU_PROTEIN_MW" {...props} />
            <CCP4i2TaskElement itemName="ASU_NUCLEICACID_MW" {...props} />
          </>
        )}
        <CCP4i2TaskElement itemName="ENSEMBLES" {...props} qualifiers={{ guiLabel: "Search model(s) - click \"Show list\" if more than one copy or more than one search model" }} />
        <CCP4i2TaskElement itemName="SEARCHMODE" {...props} qualifiers={{ guiLabel: "Use the ensembles above as" }} />
      </CCP4i2ContainerElement>

      {/* Important Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Important Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="SGALT_SELECT" {...props} />
        {sgaltSelect === "LIST" && (
          <CCP4i2TaskElement itemName="SGALT_TEST" {...props} />
        )}
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
