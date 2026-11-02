import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";

/**
 * Polish geometry for deposition.
 *
 * Two tabs: the files (model, maps, optional dictionary) and the knobs. The
 * knobs are grouped by stage so a user sees which switch disables which part
 * of the polish; the defaults are the ones that worked.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model and maps" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
            <CCP4i2TaskElement itemName="FPHIOUT" {...props} />
            <CCP4i2TaskElement itemName="DIFFPHIOUT" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Optional" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="DICT" {...props} />
            <CCP4i2TaskElement itemName="CHAIN_SELECTION" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="options" label="Options">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Stages" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="DO_ROTAMERS" {...props} />
            <CCP4i2TaskElement itemName="DO_RAMA" {...props} />
            <CCP4i2TaskElement itemName="DO_PEPTIDE_FLIPS" {...props} />
            <CCP4i2TaskElement itemName="DO_CLASHES" {...props} />
            <CCP4i2TaskElement itemName="USE_NCS" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Acceptance" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="DENSITY_GUARD" {...props} />
            <CCP4i2TaskElement itemName="MIN_DENSITY" {...props} />
            <CCP4i2TaskElement itemName="REFINE_CYCLES" {...props} />
            <CCP4i2TaskElement itemName="MIN_CLASH_OVERLAP" {...props} />
            <CCP4i2TaskElement itemName="MAX_CLASH_ROTAMERS" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
