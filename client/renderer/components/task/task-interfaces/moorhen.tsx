import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

/** The recorded Moorhen session: the same inputs as Coot 1 (with a list of
 *  dictionaries and optional real-space maps). Run opens the session window;
 *  models saved there become this job's outputs. */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <Typography variant="body2" color="text.secondary" sx={{ px: 1 }}>
        Running this job opens Moorhen on the data below. Models you save from
        the session window become the job&apos;s outputs; finish the session
        from the window or the job menu.
      </Typography>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Coordinates" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="XYZIN_LIST" {...props} />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Electron density maps" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="FPHIIN_LIST" {...props} />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Difference density maps" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="DELFPHIIN_LIST" {...props} />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Anomalous difference maps" }}
        containerHint="FolderLevel"
        initiallyOpen={false}
      >
        <CCP4i2TaskElement itemName="DELFPHIINANOM_LIST" {...props} />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Real-space maps and masks" }}
        containerHint="FolderLevel"
        initiallyOpen={false}
      >
        <CCP4i2TaskElement itemName="MAPIN_LIST" {...props} />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Ligand geometry" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="DICT_LIST" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
