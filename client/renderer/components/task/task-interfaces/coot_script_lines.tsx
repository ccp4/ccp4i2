import { useCallback } from "react";
import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { COOT_SCRIPT_TEMPLATES } from "./coot_script_templates";

// Scripted Coot: choosing a starting point puts its script in the script
// box, to run as it is or to edit, as the Qt interface did. (Without this
// interface the menu set a parameter nothing read: every choice ran the
// default script.)
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem } = useJob(props.job.id);
  const { forceUpdate: setScript } = useTaskItem("SCRIPT");

  const handleStartpoint = useCallback(
    async (updatedItem: any) => {
      const choice = String(updatedItem?._value ?? updatedItem ?? "");
      const script = COOT_SCRIPT_TEMPLATES[choice];
      if (script !== undefined) await setScript(script);
    },
    [setScript]
  );

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement {...props} itemName=""
        qualifiers={{ guiLabel: "Models and maps" }} containerHint="FolderLevel">
        <CCP4i2TaskElement itemName="XYZIN" {...props}
          qualifiers={{ guiLabel: "Models (MolHandle_1, MolHandle_2, ...)" }} />
        <CCP4i2TaskElement itemName="FPHIIN" {...props}
          qualifiers={{ guiLabel: "Maps (MapHandle_1, ...)" }} />
        <CCP4i2TaskElement itemName="DELFPHIIN" {...props}
          qualifiers={{ guiLabel: "Difference maps (DifmapHandle_1, ...)" }} />
        <CCP4i2TaskElement itemName="DICT" {...props}
          qualifiers={{ guiLabel: "Ligand dictionary" }} />
      </CCP4i2ContainerElement>
      <CCP4i2ContainerElement {...props} itemName=""
        qualifiers={{ guiLabel: "Script" }} containerHint="FolderLevel">
        <CCP4i2TaskElement itemName="STARTPOINT" {...props}
          qualifiers={{ guiLabel: "Start from" }} onChange={handleStartpoint} />
        <Typography variant="body2" color="text.secondary">
          Choosing a starting point replaces the script below with it; edit it before running
          if you wish. The script is Python run in Coot; it must write its result into dropDir.
        </Typography>
        <CCP4i2TaskElement itemName="SCRIPT" {...props}
          qualifiers={{ guiLabel: "Coot script", guiMode: "multiLine" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
