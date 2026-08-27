import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const perform = useBoolToggle(useTaskItem, "PERFORM");
  const { value: HIGH_PATH_VAR } = useTaskItem("HIGH_PATH_VAR");

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
        <CCP4i2TaskElement itemName="PERFORM" {...props} />
        {/* No F_SIGF field here: this task's def.xml declares none. The
            molrep wrappers share one plugin, and its F_SIGF branch is the
            one PERFORM != 'den' takes -- which is never, here. */}
        {perform.value && (
          <CCP4i2TaskElement itemName="XYZIN_FIX" {...props} />
        )}
        {!perform.value && (
          <CCP4i2TaskElement itemName="XYZIN" {...props} />
        )}
        {!perform.value && (
          <CCP4i2TaskElement itemName="ASUIN" {...props} />
        )}
        {!perform.value && (
          <CCP4i2TaskElement itemName="NMON" {...props} qualifiers={{ guiLabel: "The number of monomers to search for" }} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="F_PHI_MAP" {...props} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="XYZIN_FIX" {...props} />
        )}
      </CCP4i2ContainerElement>

      {/* Basic Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Basic Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="SEQ" {...props} />
        <CCP4i2TaskElement itemName="SURF" {...props} />
        <CCP4i2TaskElement itemName="NP" {...props} qualifiers={{ guiLabel: "Number of Rotation Function peaks" }} />
        <CCP4i2TaskElement itemName="NPT" {...props} qualifiers={{ guiLabel: "Number of Translation Function peaks" }} />
      </CCP4i2ContainerElement>

      {/* Advanced Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Advanced Options" }}
        containerHint="FolderLevel"
      >
        {perform.value && (
          <CCP4i2TaskElement itemName="PRF" {...props} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="SCORE" {...props} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="NMON_EXP" {...props} qualifiers={{ guiLabel: "Expected number of copies (for contrast calculation only)" }} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="ANISO" {...props} />
        )}
        <CCP4i2TaskElement itemName="HIGH_PATH_VAR" {...props} />
        {/* HIGH_PATH_VAR chooses how B-add is arrived at, and each route needs
            its own value. None of the three was here, so the choice had no
            effect a user could see. */}
        <CCP4i2TaskElement
          itemName="SIM"
          {...props}
          qualifiers={{
            guiLabel: "Identity between search and target sequences (0 to 1)",
          }}
          visibility={() => HIGH_PATH_VAR === "i"}
        />
        <CCP4i2TaskElement
          itemName="RESMAX"
          {...props}
          qualifiers={{ guiLabel: "High resolution limit (Å)" }}
          visibility={() => HIGH_PATH_VAR === "r"}
        />
        <CCP4i2TaskElement
          itemName="BADD"
          {...props}
          qualifiers={{ guiLabel: "B-add" }}
          visibility={() => HIGH_PATH_VAR === "b"}
        />
        <CCP4i2TaskElement itemName="LOW_PATH_VAR" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
