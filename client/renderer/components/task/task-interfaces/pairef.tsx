import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: shType } = useTaskItem("SH_TYPE");
  const usePreref = useBoolToggle(useTaskItem, "USE_PREREF");
  const useShake = useBoolToggle(useTaskItem, "USE_SHAKE");
  const autoWgt = useBoolToggle(useTaskItem, "AUTO_WGT");
  const fixedTls = useBoolToggle(useTaskItem, "FIXED_TLS");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input data" }}
        containerHint="FolderLevel"
      >
        {/* The reflections must extend beyond the resolution the model
            was refined at: the shells to test are the ones past it. */}
        <CCP4i2TaskElement itemName="F_SIGF" {...props}
          qualifiers={{ guiLabel: "Reflections, to beyond the resolution the model was refined at" }} />
        <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <CCP4i2TaskElement itemName="DICT" {...props}
          qualifiers={{ guiLabel: "Restraint dictionary (the one the model was refined with)" }} />
        <CCP4i2TaskElement itemName="UNMERGED" {...props}
          qualifiers={{ guiLabel: "Unmerged reflections (optional: adds data statistics per shell)" }} />
        <CCP4i2TaskElement itemName="REFMAC_KEYWORD_FILE" {...props} />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Resolution shells" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="SH_TYPE" {...props} qualifiers={{ guiLabel: "Shells to add" }} />
        <CCP4i2TaskElement itemName="INIRES" {...props}
          qualifiers={{ guiLabel: "Starting resolution (Å; 0 reads it from the model)" }} />
        {shType === "semi" && (
          <>
            <CCP4i2TaskElement itemName="NSHELL" {...props} qualifiers={{ guiLabel: "Number of shells" }} />
            <CCP4i2TaskElement itemName="WSHELL" {...props} qualifiers={{ guiLabel: "Width of each shell (Å)" }} />
          </>
        )}
        {shType === "manual" && (
          <CCP4i2TaskElement itemName="MANSHELL" {...props}
            qualifiers={{ guiLabel: "Shell limits (Å), highest first, comma-separated: 2.1,2.0,1.9" }} />
        )}
        <CCP4i2TaskElement itemName="COMPLETE" {...props}
          qualifiers={{ guiLabel: "Complete cross-validation (repeat for every free set; slow)" }} />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Refinement" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="NCYCLES" {...props}
          qualifiers={{ guiLabel: "Refinement cycles at each resolution" }} />
        <CCP4i2TaskElement itemName="AUTO_WGT" {...props} qualifiers={{ guiLabel: "Automatic weighting" }} />
        {!autoWgt.value && (
          <CCP4i2TaskElement itemName="WGT_TRM" {...props} qualifiers={{ guiLabel: "Weight" }} />
        )}
        <CCP4i2TaskElement itemName="USE_PREREF" {...props}
          qualifiers={{ guiLabel: "Refine the model at the starting resolution first" }} />
        {usePreref.value && (
          <>
            <CCP4i2TaskElement itemName="NPRECYCLES" {...props}
              qualifiers={{ guiLabel: "Pre-refinement cycles" }} />
            <CCP4i2TaskElement itemName="RESETBFAC" {...props}
              qualifiers={{ guiLabel: "Reset B-factors to their mean before pre-refinement" }} />
            <CCP4i2TaskElement itemName="USE_SHAKE" {...props}
              qualifiers={{ guiLabel: "Randomise coordinates before pre-refinement" }} />
            {useShake.value && (
              <CCP4i2TaskElement itemName="SHAKE" {...props}
                qualifiers={{ guiLabel: "Mean coordinate shift (Å)" }} />
            )}
          </>
        )}
        <CCP4i2TaskElement itemName="FIXED_TLS" {...props}
          qualifiers={{ guiLabel: "Refine TLS" }} />
        {fixedTls.value && (
          <>
            <CCP4i2TaskElement itemName="TLSIN" {...props} />
            <CCP4i2TaskElement itemName="TLSCYC" {...props}
              qualifiers={{ guiLabel: "TLS refinement cycles" }} />
          </>
        )}
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
