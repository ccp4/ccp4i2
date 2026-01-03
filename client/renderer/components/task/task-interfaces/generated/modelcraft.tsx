import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: BASIC } = useTaskItem("BASIC");
  const { value: USE_MODEL_PHASES } = useTaskItem("USE_MODEL_PHASES");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Reflection data
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
            <CCP4i2TaskElement itemName="USE_MODEL_PHASES" {...props} qualifiers={{ guiLabel: "Get initial phases from refining the starting model (uncheck to specify starting phases, e.g. from experimental phasing)" }} />
            {(USE_MODEL_PHASES === false) && (
              <CCP4i2TaskElement itemName="PHASES" {...props} />
            )}
            {(USE_MODEL_PHASES === false) && (
              <CCP4i2TaskElement itemName="UNBIASED" {...props} qualifiers={{ guiLabel: "Phases are unbiased and should be used as refinement restraints when the model is poor" }} />
            )}
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Asymmetric unit contents
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="ASUIN" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Starting model
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="BASIC" {...props} qualifiers={{ guiLabel: "Run a quicker basic pipeline" }} />
            <CCP4i2TaskElement itemName="CYCLES" {...props} qualifiers={{ guiLabel: "cycles" }} />
            <CCP4i2TaskElement itemName="STOP_CYCLES" {...props} qualifiers={{ guiLabel: "cycles" }} />
            <CCP4i2TaskElement itemName="SELENOMET" {...props} qualifiers={{ guiLabel: "Build selenomethionine (MSE) instead of methionine (MET)" }} />
            <CCP4i2TaskElement itemName="TWINNED" {...props} qualifiers={{ guiLabel: "Use twinned refinement" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Optional pipeline steps
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Uncheck to disable the following optional steps:
            </Typography>
            <CCP4i2TaskElement itemName="SHEETBEND" {...props} qualifiers={{ guiLabel: "Preliminary low-resolution refinement with Sheetbend" }} />
            {(BASIC === false) && (
              <CCP4i2TaskElement itemName="PRUNING" {...props} qualifiers={{ guiLabel: "Residue and chain pruning" }} />
            )}
            <CCP4i2TaskElement itemName="PARROT" {...props} qualifiers={{ guiLabel: "Classical density modification with Parrot" }} />
            {(BASIC === false) && (
              <CCP4i2TaskElement itemName="DUMMY_ATOMS" {...props} qualifiers={{ guiLabel: "Phase improvement through addition and refinement of dummy atoms" }} />
            )}
            {(BASIC === false) && (
              <CCP4i2TaskElement itemName="WATERS" {...props} qualifiers={{ guiLabel: "Addition of waters" }} />
            )}
            {(BASIC === false) && (
              <CCP4i2TaskElement itemName="SIDE_CHAIN_FIXING" {...props} qualifiers={{ guiLabel: "Final side-chain fixing" }} />
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;