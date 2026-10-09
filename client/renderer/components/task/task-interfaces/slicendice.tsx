import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container, useTaskItem } = useJob(props.job.id);
  const { value: BFACTOR_TREATMENT } = useTaskItem("BFACTOR_TREATMENT");

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
        <CCP4i2TaskElement itemName="NO_MOLS" {...props} qualifiers={{ guiLabel: "Copies of the model in the AU (each piece is searched for this many times)" }} />
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <CCP4i2TaskElement itemName="BFACTOR_TREATMENT" {...props}
          qualifiers={{ guiLabel: "What the model's B-factor column holds" }} />
        {/* The threshold each treatment needs. Neither was here, so choosing a
            treatment left its cut-off at whatever the default happened to be. */}
        <CCP4i2TaskElement
          itemName="PLDDT_THRESHOLD"
          {...props}
          qualifiers={{ guiLabel: "pLDDT threshold" }}
          visibility={() => BFACTOR_TREATMENT === "plddt"}
        />
        <CCP4i2TaskElement
          itemName="RMS_THRESHOLD"
          {...props}
          qualifiers={{ guiLabel: "R.m.s. threshold" }}
          visibility={() => BFACTOR_TREATMENT === "rms"}
        />
      </CCP4i2ContainerElement>

      {/* Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="MIN_SPLITS" {...props}
          qualifiers={{ guiLabel: "Fewest pieces to split the model into" }} />
        <CCP4i2TaskElement itemName="MAX_SPLITS" {...props}
          qualifiers={{ guiLabel: "Most pieces to split the model into" }} />
        <CCP4i2TaskElement itemName="SGALTERNATIVE" {...props}
          qualifiers={{ guiLabel: "Space groups Phaser tests" }} />
        <CCP4i2TaskElement itemName="NPROC" {...props}
          qualifiers={{ guiLabel: "Splits to run at once (processors)" }} />
        <CCP4i2TaskElement itemName="NCYC" {...props}
          qualifiers={{ guiLabel: "Refmac cycles after each placement" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
