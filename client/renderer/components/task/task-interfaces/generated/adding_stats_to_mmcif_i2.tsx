import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* ASU content */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Asymmetric unit content (i.e. sequences)" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="ASUCONTENT"
          {...props}
          qualifiers={{ toolTip: "ASU Content" }}
        />
      </CCP4i2ContainerElement>

      {/* Coordinates */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Coordinates" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="XYZIN"
          {...props}
          qualifiers={{ toolTip: "Input model" }}
        />
      </CCP4i2ContainerElement>

      {/* Task control */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Task control" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="SENDTOVALIDATIONSERVER"
          {...props}
          qualifiers={{
            guiLabel: "Use validation server",
            toolTip: "Send to validation server - requires internet access",
          }}
        />
        <CCP4i2TaskElement
          itemName="USEAIMLESSXML"
          {...props}
          qualifiers={{
            guiLabel: "Add statistics from Aimless",
            toolTip:
              "If you didn't run data reduction in CCP4i2 you won't have this",
          }}
        />
        <CCP4i2TaskElement
          itemName="INCLUDEUNMERGED"
          {...props}
          qualifiers={{
            guiLabel: "Include unmerged from scaling job",
            toolTip:
              "If you didn't run data reduction in CCP4i2 you won't have this",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Related files */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Related files" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="USEANOMALOUS"
          {...props}
          qualifiers={{ guiLabel: "Anomalous used in refmac job" }}
        />
        <CCP4i2TaskElement
          itemName="USE_TWIN"
          {...props}
          qualifiers={{
            guiLabel: "Twinning (i.e. intensities) used in refmac job",
          }}
        />
        <CCP4i2TaskElement
          itemName="F_SIGF"
          {...props}
          qualifiers={{ toolTip: "Input reflections" }}
        />
        <CCP4i2TaskElement itemName="SCALEDUNMERGED" {...props} />
        <CCP4i2TaskElement itemName="AIMLESSXML" {...props} />
        <CCP4i2TaskElement itemName="REFMACINPUTPARAMSXML" {...props} />
        <CCP4i2TaskElement
          itemName="FREERFLAG"
          {...props}
          qualifiers={{ toolTip: "FreeR flag" }}
        />
        <CCP4i2TaskElement
          itemName="TLSIN"
          {...props}
          qualifiers={{ toolTip: "Input TLS" }}
        />
        <CCP4i2TaskElement
          itemName="DICT_LIST"
          {...props}
          qualifiers={{ toolTip: "Dictionary list" }}
        />
        <CCP4i2TaskElement itemName="FPHIOUT" {...props} />
        <CCP4i2TaskElement itemName="DIFFPHIOUT" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
