import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { InlineField } from "../../task-elements/inline-field";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Main inputs */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Main inputs" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="XYZIN"
          {...props}
          qualifiers={{ toolTip: "Input structure" }}
        />
        <CCP4i2TaskElement
          itemName="FPHIIN"
          {...props}
          qualifiers={{ toolTip: "Input 2mFo-DFc coefficients" }}
        />
      </CCP4i2ContainerElement>

      {/* Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <InlineField label="Local radius">
          <CCP4i2TaskElement
            itemName="LOCAL_RADIUS"
            {...props}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
        <InlineField label="GM alpha">
          <CCP4i2TaskElement
            itemName="GM_ALPHA"
            {...props}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
        <InlineField label="Blur B-factor">
          <CCP4i2TaskElement
            itemName="BLUR_B_FACTOR"
            {...props}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
