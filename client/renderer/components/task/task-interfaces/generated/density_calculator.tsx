import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: BLUR_MODE } = useTaskItem("BLUR_MODE");
  const { value: FORM_FACTOR } = useTaskItem("FORM_FACTOR");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Structure
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
            <CCP4i2TaskElement itemName="FORM_FACTOR" {...props} qualifiers={{ guiLabel: "Scattering form factor" }} />
            <CCP4i2TaskElement itemName="D_MIN" {...props} qualifiers={{ guiLabel: "Resolution limit / &#8491;" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Advanced Options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="RATE" {...props} qualifiers={{ guiLabel: "Oversampling rate" }} />
            <CCP4i2TaskElement itemName="BLUR_MODE" {...props} qualifiers={{ guiLabel: "Blurring mode" }} />
            {(true) && (
              <CCP4i2TaskElement itemName="BLUR" {...props} qualifiers={{ guiLabel: "Blurring B-factor / &#8491;2" }} />
            )}
            {(true) && (
              <CCP4i2TaskElement itemName="UNBLUR" {...props} qualifiers={{ guiLabel: "Unblur when calculating the reciprocal-space map coefficients" }} />
            )}
            <CCP4i2TaskElement itemName="CUTOFF" {...props} qualifiers={{ guiLabel: "Density cutoff" }} />
            {(true) && (
              <CCP4i2TaskElement itemName="MOTT_BETHE" {...props} qualifiers={{ guiLabel: "Approximate electron scattering factors using the Mott-Bethe formula" }} />
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;