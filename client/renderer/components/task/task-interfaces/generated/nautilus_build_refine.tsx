import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: XYZIN_MODE } = useTaskItem("XYZIN_MODE");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Select experimental data
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="ABCD" {...props} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Enter the AU content containing the structure sequence(s)
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="ASUIN" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2TaskElement itemName="XYZIN_MODE" {...props} qualifiers={{ guiLabel: "Start from a partially built model", toolTip: "Your model should contain confidently-built metals, ligands and/or protein parts" }} />
          {(XYZIN_MODE === true) && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="XYZIN" {...props} />
            </CCP4i2ContainerElement>
          )}
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Options">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Pipeline control
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="REFMAC_CYCLES" {...props} qualifiers={{ guiLabel: "of refinement" }} />
          </CCP4i2ContainerElement>
          <CCP4i2TaskElement itemName="NAUTILUS_ANISOTROPY_CORRECTION" {...props} qualifiers={{ guiLabel: "Apply anisotropy correction" }} />
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Advanced Nautilus Options">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="NAUTILUS_RESOLUTION" {...props} qualifiers={{ guiLabel: "&#8491; resolution limit" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;