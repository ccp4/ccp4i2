import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Input reflections
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF" {...props} qualifiers={{ toolTip: "Input reflections" }} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Basic Options">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Parameters to refine
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="REFINE_COORD" {...props} qualifiers={{ guiLabel: "coordinates" }} />
            <CCP4i2TaskElement itemName="REFINE_U_ISO" {...props} qualifiers={{ guiLabel: "isotropic B factors" }} />
            <CCP4i2TaskElement itemName="REFINE_U_ANISO" {...props} qualifiers={{ guiLabel: "anisotropic B factors" }} />
          </CCP4i2ContainerElement>
          <CCP4i2TaskElement itemName="CYCLES" {...props} qualifiers={{ guiLabel: "Number of cycles" }} />
          <CCP4i2TaskElement itemName="RESOLUTION" {...props} qualifiers={{ guiLabel: "Angstroms", toolTip: "Enter a single number or a comma separated list" }} />
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Advanced Options">
          <CCP4i2TaskElement itemName="RADIUS_SCALE" {...props} qualifiers={{ guiLabel: "times resolution" }} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Refine-regularize macro cycles
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="PSEUDO_REGULARIZE" {...props} qualifiers={{ guiLabel: "Perform pseudo-regularization" }} />
            <CCP4i2TaskElement itemName="REFINE_REGULARIZE_CYCLES" {...props} qualifiers={{ guiLabel: "Number of refine-regularize macro-cycles" }} />
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              Additional parameters to refine before regularization
            </Typography>
            <CCP4i2TaskElement itemName="POSTREFINE_COORD" {...props} qualifiers={{ guiLabel: "coordinates" }} />
            <CCP4i2TaskElement itemName="POSTREFINE_U_ISO" {...props} qualifiers={{ guiLabel: "isotropic B factors" }} />
            <CCP4i2TaskElement itemName="POSTREFINE_U_ANISO" {...props} qualifiers={{ guiLabel: "anisotropic B factors" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;