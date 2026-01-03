import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: MAP_OUTPUT } = useTaskItem("MAP_OUTPUT");
  const { value: SCALE } = useTaskItem("SCALE");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Type of map to create
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="MAPTYPE" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            First dataset
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF1" {...props} />
            <CCP4i2TaskElement itemName="ABCD1" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Second dataset
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF2" {...props} />
            <CCP4i2TaskElement itemName="ABCD2" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Basic operations on map coefficients
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="B_VALUE" {...props} qualifiers={{ guiLabel: "Sharpen or blur the output map using an isotropic B-factor:", toolTip: "Use a negative B-factor to sharpen the map, or a positive one to blur it" }} />
            <CCP4i2TaskElement itemName="SCALE" {...props} qualifiers={{ guiLabel: "scale the observations" }} />
            {(SCALE === true) && (
              <CCP4i2TaskElement itemName="F1_TO_F2" {...props} />
            )}
            <CCP4i2TaskElement itemName="RESOLUTION" {...props} qualifiers={{ guiLabel: "Choose a different resolution limit" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Output options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="MAP_OUTPUT" {...props} qualifiers={{ guiLabel: "produce a map file too" }} />
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Grid sampling should not be changed unless you know what you are doing.
            </Typography>
            {(MAP_OUTPUT === true) && (
              <CCP4i2TaskElement itemName="INDEX_W" {...props} qualifiers={{ guiLabel: "w" }} />
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;