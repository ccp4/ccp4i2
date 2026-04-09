import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: PERFORM } = useTaskItem("PERFORM");
  const { value: SIMBAD_SEARCH_LEVEL } = useTaskItem("SIMBAD_SEARCH_LEVEL");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          {(PERFORM === "den") && (
            <>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                Input reflections
              </Typography>
              <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            </>
          )}
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Search level:
          </Typography>
          <CCP4i2TaskElement itemName="SIMBAD_SEARCH_LEVEL" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Organism:
          </Typography>
          <CCP4i2TaskElement itemName="SIMBAD_ORGANISM" {...props} />
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              x
            </Typography>
            <CCP4i2TaskElement itemName="SIMBAD_NPROC" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="inputData" label="Advanced Options">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            label
          </Typography>
          <CCP4i2TaskElement itemName="SIMBAD_ROT_PROGRAM" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            label
          </Typography>
          <CCP4i2TaskElement itemName="SIMBAD_MR_PROGRAM" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            label
          </Typography>
          <CCP4i2TaskElement itemName="SIMBAD_NMOL" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            label
          </Typography>
          <CCP4i2TaskElement itemName="SIMBAD_PROCESS_ALL" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            label
          </Typography>
          <CCP4i2TaskElement itemName="SIMBAD_SGALTERNATIVE" {...props} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;