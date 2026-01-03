import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: PROVIDEDICT } = useTaskItem("PROVIDEDICT");
  const { value: PROVIDESEQUENCES } = useTaskItem("PROVIDESEQUENCES");
  const { value: PROVIDETLS } = useTaskItem("PROVIDETLS");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <CCP4i2TaskElement itemName="USINGIORF" {...props} qualifiers={{ guiLabel: "Refinement in final step used reflection" }} />
          <CCP4i2TaskElement itemName="PROVIDESEQUENCES" {...props} qualifiers={{ guiLabel: "Provide sequences of crystallized species" }} />
          {(PROVIDESEQUENCES === true) && (
            <CCP4i2TaskElement itemName="ASUIN" {...props} />
          )}
          <CCP4i2TaskElement itemName="PROVIDETLS" {...props} qualifiers={{ guiLabel: "Provide Refined TLS parameters" }} />
          {(PROVIDETLS === true) && (
            <CCP4i2TaskElement itemName="TLSIN" {...props} />
          )}
          <CCP4i2TaskElement itemName="PROVIDEDICT" {...props} qualifiers={{ guiLabel: "Provide DICT file for ligand in structure" }} />
          {(PROVIDEDICT === true) && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="DICT" {...props} />
            </CCP4i2ContainerElement>
          )}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "NEW: Allow ANISO use in zero cycle refinement", initiallyOpen: true }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="B_REFINEMENT_MODE" {...props} qualifiers={{ guiLabel: "Temperature factors" }} />
          </CCP4i2ContainerElement>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Choose a directory in which to put "Reflections.cif" and "Coordinates.cif"
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              which can be uploaded to the deposition service
            </Typography>
            <CCP4i2TaskElement itemName="OUTPUT_DIRECTORY" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;