import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: GLYTOUCAN } = useTaskItem("GLYTOUCAN");
  const { value: NEW_SUGAR } = useTaskItem("NEW_SUGAR");
  const { value: RING_TYPE } = useTaskItem("RING_TYPE");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Model and experimental data
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="RADIUSIN" {...props} qualifiers={{ guiLabel: "Angstroems" }} />
            <CCP4i2TaskElement itemName="NEW_SUGAR" {...props} qualifiers={{ guiLabel: "A sugar I want to validate is not yet part of the Chemical Component Dictionary" }} />
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              Define a new sugar
            </Typography>
            {(NEW_SUGAR === true) && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ initiallyOpen: true }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement itemName="RING_TYPE" {...props} qualifiers={{ guiLabel: "and type" }} />
                {(true) && (
                  <CCP4i2TaskElement itemName="CONFORMATION_PYRANOSE" {...props} qualifiers={{ guiLabel: "Expected minimal energy ring conformation:" }} />
                )}
                {(true) && (
                  <CCP4i2TaskElement itemName="CONFORMATION_FURANOSE" {...props} qualifiers={{ guiLabel: "Expected minimal energy ring conformation:" }} />
                )}
                {(true) && (
                  <CCP4i2TaskElement itemName="RING_C5" {...props} qualifiers={{ guiLabel: "Ring atoms:" }} />
                )}
                {(true) && (
                  <CCP4i2TaskElement itemName="RING_C4" {...props} qualifiers={{ guiLabel: "Ring atoms:" }} />
                )}
              </CCP4i2ContainerElement>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="folder_1" label="Settings">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Glycosylation analysis
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="VERTICAL" {...props} qualifiers={{ guiLabel: "orientation" }} />
            <CCP4i2TaskElement itemName="INVERT" {...props} qualifiers={{ guiLabel: "outlines" }} />
            <CCP4i2TaskElement itemName="EXPRESSION" {...props} qualifiers={{ guiLabel: "expression system" }} />
            <CCP4i2TaskElement itemName="GLYTOUCAN" {...props} qualifiers={{ guiLabel: "Validate Glycans with GlyTouCan and GlyConnect databases" }} />
            {(GLYTOUCAN === true) && (
              <>
                <CCP4i2TaskElement itemName="CLOSESTMATCH" {...props} qualifiers={{ guiLabel: "GlyConnect: Don't look for the closest match, if modelled glycan is not found in the database." }} />
                <CCP4i2TaskElement itemName="ALLPERMUTATIONS" {...props} qualifiers={{ guiLabel: "Conduct all possible permutation combinations(WARNING: Will take a long time to finish, should only be really used for O-Glycans)." }} />
              </>
            )}
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              Performance settings
            </Typography>
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="NUMTHREADS" {...props} qualifiers={{ guiLabel: "threads." }} />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;