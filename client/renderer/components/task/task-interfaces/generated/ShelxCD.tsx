import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: MODE } = useTaskItem("MODE");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <CCP4i2TaskElement itemName="The Shelx programs have not been found. They are not part of CCP4 but you can get them from
http://shelx.uni-ac.gwdg.de/SHELX/download.php
If you already have them make sure they are on the search path
OR specify where they are in the Preferences window - under Other Software." {...props} />
          <CCP4i2TaskElement itemName="MODE" {...props} qualifiers={{ guiLabel: "Experiment type", toolTip: "What sort of data  are available for finding heavy atoms" }} />
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="SFAC" {...props} qualifiers={{ guiLabel: "Atom type to seek" }} />
            <CCP4i2TaskElement itemName="FIND" {...props} qualifiers={{ guiLabel: "Number to find" }} />
            <CCP4i2TaskElement itemName="NTRY" {...props} qualifiers={{ guiLabel: "Number of trys" }} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Anomalous (SAD) dataset
            </Typography>
            {(MODE === "SAD") && (
              <CCP4i2TaskElement itemName="SAD" {...props} />
            )}
            {(MODE === "MAD") && (
              <>
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  High energy remote dataset
                </Typography>
                <CCP4i2TaskElement itemName="HREM" {...props} />
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Low energy remote dataset
                </Typography>
                <CCP4i2TaskElement itemName="LREM" {...props} />
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Inflection point dataset
                </Typography>
                <CCP4i2TaskElement itemName="INFL" {...props} />
                <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                  Peak dataset
                </Typography>
                <CCP4i2TaskElement itemName="PEAK" {...props} />
              </>
            )}
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Derivative dataset
            </Typography>
            {(MODE === "SIR") && (
              <CCP4i2TaskElement itemName="SIR" {...props} />
            )}
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Anomalous derivative dataset
            </Typography>
            {(MODE === "SIRAS") && (
              <CCP4i2TaskElement itemName="SIRA" {...props} />
            )}
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Native dataset
            </Typography>
            <CCP4i2TaskElement itemName="NAT" {...props} />
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Radiation damaged dataset
            </Typography>
            {(MODE === "RIP") && (
              <CCP4i2TaskElement itemName="RIP" {...props} />
            )}
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Radiation damaged anomalous dataset
            </Typography>
            {(MODE === "RIPAS") && (
              <CCP4i2TaskElement itemName="RIPA" {...props} />
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Keywords">
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;