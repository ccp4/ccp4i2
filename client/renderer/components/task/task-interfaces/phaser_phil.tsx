import { useState } from "react";
import {
  Box,
  LinearProgress,
  Paper,
  ToggleButton,
  ToggleButtonGroup,
  Typography,
} from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { ExpertLevelContext } from "../task-elements/expert-level-context";
import { useJob } from "../../../utils";

/**
 * Phaser, driven from its PHIL tree (~270 generated parameters).
 *
 * As in the Qt interface, only the two file inputs are laid out by hand; the
 * parameter body is rendered from the container tree and filtered by expert
 * level rather than enumerated.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);
  const [expertLevel, setExpertLevel] = useState(0);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          {/* --- Input coordinates --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input coordinates" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="XYZIN" {...props} />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* --- Input reflections --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input reflections" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* --- Basic parameters: the whole PHIL tree at expert level 0 --- */}
          <ExpertLevelContext.Provider value={0}>
            <CCP4i2ContainerElement
              {...props}
              itemName="controlParameters"
              qualifiers={{ guiLabel: "Basic parameters" }}
              containerHint="FolderLevel"
            />
          </ExpertLevelContext.Provider>
        </CCP4i2Tab>

        <CCP4i2Tab key="controlParameters" label="Advanced parameters">
          <Box
            sx={{ mx: 2, mb: 1, display: "flex", alignItems: "center", gap: 1 }}
          >
            <Typography variant="body2" sx={{ fontWeight: 600 }}>
              Expert level:
            </Typography>
            <ToggleButtonGroup
              value={expertLevel}
              exclusive
              onChange={(_, value) => value !== null && setExpertLevel(value)}
              size="small"
            >
              <ToggleButton value={0}>Basic</ToggleButton>
              <ToggleButton value={1}>Advanced</ToggleButton>
              <ToggleButton value={10}>All</ToggleButton>
            </ToggleButtonGroup>
          </Box>
          <ExpertLevelContext.Provider value={expertLevel}>
            <CCP4i2ContainerElement
              {...props}
              itemName="controlParameters"
              qualifiers={{ guiLabel: "Phaser parameters" }}
              containerHint="FolderLevel"
            />
          </ExpertLevelContext.Provider>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
