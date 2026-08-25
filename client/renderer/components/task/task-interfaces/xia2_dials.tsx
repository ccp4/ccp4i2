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
 * xia2 / DIALS.
 *
 * The def.xml is generated from xia2's master_phil and declares ~280
 * parameters. Enumerating them by hand is hopeless, so — as the Qt interface
 * did with `nestedAutoGenerate` — the parameter body is rendered from the
 * container tree and filtered by expert level: level 0 inline under "Basic
 * parameters", the rest behind the Advanced tab.
 *
 * Only "Locate datasets" is laid out by hand, because CXia2ImageSelectionList
 * and the directory picker need prose around them to explain that they are
 * alternatives.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);

  // Start at Basic and let the user opt upward. Advanced is 179 of xia2's
  // 281 parameters and All is all of them, which is a lot of widgets to build
  // before anyone has asked for them.
  const [expertLevel, setExpertLevel] = useState(0);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          {/* --- Locate datasets --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Locate datasets" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <Typography variant="body2" color="text.secondary">
                Choose one image from each dataset...
              </Typography>
              <CCP4i2TaskElement itemName="IMAGE_FILE" {...props} />
              <Typography variant="body2" color="text.secondary">
                ...Or let xia2 find datasets under a parent directory
              </Typography>
              <CCP4i2TaskElement itemName="IMAGE_DIRECTORY" {...props} />
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
              qualifiers={{ guiLabel: "xia2 parameters" }}
              containerHint="FolderLevel"
            />
          </ExpertLevelContext.Provider>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
