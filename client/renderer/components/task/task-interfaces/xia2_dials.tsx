import {
  LinearProgress,
  Paper,
  Typography,
} from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { ExpertLevelContext } from "../task-elements/expert-level-context";
import {
  EXCLUDE_EXPERT_LEVEL,
  PhilExpertLevelSelector,
  usePhilExpertLevel,
} from "../task-elements/phil-expert-level";
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
  const { container, useTaskItem } = useJob(props.job.id);

  const { expertLevel, changeExpertLevel } = usePhilExpertLevel(props.job);

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
              excludeItems={EXCLUDE_EXPERT_LEVEL}
            />
          </ExpertLevelContext.Provider>
        </CCP4i2Tab>

        <CCP4i2Tab key="controlParameters" label="Advanced parameters">
          <PhilExpertLevelSelector
            expertLevel={expertLevel}
            onChange={changeExpertLevel}
          />
          <ExpertLevelContext.Provider value={expertLevel}>
            <CCP4i2ContainerElement
              {...props}
              itemName="controlParameters"
              qualifiers={{ guiLabel: "xia2 parameters" }}
              containerHint="FolderLevel"
              excludeItems={EXCLUDE_EXPERT_LEVEL}
            />
          </ExpertLevelContext.Provider>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
