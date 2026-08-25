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
 * Phaser, driven from its PHIL tree (~270 generated parameters).
 *
 * As in the Qt interface, only the two file inputs are laid out by hand; the
 * parameter body is rendered from the container tree and filtered by expert
 * level rather than enumerated.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);
  const { expertLevel, changeExpertLevel } = usePhilExpertLevel(props.job);

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
              qualifiers={{ guiLabel: "Phaser parameters" }}
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
