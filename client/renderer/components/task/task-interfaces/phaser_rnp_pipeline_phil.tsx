import { LinearProgress, Paper } from "@mui/material";
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
 * Phaser rigid-body refinement pipeline over Phaser's own PHIL: a parent
 * model cut into rigid bodies by atom selections, the composition, a Free-R
 * set and the follow-on steps. The typed inputs are laid out by hand; the
 * composition source shows only the input it needs; the parameters are
 * rendered from the container tree, filtered by expert level.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container, useTaskItem } = useJob(props.job.id);
  const { expertLevel, changeExpertLevel } = usePhilExpertLevel(props.job);
  const { value: compBy } = useTaskItem("COMP_BY");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Reflections" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Rigid bodies" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN_PARENT" {...props} />
            <CCP4i2TaskElement itemName="SELECTIONS" {...props} />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Composition of the asymmetric unit" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="COMP_BY" {...props} />
            <CCP4i2TaskElement
              itemName="ASUFILE"
              {...props}
              visibility={() => compBy === "ASU"}
            />
            <CCP4i2TaskElement
              itemName="SEQUENCES"
              {...props}
              visibility={() => compBy === "SEQUENCES"}
            />
            <CCP4i2TaskElement
              itemName="ASU_PROTEIN_MW"
              {...props}
              visibility={() => compBy === "MW"}
            />
            <CCP4i2TaskElement
              itemName="ASU_NUCLEICACID_MW"
              {...props}
              visibility={() => compBy === "MW"}
            />
            <CCP4i2TaskElement
              itemName="SOLVENT_FRACTION"
              {...props}
              visibility={() => compBy === "SOLVENT"}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="afterwards" label="After Phaser">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Refinement of the solution" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="RUNSHEETBEND" {...props} />
            <CCP4i2TaskElement itemName="RUNREFMAC" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="controlParameters" label="Phaser parameters">
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
