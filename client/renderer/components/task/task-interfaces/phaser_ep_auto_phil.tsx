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

/** Phaser EP_AUTO over its own PHIL; see phaser_mr_auto_phil.tsx. */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container, useTaskItem } = useJob(props.job.id);
  const { expertLevel, changeExpertLevel } = usePhilExpertLevel(props.job);
  const { value: compBy } = useTaskItem("COMP_BY");
  const { value: partialBy } = useTaskItem("PARTIAL_BY");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Anomalous data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="WAVELENGTH" {...props} />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Substructure" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN_HA" {...props} />
            <CCP4i2TaskElement itemName="ELEMENTS" {...props} />
            <CCP4i2TaskElement itemName="LLGC_CYCLES" {...props} />
            <CCP4i2TaskElement itemName="PURE_ANOMALOUS" {...props} />
            <CCP4i2TaskElement itemName="PARTIAL_BY" {...props} />
            <CCP4i2TaskElement
              itemName="XYZIN_PARTIAL"
              {...props}
              visibility={() => partialBy === "MODEL"}
            />
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
              itemName="SOLVENT_FRACTION"
              {...props}
              visibility={() => compBy === "SOLVENT"}
            />
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
