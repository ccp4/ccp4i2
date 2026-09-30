import { useCallback } from "react";
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
  const { container, useTaskItem, fetchDigest } = useJob(props.job.id);
  const { expertLevel, changeExpertLevel } = usePhilExpertLevel(props.job);
  const { value: compBy } = useTaskItem("COMP_BY");
  const { value: partialBy } = useTaskItem("PARTIAL_BY");

  // The wavelength is in the reflection file: take it from there when the
  // data are chosen, as the classic phaser_EP interface did, rather than ask
  // for what the page already knows. Still editable, e.g. for a remote edge.
  const { item: fSigFItem } = useTaskItem("F_SIGF");
  const { forceUpdate: setWavelength } = useTaskItem("WAVELENGTH");
  const takeWavelengthFromData = useCallback(async () => {
    if (!fSigFItem?._objectPath) return;
    const digest = await fetchDigest(fSigFItem._objectPath);
    const wavelength = digest?.wavelengths?.at(-1);
    if (wavelength && wavelength > 0 && wavelength < 9) {
      await setWavelength(wavelength);
    }
  }, [fSigFItem?._objectPath, fetchDigest, setWavelength]);

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
            <CCP4i2TaskElement
              itemName="F_SIGF"
              {...props}
              onChange={takeWavelengthFromData}
            />
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
