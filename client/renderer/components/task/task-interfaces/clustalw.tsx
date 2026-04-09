import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { useApi } from "../../../api";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const api = useApi();
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: SEQUENCELISTORALIGNMENT } = useTaskItem(
    "SEQUENCELISTORALIGNMENT"
  );

  if (!container) return <LinearProgress />;

  const showSequenceList = useMemo(
    () => SEQUENCELISTORALIGNMENT === "SEQUENCELIST",
    [SEQUENCELISTORALIGNMENT]
  );
  const showAlignment = useMemo(
    () => SEQUENCELISTORALIGNMENT === "ALIGNMENT",
    [SEQUENCELISTORALIGNMENT]
  );
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="Input data" label="Input data">
          <CCP4i2TaskElement
            itemName="SEQUENCELISTORALIGNMENT"
            {...props}
            qualifiers={{ guiLabel: "Input type" }}
          />
          {showSequenceList && (
            <CCP4i2TaskElement
              itemName="SEQIN"
              {...props}
              qualifiers={{ guiLabel: "Sequences" }}
            />
          )}
          {showAlignment && (
            <CCP4i2TaskElement
              itemName="ALIGNMENTIN"
              {...props}
              qualifiers={{ guiLabel: "Alignment file" }}
              visibility={showAlignment}
            />
          )}
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
