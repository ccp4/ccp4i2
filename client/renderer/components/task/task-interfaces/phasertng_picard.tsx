import { useState, useEffect } from "react";
import {
  Box,
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

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { item: expertLevelItem, update: updateExpertLevel } = useTaskItem(
    "controlParameters.PHIL_EXPERT_LEVEL"
  );

  const [expertLevel, setExpertLevel] = useState(0);

  // Sync from server/container value
  useEffect(() => {
    if (
      expertLevelItem?._value !== undefined &&
      expertLevelItem._value !== null
    ) {
      setExpertLevel(Number(expertLevelItem._value));
    }
  }, [expertLevelItem?._value]);

  const handleExpertLevelChange = (
    _: React.MouseEvent<HTMLElement>,
    newValue: number | null
  ) => {
    if (newValue !== null) {
      setExpertLevel(newValue);
      updateExpertLevel(newValue);
    }
  };

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Reflection data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="HKLIN" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Search models" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
            <CCP4i2TaskElement itemName="ASUIN" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Refine and continue (Riker)",
              initiallyOpen: false,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="FIXED" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Additional input", initiallyOpen: false }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="DICT" {...props} />
            <CCP4i2TaskElement itemName="DAGIN" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="keywords" label="Options">
          <Box
            sx={{
              mx: 2,
              mb: 1,
              display: "flex",
              alignItems: "center",
              gap: 1,
            }}
          >
            <Typography variant="body2" sx={{ fontWeight: 600 }}>
              Expert level:
            </Typography>
            <ToggleButtonGroup
              value={expertLevel}
              exclusive
              onChange={handleExpertLevelChange}
              size="small"
            >
              <ToggleButton value={0}>Basic</ToggleButton>
              <ToggleButton value={1}>Standard</ToggleButton>
              <ToggleButton value={3}>Expert</ToggleButton>
              <ToggleButton value={10}>All</ToggleButton>
            </ToggleButtonGroup>
          </Box>
          <ExpertLevelContext.Provider value={expertLevel}>
            <CCP4i2ContainerElement
              {...props}
              itemName="controlParameters"
              qualifiers={{ guiLabel: "PhaserTNG Parameters" }}
              excludeItems={["PHIL_EXPERT_LEVEL"]}
            />
          </ExpertLevelContext.Provider>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
