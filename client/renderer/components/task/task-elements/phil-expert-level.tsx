import { useEffect, useState } from "react";
import { Box, ToggleButton, ToggleButtonGroup, Typography } from "@mui/material";

import { Job } from "../../../types/models";
import { useJob } from "../../../utils";

/**
 * The expert-level control for a PHIL-driven task.
 *
 * PhilPluginScript injects PHIL_EXPERT_LEVEL into controlParameters, and it is
 * more than a display filter: the wrapper also uses it to decide which
 * parameters are written into working.phil. Keeping the value in the container
 * rather than in component state means the interface and the phil cannot
 * disagree, and the choice persists with the job.
 *
 * Level 0 is the default because a full PHIL tree is hundreds of widgets, and
 * building them before anyone has asked is what made these tasks slow.
 */
export const usePhilExpertLevel = (job: Job) => {
  const { useTaskItem } = useJob(job.id);
  const { item, update } = useTaskItem("controlParameters.PHIL_EXPERT_LEVEL");
  const [expertLevel, setExpertLevel] = useState(0);

  useEffect(() => {
    const value = item?._value;
    if (value !== undefined && value !== null) setExpertLevel(Number(value));
  }, [item?._value]);

  const changeExpertLevel = (value: number) => {
    setExpertLevel(value);
    update(value);
  };

  return { expertLevel, changeExpertLevel };
};

/** Keep PHIL_EXPERT_LEVEL out of the rendered parameter tree — it is a control. */
export const EXCLUDE_EXPERT_LEVEL = ["PHIL_EXPERT_LEVEL"];

interface PhilExpertLevelSelectorProps {
  expertLevel: number;
  onChange: (value: number) => void;
}

export const PhilExpertLevelSelector: React.FC<
  PhilExpertLevelSelectorProps
> = ({ expertLevel, onChange }) => (
  <Box sx={{ mx: 2, mb: 1, display: "flex", alignItems: "center", gap: 1 }}>
    <Typography variant="body2" sx={{ fontWeight: 600 }}>
      Expert level:
    </Typography>
    <ToggleButtonGroup
      value={expertLevel}
      exclusive
      onChange={(_, value) => value !== null && onChange(value)}
      size="small"
    >
      <ToggleButton value={0}>Basic</ToggleButton>
      <ToggleButton value={1}>Advanced</ToggleButton>
      <ToggleButton value={10}>All</ToggleButton>
    </ToggleButtonGroup>
  </Box>
);
