import { Box, Typography } from "@mui/material";
import { Job } from "../types/models";
import { CCP4i2JobAvatar } from "./job-avatar";
import { useUiPreference } from "../lib/ui-preferences";

/**
 * What follows the cursor while a job is dragged out of the job list.
 *
 * It shows the job as a compact, row-like chip rather than a generic offset
 * icon, so what moves reads as the row itself — which is what the drag does.
 * The icon is the job list's own CCP4i2JobAvatar, so the chip carries the
 * task's icon at the row's size instead of a generic CCP4i2 diamond squeezed
 * into a 20px circle (which clipped it top and bottom), and it follows the
 * same "show job icons" preference as the rows: with icons off the chip is
 * text only (Paul Bond, issue #432).
 */
export const JobDragChip: React.FC<{ job: Job }> = ({ job }) => {
  const [showJobIcons] = useUiPreference("showJobIcons");

  return (
    <Box
      sx={{
        display: "inline-flex",
        alignItems: "center",
        gap: 1,
        px: 1.5,
        py: 0.5,
        maxWidth: 320,
        bgcolor: "background.paper",
        border: 1,
        borderColor: "divider",
        borderRadius: 1,
        boxShadow: 4,
        opacity: 0.95,
        pointerEvents: "none",
        cursor: "grabbing",
      }}
    >
      {showJobIcons ? <CCP4i2JobAvatar job={job} /> : null}
      <Typography
        variant="body2"
        sx={{
          fontWeight: 600,
          whiteSpace: "nowrap",
          overflow: "hidden",
          textOverflow: "ellipsis",
        }}
      >
        {job.number ? `${job.number}. ` : ""}
        {job.title || job.task_name || "Job"}
      </Typography>
    </Box>
  );
};
