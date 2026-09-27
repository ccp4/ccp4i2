import { Avatar } from "@mui/material";
import { Job } from "../types/models";
import { forwardRef, useMemo } from "react";
import { useDraggable } from "@dnd-kit/core";
import { useTheme } from "../theme/theme-provider";
import { getStatusColour, isWorking } from "../lib/job-status-colour";

interface CCP4i2JobAvatarProps {
  job: Job;
}
export const CCP4i2JobAvatar = forwardRef<HTMLDivElement, CCP4i2JobAvatarProps>(
  ({ job, ...props }, ref) => {
    const { customColors } = useTheme();
    const bgColor = useMemo(() => getStatusColour(job?.status ?? 0), [job]);

    // Animation for a job that is doing something -- including one whose
    // program runs on a Batch node, which is no less running for being
    // somewhere else. Gated on RUNNING alone, a dispatched job sat looking
    // inert for the hours its run took.
    const runningAnimation =
      isWorking(job?.status ?? 0)
        ? {
            animation: "avatarPulse 2s infinite",
            "@keyframes avatarPulse": {
              "0%": { boxShadow: "0 0 0 0 rgba(25,118,210,0.7)" },
              "70%": { boxShadow: "0 0 0 10px rgba(25,118,210,0)" },
              "100%": { boxShadow: "0 0 0 0 rgba(25,118,210,0)" },
            },
          }
        : {};

    return (
      <Avatar
        {...props}
        ref={ref}
        sx={{
          width: "2rem",
          height: "2rem",
          backgroundColor: bgColor,
          border: `2px dashed ${customColors.ui.lightBlue}`,
          padding: "4px",
          cursor: "grab",
          transition: "box-shadow 0.2s ease",
          "&:hover": {
            boxShadow: "0 0 0 3px rgba(25, 118, 210, 0.5)",
          },
          ...runningAnimation,
        }}
        src={`/svgicons/${job.task_name}.svg`}
        alt={job.task_name}
      >
        {job.task_name?.[0]?.toUpperCase()}
      </Avatar>
    );
  }
);
