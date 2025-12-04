import { Avatar } from "@mui/material";
import { Job } from "../types/models";
import { forwardRef, useMemo } from "react";
import { useDraggable } from "@dnd-kit/core";
import { useTheme } from "../theme/theme-provider";

interface CCP4i2JobAvatarProps {
  job: Job;
}
export const CCP4i2JobAvatar = forwardRef<HTMLDivElement, CCP4i2JobAvatarProps>(
  ({ job, ...props }, ref) => {
    const { customColors } = useTheme();
    const bgColor = useMemo(() => {
      switch (job?.status) {
        case 0:
          return "#AAA";
        case 1:
          return "#FFF";
        case 2:
          return "#FFA";
        case 3:
          return "#AAF";
        case 4:
          return "#FDA";
        case 5:
          return "#FAA";
        case 6:
          return "#AFA";
        default:
          return "#AAA";
      }
    }, [job]);

    // Animation for running status
    const runningAnimation =
      job?.status === 3
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
        src={`/api/proxy/djangostatic/svgicons/${job.task_name}.svg`}
        alt={job.task_name}
      />
    );
  }
);

// Add global keyframes for the animation
// Place this in your main CSS or a global style file:
/*
@keyframes avatarPulse {
  0% {
    box-shadow: 0 0 0 0 rgba(25,118,210,0.7);
  }
  70% {
    box-shadow: 0 0 0 10px rgba(25,118,210,0);
  }
  100% {
    box-shadow: 0 0 0 0 rgba(25,118,210,0);
  }
}
*/
