import { Avatar, Button, Toolbar, Tooltip, Typography } from "@mui/material";
import { useApi } from "../../api";
import { useCCP4i2Window } from "../../app-context";
import { useJob, useProject } from "../../utils";
import { Job } from "../../types/models";
import {
  useCallback,
  useEffect,
  useLayoutEffect,
  useRef,
  useState,
} from "react";
import { useTheme } from "../../theme/theme-provider";
import { useRouter } from "next/navigation";
import { usePopcorn } from "../../providers/popcorn-provider";

// Label tiers, widest first. The bar is pinned so it cannot scroll vertically
// out of reach; when it is too narrow for the full list we shed label width
// rather than hiding suggestions, so every one stays a single click away.
const TIER_FULL = 0; // descriptive title, e.g. "Model building - COOT 0.9"
const TIER_SHORT = 1; // compact title, e.g. "Coot 0.9"
const TIER_ICON = 2; // icon only, identified by tooltip
const NARROWEST_TIER = TIER_ICON;

// Even with room to spare, a pathologically long title (ModelCraft's runs to 51
// characters) would drag every other label down a tier, so clamp it on its own.
const MAX_LABEL_CHARS = 30;

const useIsomorphicLayoutEffect =
  typeof window !== "undefined" ? useLayoutEffect : useEffect;

export const CCP4i2WhatNext = () => {
  const api = useApi();
  const { jobId } = useCCP4i2Window();
  const { job } = useJob(jobId);
  const { mutateJobs } = useProject(job?.project);
  const { customColors } = useTheme();
  const router = useRouter();
  const { setMessage } = usePopcorn();
  const { data: what_next } = api.get_endpoint<any>({
    type: "jobs",
    id: job?.id,
    endpoint: "what_next",
  });

  const handleTaskSelect = useCallback(
    async (task_name: string) => {
      if (!job) return;
      try {
        const created_job_result: any = await api.post(
          `projects/${job.project}/create_task/`,
          {
            task_name,
            context_job_uuid: job.uuid,
          }
        );
        if (created_job_result?.success && created_job_result.data?.new_job) {
          const created_job: Job = created_job_result.data.new_job;
          mutateJobs();
          router.push(`/ccp4i2/project/${job.project}/job/${created_job.id}`);
        } else {
          const errorMessage = created_job_result?.error || "Failed to create task";
          setMessage(`Failed to create task: ${errorMessage}`, "error");
        }
      } catch (error) {
        setMessage(`Error creating task: ${error instanceof Error ? error.message : String(error)}`, "error");
      }
    },
    [job, mutateJobs, api, router, setMessage]
  );

  // Extract result from new API format: {success: true, data: {result: [...]}}
  const whatNextResult = what_next?.success ? what_next.data?.result : null;

  const barRef = useRef<HTMLDivElement>(null);
  const [tier, setTier] = useState(TIER_FULL);
  // Bumped on every resize so the fitting pass below re-runs even when the tier
  // it lands on is the one already in state.
  const [resizeGeneration, setResizeGeneration] = useState(0);

  // On any width change, optimistically return to the widest tier; the fitting
  // pass then walks back down. This is what lets the bar recover its full
  // labels when the pane is widened again.
  useIsomorphicLayoutEffect(() => {
    const el = barRef.current;
    if (!el || typeof ResizeObserver === "undefined") return;
    const observer = new ResizeObserver(() => {
      setTier(TIER_FULL);
      setResizeGeneration((generation) => generation + 1);
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, [whatNextResult]);

  // Step down one tier per render until the content fits (or we run out of
  // tiers). Converges because the tier only ever decreases in width.
  useIsomorphicLayoutEffect(() => {
    const el = barRef.current;
    if (!el || tier >= NARROWEST_TIER) return;
    if (el.scrollWidth > el.clientWidth + 1) setTier(tier + 1);
  }, [tier, resizeGeneration, whatNextResult]);

  const labelFor = (task: any) => {
    if (tier >= TIER_ICON) return null;
    const shortLabel = task?.shortTitle || task?.taskName;
    if (tier >= TIER_SHORT) return shortLabel;
    return task?.title && task.title.length <= MAX_LABEL_CHARS
      ? task.title
      : shortLabel;
  };

  return (
    whatNextResult &&
    whatNextResult.length > 0 &&
    job?.status == 6 && (
      <Toolbar
        ref={barRef}
        variant="dense"
        sx={{
          width: "100%",
          minHeight: "48px",
          // Pinned to the bottom of the report pane: never shrink, and never
          // wrap. The overflow is a backstop only -- the tier logic above
          // should have shed enough label width to avoid ever needing it.
          flexShrink: 0,
          flexWrap: "nowrap",
          overflowX: "auto",
          scrollbarWidth: "thin",
          backgroundColor: customColors.ui.veryLightGray,
          borderTop: `1px solid ${customColors.ui.mediumGray}`,
          px: 2,
        }}
      >
        <Typography
          variant="h6"
          sx={{ fontWeight: "bold", mr: 3, whiteSpace: "nowrap" }}
        >
          What next:
        </Typography>
        {whatNextResult.map((task: any, iTask: number) => {
          const label = labelFor(task);
          return (
            <Tooltip key={iTask} title={task.title || task.taskName}>
              <Button
                variant="outlined"
                size="small"
                sx={{
                  ml: 1,
                  minWidth: "auto",
                  flexShrink: 0,
                  whiteSpace: "nowrap",
                  // Labels are sentence-style task titles, so leave the casing
                  // as authored rather than MUI's default uppercasing, which
                  // loses the distinctions the titles are drawing.
                  textTransform: "none",
                }}
                onClick={() => {
                  handleTaskSelect(task.taskName);
                }}
              >
                <Avatar
                  sx={{
                    width: 20,
                    height: 20,
                    mr: label ? 1 : 0,
                  }}
                  src={`/svgicons/${task.taskName}.svg`}
                  alt={task.taskName}
                >
                  {task.taskName?.[0]?.toUpperCase()}
                </Avatar>
                {label}
              </Button>
            </Tooltip>
          );
        })}
      </Toolbar>
    )
  );
};
