import { Avatar, Button, Stack, Toolbar, Typography } from "@mui/material";
import { useApi } from "../../api";
import { useCCP4i2Window } from "../../app-context";
import { useJob, useProject } from "../../utils";
import { Job } from "../../types/models";
import { useCallback } from "react";
import { useTheme } from "../../theme/theme-provider";
import { useRouter } from "next/navigation";
import { usePopcorn } from "../../providers/popcorn-provider";

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

  return (
    whatNextResult &&
    whatNextResult.length > 0 &&
    job?.status == 6 && (
      <Toolbar
        variant="dense"
        sx={{
          width: "100%",
          minHeight: "48px",
          backgroundColor: customColors.ui.veryLightGray,
          borderTop: `1px solid ${customColors.ui.mediumGray}`,
          px: 2,
        }}
      >
        <Typography variant="h6" sx={{ fontWeight: "bold", mr: 3 }}>
          What next:
        </Typography>
        {whatNextResult.map((task: any, iTask: number) => (
          <Button
            key={iTask}
            variant="outlined"
            size="small"
            sx={{ ml: 1, minWidth: "auto" }}
            onClick={() => {
              handleTaskSelect(task.taskName);
            }}
          >
            <Avatar
              sx={{
                width: 20,
                height: 20,
                mr: 1,
              }}
              src={`/svgicons/${task.taskName}.svg`}
              alt={`/qticons/${task.taskName}.png`}
            />
            {task.shortTitle}
          </Button>
        ))}
      </Toolbar>
    )
  );
};
