"use client";
import { use, useCallback } from "react";
import { CircularProgress } from "@mui/material";
import { useApi } from "../../../../api";
import { Job } from "../../../../types/models";
import { useProject } from "../../../../utils";
import { CCP4i2TaskTree } from "../../../../components/task/task-chooser";
import { useRouter } from "next/navigation";
import { usePopcorn } from "../../../../providers/popcorn-provider";

export default function DashboardPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const api = useApi();
  const { id } = use(params);
  const { project, jobs, mutateJobs } = useProject(parseInt(id));
  const router = useRouter();
  const { setMessage } = usePopcorn();

  const handleTaskSelect = useCallback(
    async (task_name: string) => {
      try {
        const created_job_result: any = await api.post(
          `projects/${id}/create_task/`,
          {
            task_name,
          }
        );
        if (created_job_result?.success && created_job_result.data?.new_job) {
          const created_job: Job = created_job_result.data.new_job;
          mutateJobs();
          router.push(`/ccp4i2/project/${id}/job/${created_job.id}`);
        } else {
          // API returned an error response
          const errorMessage = created_job_result?.error || "Failed to create task";
          console.error("Task creation failed:", errorMessage);
          setMessage(`Failed to create task: ${errorMessage}`, "error");
        }
      } catch (error) {
        // Network or other error
        console.error("Task creation error:", error);
        setMessage(`Error creating task: ${error instanceof Error ? error.message : String(error)}`, "error");
      }
    },
    [id, mutateJobs, api, router, setMessage]
  );

  return project ? (
    <CCP4i2TaskTree onTaskSelect={handleTaskSelect} />
  ) : (
    <CircularProgress variant="indeterminate" />
  );
}
