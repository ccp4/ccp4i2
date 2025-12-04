"use client";
import { use, useCallback } from "react";
import { CircularProgress } from "@mui/material";
import { useApi } from "../../../api";
import { Job } from "../../../types/models";
import { useProject } from "../../../utils";
import { CCP4i2TaskTree } from "../../../components/task/task-chooser";
import { useRouter } from "next/navigation";

export default function DashboardPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const api = useApi();
  const { id } = use(params);
  const { project, jobs, mutateJobs } = useProject(parseInt(id));
  const router = useRouter();
  const handleTaskSelect = useCallback(
    async (task_name: string) => {
      const created_job_result: any = await api.post(
        `projects/${id}/create_task/`,
        {
          task_name,
        }
      );
      if (created_job_result?.success && created_job_result.data?.new_job) {
        const created_job: Job = created_job_result.data.new_job;
        mutateJobs();
        router.push(`/project/${id}/job/${created_job.id}`);
      }
    },
    [id, mutateJobs]
  );

  return project ? (
    <CCP4i2TaskTree onTaskSelect={handleTaskSelect} />
  ) : (
    <CircularProgress variant="indeterminate" />
  );
}
