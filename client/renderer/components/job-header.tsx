/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import {
  Autocomplete,
  Button,
  LinearProgress,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import { Job } from "../types/models";
import EditableTypography from "./editable-typography";
import { useState } from "react";
import { CCP4i2JobAvatar } from "./job-avatar";
import { JobMenu, useJobMenu } from "../providers/job-context-menu";
import { Menu } from "@mui/icons-material";
import { useJob, useProjectJobs } from "../utils";
import { useApi } from "../api";
import { useDroppable } from "@dnd-kit/core";
import { useTheme } from "../theme/theme-provider";

interface JobHeaderProps {
  job: Job;
  mutateJobs: () => void;
  mutateJob?: () => void;
}
export const JobHeader: React.FC<JobHeaderProps> = ({ job, mutateJobs }) => {
  const { customColors } = useTheme();
  const [contextJob, setContextJob] = useState<Job | null>(null);
  const { setJobMenuAnchorEl, setJob } = useJobMenu();
  const api = useApi();

  const { container, mutateContainer } = useJob(job.id);

  const { jobs: project_jobs } = useProjectJobs(job.project);

  const { isOver, setNodeRef } = useDroppable({
    id: `job_${job.id}`,
    data: { job },
  });

  const handleContextJobChange = async (
    event: React.SyntheticEvent,
    value: Job | null
  ) => {
    setContextJob(value);
    const formData = { context_job_uuid: value ? value.uuid : null };
    const result = await api
      .post<Job>(`jobs/${job.id}/set_context_job`, formData)
      .then(() => {
        if (mutateContainer) {
          mutateContainer();
        }
      });
  };

  if (!job) return <LinearProgress />;

  return (
    <>
      <Stack
        direction="row"
        spacing={2}
        alignItems="center"
        ref={setNodeRef}
        sx={{
          backgroundColor: isOver
            ? customColors.ui.mediumGray
            : customColors.ui.veryLightGray,
          border: isOver ? `2px dashed ${customColors.ui.lightBlue}` : "none",
          transition: "background-color 0.3s, border 0.3s",
          px: 3,
          py: 1,
        }}
      >
        <CCP4i2JobAvatar job={job} />
        <Typography variant="h5">
          {job.number}
        </Typography>
        <EditableTypography
          variant="h5"
          text={job.title}
          onDelay={async (name) => {
            const formData = new FormData();
            formData.set("title", name);
            await api.patch(`jobs/${job.id}`, formData);
            mutateJobs();
          }}
          sx={{ flex: "auto" }}
        />
        {project_jobs && job?.status == 1 && (
          <Autocomplete
            onChange={handleContextJobChange}
            disabled={job.status !== 1}
            options={project_jobs.filter((j: Job) => j.parent == null)}
            value={contextJob}
            renderInput={(params) => (
              <TextField
                {...params}
                label="Take context from: "
                size="small"
              />
            )}
            getOptionLabel={(option: Job) =>
              `${option.number}:${option.title || option.task_name}`
            }
            sx={{ flex: "auto" }}
          />
        )}
        <Button
          variant="outlined"
          onClick={(ev) => {
            ev.stopPropagation();
            setJobMenuAnchorEl(ev.currentTarget);
            setJob(job);
          }}
        >
          <Menu />
        </Button>
      </Stack>
      <JobMenu />
    </>
  );
};
