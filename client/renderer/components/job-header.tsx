import {
  Autocomplete,
  Button,
  LinearProgress,
  Stack,
  TextField,
  Toolbar,
  Typography,
} from "@mui/material";
import { Job } from "../types/models";
import EditableTypography from "./editable-typography";
import { useApi } from "../api";
import { KeyedMutator } from "swr";
import { useState } from "react";
import { CCP4i2JobAvatar } from "./job-avatar";
import { JobMenu, useJobMenu } from "../providers/job-context-menu";
import { Menu } from "@mui/icons-material";
import { useJob } from "../utils";
import { useDroppable } from "@dnd-kit/core";
import { useTheme } from "../theme/theme-provider";

interface JobHeaderProps {
  job: Job;
  mutateJobs: KeyedMutator<Job[]>;
  mutateJob?: KeyedMutator<Job>;
}
export const JobHeader: React.FC<JobHeaderProps> = ({ job, mutateJobs }) => {
  const { customColors } = useTheme();
  const [contextJob, setContextJob] = useState<Job | null>(null);
  const { setJobMenuAnchorEl, setJob } = useJobMenu();
  const api = useApi();

  const { container, mutateContainer } = useJob(job.id);

  const { data: project_jobs } = api.get_endpoint<Job[]>({
    type: "projects",
    id: job.project,
    endpoint: "jobs",
  });

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
    console.log({ context_job_response: result });
  };

  if (!job) return <LinearProgress />;

  return (
    <>
      <Toolbar
        variant="regular"
        ref={setNodeRef}
        sx={{
          backgroundColor: isOver
            ? customColors.ui.mediumGray
            : customColors.ui.veryLightGray,
          border: isOver ? `2px dashed ${customColors.ui.lightBlue}` : "none",
          transition: "background-color 0.3s, border 0.3s",
        }}
      >
        <CCP4i2JobAvatar job={job} />
        <Stack direction="column" spacing={1} sx={{ ml: 2 }}>
          <Stack direction="row" spacing={1} alignItems="center">
            <Typography variant="h5" sx={{ ml: 2, mr: 2 }}>
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
            />
          </Stack>
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
                  sx={{
                    "& .MuiInputBase-root": {
                      padding: "2px",
                      fontSize: "0.875rem",
                    },
                    "& .MuiInputLabel-root": {
                      fontSize: "0.75rem",
                      lineHeight: "1rem",
                    },
                  }}
                />
              )}
              getOptionLabel={(option: Job) =>
                `${option.number}:${option.title || option.task_name}`
              }
              sx={{
                "& .MuiAutocomplete-inputRoot": {
                  padding: "2px !important",
                  minHeight: "30px",
                },
                "& .MuiAutocomplete-endAdornment": {
                  top: "calc(50% - 12px)",
                },
              }}
            />
          )}
        </Stack>
        <Typography sx={{ flexGrow: 1 }} />
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
      </Toolbar>
      <JobMenu />
    </>
  );
};
