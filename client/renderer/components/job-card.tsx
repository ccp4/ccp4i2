import { useApi } from "../api";
import {
  Button,
  Card,
  CardActions,
  CardContent,
  CardHeader,
  Chip,
  Collapse,
  Grid2,
  List,
  ListItem,
  Menu,
  MenuItem,
  Paper,
  styled,
  Toolbar,
} from "@mui/material";
import { Job } from "../types/models";
import { useCallback, useMemo, useState } from "react";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import { MyExpandMore } from "./expand-more";
import { JobsGrid } from "./jobs-grid";
import FilesTable from "./files-table";
import {
  CopyAll,
  Delete,
  Menu as MenuIcon,
  RunCircle,
} from "@mui/icons-material";
import { useRouter } from "next/navigation";
import { JobHeader } from "./job-header";
import { useDeleteDialog } from "../providers/delete-dialog";
import { CCP4i2JobAvatar } from "./job-avatar";
import { useProject } from "../utils";
import { usePopcorn } from "../providers/popcorn-provider";
import { useRunCheck } from "../providers/run-check-provider";
import { useCCP4i2Window } from "../app-context";
import { useRecentlyStartedJobs } from "../providers/recently-started-jobs-context";

const MyCard = styled(Card)(({ theme }) => ({
  padding: theme.spacing(1),
  borderRadius: 20,
}));

/** Embedded KPIs structure from job_tree endpoint */
interface EmbeddedKPIs {
  float_values: Record<string, number>;
  char_values: Record<string, string>;
}

interface JobCardProps {
  job: Job;
  /** Embedded KPIs from job_tree endpoint (preferred) */
  kpis?: EmbeddedKPIs;
  withSubtitle?: boolean;
}
export const JobCard: React.FC<JobCardProps> = ({
  job,
  kpis,
  withSubtitle = false,
}) => {
  const api = useApi();
  const router = useRouter();
  const deleteDialog = useDeleteDialog();
  const { jobId, setJobId } = useCCP4i2Window();
  const projectId = useMemo(() => {
    return job.project;
  }, [job]);
  const {
    jobs,
    mutateJobs,
    mutateAllJobs,
    files,
    mutateFiles,
  } = useProject(job.project);

  const { confirmTaskRun } = useRunCheck();
  const { markJobAsStarting } = useRecentlyStartedJobs();

  const subJobs: any[] | undefined = useMemo(() => {
    return jobs?.filter((aJob) => aJob.parent === job.id);
  }, [jobs, job]);

  const jobFiles: any[] | undefined = useMemo(() => {
    return files?.filter((aFile) => aFile.job === job.id);
  }, [files, job]);

  const dependentJobs = useMemo(() => {
    if (jobs && job) {
      return jobs.filter(
        (possible_child: Job) => possible_child.parent == job.id,
        []
      );
    }
  }, [jobs]);

  const kpiContent = useMemo(() => {
    if (!kpis) return null;

    const charChips = Object.entries(kpis.char_values || {}).map(
      ([key, value]) => (
        <Chip
          key={`char_${key}`}
          avatar={<div style={{ width: "5rem" }}>{key}</div>}
          label={value}
        />
      )
    );

    const floatChips = Object.entries(kpis.float_values || {}).map(
      ([key, value]) => (
        <Chip
          key={`float_${key}`}
          sx={{ backgroundColor: "#DFD" }}
          avatar={<div style={{ width: "5rem" }}>{key}</div>}
          label={value.toPrecision(3)}
        />
      )
    );

    return (
      <>
        {charChips}
        {floatChips}
      </>
    );
  }, [kpis]);

  const [jobsExpanded, setJobsExpanded] = useState(false);
  const [filesExpanded, setFilesExpanded] = useState(false);
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const isMenuOpen = Boolean(anchorEl);
  const { setMessage } = usePopcorn();
  const handleMenuOpen = (event: React.MouseEvent<HTMLElement>) => {
    event.stopPropagation();
    setAnchorEl(event.currentTarget);
  };
  const handleMenuClose = () => {
    setAnchorEl(null);
  };
  const handleClone = async () => {
    try {
      const cloneResult: any = await api.post(`jobs/${job.id}/clone/`);
      if (cloneResult?.success === false) {
        setMessage(`Failed to clone job: ${cloneResult?.error || "Unknown error"}`, "error");
        return;
      }
      if (cloneResult?.id) {
        mutateAllJobs();
        setAnchorEl(null);
        router.push(`/project/${projectId}/job/${cloneResult.id}`);
      }
    } catch (error) {
      setMessage(`Error cloning job: ${error instanceof Error ? error.message : String(error)}`, "error");
    }
  };
  const handleRun = async () => {
    setAnchorEl(null);
    const confirmed = await confirmTaskRun(job.id);
    if (!confirmed) {
      return;
    }
    try {
      const runResult: any = await api.post(`jobs/${job.id}/run/`);
      if (runResult?.success === false) {
        setMessage(`Failed to run job: ${runResult?.error || "Unknown error"}`, "error");
        return;
      }
      setMessage(`Submitted job ${runResult?.number}: ${runResult?.task_name}`, "success");
      if (runResult?.id) {
        // Mark job as recently started to trigger grace period polling
        markJobAsStarting(
          runResult.id,
          projectId,
          runResult.number,
          runResult.task_name || job.title
        );
        mutateAllJobs();
        router.push(`/project/${projectId}/job/${runResult.id}`);
      }
    } catch (error) {
      setMessage(`Error running job: ${error instanceof Error ? error.message : String(error)}`, "error");
    }
  };

  const handleDelete = useCallback(() => {
    if (deleteDialog)
      deleteDialog({
        type: "show",
        what: `${job.number}: ${job.title}`,
        onDelete: () => {
          api.delete(`jobs/${job.id}`).then(() => {
            mutateAllJobs();
            if (setJobId && jobId === job.id) setJobId(null);
          });
        },
        children: [
          <Paper sx={{ maxHeight: "10rem", overflowY: "auto" }}>
            {dependentJobs && dependentJobs?.length > 0 && (
              <>
                The following {dependentJobs.length} dependent jobs would be
                deleted
                <List dense>
                  {dependentJobs &&
                    dependentJobs.map((dependentJob: Job) => {
                      return (
                        <ListItem key={dependentJob.uuid}>
                          <Toolbar>
                            <CCP4i2JobAvatar job={dependentJob} />
                            {`${dependentJob.number}: ${dependentJob.title}`}
                          </Toolbar>
                        </ListItem>
                      );
                    })}
                </List>
              </>
            )}
          </Paper>,
        ],
        deleteDisabled: !(
          (dependentJobs && dependentJobs?.length == 0) ||
          (dependentJobs &&
            dependentJobs.some((dependentJob: Job) => dependentJob.status == 6))
        ),
      });
  }, [dependentJobs, mutateAllJobs]);

  const renderMenu = (
    <Menu
      anchorEl={anchorEl}
      anchorOrigin={{
        vertical: "top",
        horizontal: "right",
      }}
      id={`context-menu-${job.id}`}
      keepMounted
      transformOrigin={{
        vertical: "top",
        horizontal: "right",
      }}
      open={isMenuOpen}
      onClose={handleMenuClose}
    >
      <MenuItem key="Clone" onClick={handleClone}>
        <CopyAll /> Clone
      </MenuItem>
      <MenuItem key="Run" onClick={handleRun}>
        <RunCircle /> Run
      </MenuItem>
      <MenuItem key="Delete" onClick={handleDelete}>
        <Delete /> Delete
      </MenuItem>
    </Menu>
  );

  const handleExpandJobsClick = (ev: any) => {
    ev.stopPropagation();
    setJobsExpanded(!jobsExpanded);
  };
  const handleExpandFilesClick = (ev: any) => {
    ev.stopPropagation();
    setFilesExpanded(!filesExpanded);
  };
  return (
    <>
      <MyCard
        key={job.number}
        variant="elevation"
        onClick={(ev) => {
          ev.stopPropagation();
          router.push(`/project/${job.project}/job/${job.id}`);
        }}
      >
        <CardHeader
          action={
            <Button onClick={handleMenuOpen}>
              <MenuIcon />
            </Button>
          }
          sx={{ my: 0, mx: 0, px: 0, py: 0 }}
          title={<JobHeader job={job} mutateJobs={mutateJobs} />}
          subheader={kpiContent}
        />
        <CardActions sx={{ p: 0.5 }}>
          <Grid2 container>
            {subJobs && subJobs.length > 0 && (
              <Grid2 size={{ xs: 6 }}>
                Child jobs
                <MyExpandMore
                  expand={jobsExpanded}
                  onClick={handleExpandJobsClick}
                  aria-expanded={jobsExpanded}
                  aria-label="Show child jobs"
                >
                  <ExpandMoreIcon />
                </MyExpandMore>
              </Grid2>
            )}
            {jobFiles && jobFiles.length > 0 && (
              <Grid2 size={{ xs: 6 }}>
                Files
                <MyExpandMore
                  expand={filesExpanded}
                  onClick={handleExpandFilesClick}
                  aria-expanded={filesExpanded}
                  aria-label="Show files"
                >
                  <ExpandMoreIcon />
                </MyExpandMore>
              </Grid2>
            )}
          </Grid2>
        </CardActions>
        <Collapse
          key="ChildJobs"
          in={jobsExpanded}
          timeout="auto"
          unmountOnExit
        >
          <CardContent sx={{ p: 0.5 }}>
            {subJobs?.length && subJobs.length > 0 && (
              <JobsGrid
                projectId={projectId}
                size={{ xs: 12 }}
                parent={job.id}
                withSubtitles={false}
              />
            )}
          </CardContent>
        </Collapse>
        <Collapse
          key="JobFiles"
          in={filesExpanded}
          timeout="auto"
          unmountOnExit
        >
          <CardContent sx={{ p: 0.5 }}>
            {jobFiles && jobFiles.length > 0 && (
              <FilesTable files={jobFiles} mutate={mutateFiles} />
            )}
          </CardContent>
        </Collapse>
      </MyCard>
      {renderMenu}
    </>
  );
};
