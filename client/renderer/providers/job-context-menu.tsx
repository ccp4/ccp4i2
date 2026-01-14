import React, {
  PropsWithChildren,
  SyntheticEvent,
  useCallback,
  useContext,
  useMemo,
  useState,
  useRef,
  useEffect,
} from "react";
import { Job, File as DjangoFile } from "../types/models";
import { doDownload, useApi } from "../api";
import { useRouter } from "next/navigation";
import { useDeleteDialog } from "./delete-dialog";
import {
  List,
  ListItem,
  Menu,
  MenuItem,
  Paper,
  Toolbar,
  Popper,
  TextField,
  ClickAwayListener,
  Box,
} from "@mui/material";
import { CCP4i2JobAvatar } from "../components/job-avatar";
import {
  CopyAll,
  Delete,
  FireExtinguisherRounded,
  Pending,
  Preview,
  RunCircle,
  SmsFailed,
  TerminalOutlined,
  Edit,
  Download,
  Refresh,
  FolderCopy,
} from "@mui/icons-material";
import { createContext } from "react";
import { usePopcorn } from "./popcorn-provider";
import { useRunCheck } from "./run-check-provider";
import { CCP4i2MoorhenIcon } from "../components/General/CCP4i2Icons";
import { useJob } from "../utils";
import ExportJobMenu from "../components/export-job-file-menu";
import { useRecentlyStartedJobs } from "./recently-started-jobs-context";
import { useSWRConfig } from "swr";

interface JobMenuContextDataProps {
  jobMenuAnchorEl: HTMLElement | null;
  setJobMenuAnchorEl: (element: HTMLElement | null) => void;
  job: Job | null;
  setJob: (job: Job | null) => void;
  fileExportJobId?: number | null;
  setFileExportJobId?: (jobId: number | null) => void;
}

export const JobMenuContext = createContext<JobMenuContextDataProps>({
  jobMenuAnchorEl: null,
  setJobMenuAnchorEl: () => {},
  job: null,
  setJob: () => {},
});

export const JobMenuProvider: React.FC<PropsWithChildren> = ({ children }) => {
  const [jobMenuAnchorEl, setJobMenuAnchorEl] = useState<HTMLElement | null>(
    null
  );
  const [job, setJob] = useState<Job | null>(null);
  const [fileExportJobId, setFileExportJobId] = useState<number | null>(null);

  return (
    <>
      <JobMenuContext.Provider
        value={{
          jobMenuAnchorEl,
          setJobMenuAnchorEl,
          job,
          setJob,
          fileExportJobId,
          setFileExportJobId,
        }}
      >
        {children}
        <JobMenu />
        <ExportJobMenu jobId={fileExportJobId} setJobId={setFileExportJobId} />
      </JobMenuContext.Provider>
    </>
  );
};

export interface JobWithChildren extends Job {
  children: (Job | DjangoFile)[];
}

export const JobMenu: React.FC = () => {
  const {
    jobMenuAnchorEl,
    setJobMenuAnchorEl,
    job,
    setJob,
    fileExportJobId,
    setFileExportJobId,
  } = useJobMenu();
  const api = useApi();
  const router = useRouter();
  const { setMessage } = usePopcorn();
  const { confirmTaskRun } = useRunCheck();
  const { markJobAsStarting } = useRecentlyStartedJobs();
  const { mutate: globalMutate } = useSWRConfig();

  const deleteDialog = useDeleteDialog();
  const [statusMenuAnchorEl, setStatusMenuAnchorEl] = useState<
    null | (EventTarget & Element)
  >(null);

  // State for annotation editor popper
  const [annotationPopperAnchorEl, setAnnotationPopperAnchorEl] =
    useState<HTMLElement | null>(null);
  const [annotationValue, setAnnotationValue] = useState<string>("");

  // Store the original title value when popper opens
  const originalAnnotationValue = useRef<string>("");

  // Ref to store the current timeout ID for cleanup
  const timeoutRef = useRef<NodeJS.Timeout | null>(null);

  const { data: jobs, mutate: mutateJobs } = api.get_endpoint<Job[]>({
    type: "projects",
    id: job?.project,
    endpoint: "jobs",
  });

  // Also get mutator for job_tree endpoint (used by classic-jobs-list)
  const { mutate: mutateJobTree } = api.get_endpoint<any>({
    type: "projects",
    id: job?.project,
    endpoint: "job_tree",
  });

  // Combined mutation function that invalidates both endpoints
  const mutateAllJobs = useCallback(() => {
    mutateJobs();
    mutateJobTree();
  }, [mutateJobs, mutateJobTree]);

  const { data: dependentJobs } = api.get_endpoint<Job[]>({
    type: "jobs",
    id: job?.id,
    endpoint: "dependent_jobs",
  });

  const topLevelDependentJobs = useMemo<Job[]>(() => {
    try {
      if (Array.isArray(dependentJobs) && dependentJobs.length > 0) {
        const result = dependentJobs
          ? dependentJobs.filter((job) => job.parent === null)
          : [];

        return result;
      }
    } catch (error) {
      console.log(dependentJobs);
      console.error(error);
    }
    return [];
  }, [dependentJobs]);

  // Function to save annotation immediately
  const saveAnnotation = useCallback(
    async (value: string) => {
      if (!job) return;

      try {
        console.log("Job annotation updated:", {
          jobId: job.id,
          jobNumber: job.number,
          jobTitle: job.title,
          newAnnotation: value,
        });

        const formData = new FormData();
        formData.append("title", value);

        await api.patch(`jobs/${job.id}/`, formData);
        mutateAllJobs(); // Refresh the jobs list (both jobs and job_tree endpoints)
      } catch (error) {
        console.error("Failed to update job annotation:", error);
      }
    },
    [job, api, mutateAllJobs]
  );

  // Handle Enter and Escape key presses in the text field
  const handleKeyDown = useCallback(
    async (event: React.KeyboardEvent<HTMLDivElement>) => {
      if (event.key === "Enter" && !event.shiftKey) {
        event.preventDefault();

        // Clear any pending timeout
        if (timeoutRef.current) {
          clearTimeout(timeoutRef.current);
          timeoutRef.current = null;
        }

        // Save immediately and close popper
        await saveAnnotation(annotationValue);
        setAnnotationPopperAnchorEl(null);
      } else if (event.key === "Escape") {
        event.preventDefault();

        // Clear any pending timeout
        if (timeoutRef.current) {
          clearTimeout(timeoutRef.current);
          timeoutRef.current = null;
        }

        // Restore original value and close popper without saving
        await saveAnnotation(originalAnnotationValue.current);
        setAnnotationPopperAnchorEl(null);
      }
    },
    [annotationValue, saveAnnotation]
  );

  // Handle annotation text change with delay
  const handleAnnotationChange = useCallback(
    (newValue: string) => {
      setAnnotationValue(newValue);

      // Clear existing timeout
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }

      // Create a delayed function to save annotation
      timeoutRef.current = setTimeout(() => {
        saveAnnotation(newValue);
        timeoutRef.current = null;
      }, 500); // 500ms delay
    },
    [saveAnnotation]
  );

  // Handle popper close
  const handleAnnotationPopperClose = useCallback(() => {
    // Clear any pending timeout when closing
    if (timeoutRef.current) {
      clearTimeout(timeoutRef.current);
      timeoutRef.current = null;
    }
    setAnnotationPopperAnchorEl(null);
  }, []);

  // Handle click away from popper
  const handleClickAway = useCallback(
    (event: MouseEvent | TouchEvent) => {
      handleAnnotationPopperClose();
    },
    [handleAnnotationPopperClose]
  );

  // Cleanup timeout on unmount
  useEffect(() => {
    return () => {
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }
    };
  }, []);

  // Handle edit annotation menu item click
  const handleEditAnnotation = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (job) {
        // Store the original title value for potential restoration
        const currentTitle = job.title || "";
        originalAnnotationValue.current = currentTitle;

        // Set current title/annotation value
        setAnnotationValue(currentTitle);
        // Use the menu anchor element as the popper anchor
        setAnnotationPopperAnchorEl(jobMenuAnchorEl);
        // Close the menu
        setJobMenuAnchorEl(null);
      }
    },
    [job, jobMenuAnchorEl, setJobMenuAnchorEl]
  );

  const handleClone = useCallback(
    async (ev: SyntheticEvent) => {
      if (!job) return;
      ev.stopPropagation();

      // Show subtle progress message
      setMessage(`Cloning job ${job.number}...`);

      try {
        const cloneResult: Job = await api.post(`jobs/${job.id}/clone/`);
        if (cloneResult?.id) {
          mutateAllJobs();
          setJobMenuAnchorEl(null);
          setMessage(`Job ${job.number} cloned successfully`);
          router.push(`/ccp4i2/project/${job.project}/job/${cloneResult.id}`);
        }
      } catch (error) {
        setMessage(`Failed to clone job ${job.number}`);
        console.error("Clone failed:", error);
      }
    },
    [job, mutateAllJobs, setMessage]
  );

  // Add state for run loading (add this near the other useState declarations)
  const [isRunning, setIsRunning] = useState(false);

  const handleRun = useCallback(
    async (ev: SyntheticEvent) => {
      if (!job) return;

      // Set running state immediately
      setIsRunning(true);
      setMessage(`Preparing to run job ${job.number}...`);

      try {
        setJobMenuAnchorEl(null);

        const confirmed = await confirmTaskRun(job.id);
        if (!confirmed) {
          setMessage(`Run cancelled for job ${job.number}`);
          return;
        }

        ev.stopPropagation();

        // Update message for actual submission
        setMessage(`Submitting job ${job.number}...`);

        const runResult: Job = await api.post(`jobs/${job.id}/run/`);

        setMessage(
          `Submitted job ${runResult?.number}: ${runResult?.task_name}`
        );

        if (runResult?.id) {
          // Mark job as recently started to trigger grace period polling
          // even if the database status hasn't caught up yet
          markJobAsStarting(
            runResult.id,
            job.project,
            runResult.number,
            runResult.task_name || job.title
          );
          mutateAllJobs();
          router.push(`/ccp4i2/project/${job.project}/job/${runResult.id}`);
        }
      } catch (error) {
        setMessage(`Failed to run job ${job.number}`);
        console.error("Run failed:", error);
      } finally {
        setIsRunning(false);
      }
    },
    [job, mutateAllJobs, setMessage, isRunning, markJobAsStarting]
  );

  const handleTerminalInJobDirectory = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (job) {
        api.post<any>(`jobs/${job.id}/preview/`, { viewer: "terminal" });
        setJobMenuAnchorEl(null);
      }
    },
    [job, setJobMenuAnchorEl]
  );

  const handleStatusClicked = useCallback(
    (ev: SyntheticEvent) => {
      if (!job) return;
      ev.stopPropagation();
      setStatusMenuAnchorEl(ev.currentTarget);
    },
    [job, setStatusMenuAnchorEl]
  );

  const handleOpenInNewWindow = useCallback(
    (ev: SyntheticEvent) => {
      if (!job) return;
      ev.stopPropagation();
      setJobMenuAnchorEl(null);
      const path = `/ccp4i2/moorhen-page/job-by-id/${job.id}`;
      window.open(path, "_blank", "noopener,noreferrer");
    },
    [job, setJobMenuAnchorEl]
  );

  const handleDelete = useCallback(
    (ev: SyntheticEvent) => {
      if (!job) return;
      ev.stopPropagation();
      if (deleteDialog)
        deleteDialog({
          type: "show",
          what: `${job.number}: ${job.title}`,
          onDelete: async () => {
            const deleteResult = await api.delete(`jobs/${job.id}`);
            console.log(deleteResult);
            mutateAllJobs();
            setJobMenuAnchorEl(null);
            setJob(null);
            router.push(`/ccp4i2/project/${job.project}`);
          },
          onCancel: () => {
            setJobMenuAnchorEl(null);
            setJob(null);
          },
          children: [
            <Paper
              key="dependentJobs"
              sx={{ maxHeight: "10rem", overflowY: "auto" }}
            >
              {topLevelDependentJobs && topLevelDependentJobs?.length > 0 && (
                <>
                  The following {topLevelDependentJobs.length} dependent jobs
                  would be deleted
                  <List dense>
                    {topLevelDependentJobs &&
                      topLevelDependentJobs.map((dependentJob: Job) => {
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
              dependentJobs.some(
                (dependentJob: Job) => dependentJob.status == 6
              ))
          ),
        });
    },
    [dependentJobs, job, mutateAllJobs]
  );

  const handleExportJob = useCallback(
    async (ev: SyntheticEvent) => {
      if (!job) return;
      ev.stopPropagation();

      try {
        const theUrl = api.noSlashUrl(`jobs/${job.id}/export_job/`);
        doDownload(theUrl, `job_${job.number}_export.zip`);
        setJobMenuAnchorEl(null);
        setMessage(`Job ${job.number} exported successfully`);
      } catch (error) {
        console.error("Failed to export job:", error);
        setMessage(`Failed to export job ${job.number}: ${error.message}`);
      }
    },
    [job, setJobMenuAnchorEl, setMessage]
  );

  const handleExportOptions = useCallback(
    (ev: SyntheticEvent) => {
      if (!job || !setFileExportJobId) return;
      ev.stopPropagation();

      // Set the job ID for the export dialog
      setFileExportJobId(job.id);

      // Close the context menu
      setJobMenuAnchorEl(null);
    },
    [job, setFileExportJobId, setJobMenuAnchorEl]
  );

  const handleRegenerateReport = useCallback(
    async (ev: SyntheticEvent) => {
      if (!job) return;
      ev.stopPropagation();

      setMessage(`Regenerating report for job ${job.number}...`);
      setJobMenuAnchorEl(null);

      try {
        await api.post(`jobs/${job.id}/regenerate_report/`);

        // Invalidate the report_xml SWR cache to trigger refetch
        await globalMutate(`/api/proxy/ccp4i2/jobs/${job.id}/report_xml/`);

        setMessage(`Report regenerated for job ${job.number}`);
      } catch (error) {
        console.error("Failed to regenerate report:", error);
        setMessage(`Failed to regenerate report for job ${job.number}`);
      }
    },
    [job, api, globalMutate, setMessage, setJobMenuAnchorEl]
  );

  const handleCopyDirectory = useCallback(
    async (ev: SyntheticEvent) => {
      if (!job) return;
      ev.stopPropagation();

      try {
        // If directory is already available on the job object, use it
        if (job.directory) {
          await navigator.clipboard.writeText(job.directory);
          setMessage(`Copied to clipboard: ${job.directory}`);
          setJobMenuAnchorEl(null);
          return;
        }

        // Otherwise fetch the full job data to get the directory property
        const response = await fetch(`/api/proxy/ccp4i2/jobs/${job.id}/`);
        const jobData = await response.json();

        if (jobData?.directory) {
          await navigator.clipboard.writeText(jobData.directory);
          setMessage(`Copied to clipboard: ${jobData.directory}`);
        } else {
          setMessage(`Job directory not available`);
        }
      } catch (error) {
        console.error("Failed to copy directory path:", error);
        setMessage(`Failed to copy directory path`);
      }

      setJobMenuAnchorEl(null);
    },
    [job, setMessage, setJobMenuAnchorEl]
  );

  return (
    job && (
      <>
        <Menu
          open={Boolean(jobMenuAnchorEl)}
          anchorEl={jobMenuAnchorEl}
          onClose={() => setJobMenuAnchorEl(null)}
        >
          <MenuItem
            key="Clone"
            disabled={job.number.includes(".")}
            onClick={handleClone}
          >
            <CopyAll sx={{ mr: 1 }} /> Clone
          </MenuItem>
          <MenuItem
            key="Run"
            disabled={job.number.includes(".") || job.status !== 1 || isRunning}
            onClick={handleRun}
          >
            <RunCircle sx={{ mr: 1, opacity: isRunning ? 0.5 : 1 }} />{" "}
            {isRunning ? "Running..." : "Run"}
          </MenuItem>
          <MenuItem key="EditAnnotation" onClick={handleEditAnnotation}>
            <Edit sx={{ mr: 1 }} /> Edit title
          </MenuItem>
          <MenuItem key="ExportJob" onClick={handleExportJob}>
            <Download sx={{ mr: 1 }} /> Export job
          </MenuItem>
          <MenuItem
            key="ExportOptions"
            disabled={job.status !== 6}
            onClick={handleExportOptions}
          >
            <Download sx={{ mr: 1 }} /> Export options
          </MenuItem>
          <MenuItem
            key="RegenerateReport"
            disabled={job.status === 1}
            onClick={handleRegenerateReport}
          >
            <Refresh sx={{ mr: 1 }} /> Regenerate report
          </MenuItem>
          <MenuItem
            key="Delete"
            disabled={job.number.includes(".")}
            onClick={handleDelete}
          >
            <Delete sx={{ mr: 1 }} /> Delete
          </MenuItem>
          <MenuItem
            key="Moorhen"
            disabled={job.status != 6}
            onClick={handleOpenInNewWindow}
          >
            <CCP4i2MoorhenIcon sx={{ mr: 1 }} /> Moorhen
          </MenuItem>
          <MenuItem
            key="Status"
            //disabled={job.number.includes(".")} // There are cases where we want to set the status of a job that is not top level
            onClick={handleStatusClicked}
          >
            <Delete sx={{ mr: 1 }} /> Set status
          </MenuItem>
          <MenuItem
            key="Terminal"
            disabled={false}
            onClick={handleTerminalInJobDirectory}
          >
            <TerminalOutlined sx={{ mr: 1 }} /> Terminal in job directory
          </MenuItem>
          <MenuItem
            key="CopyDirectory"
            disabled={false}
            onClick={handleCopyDirectory}
          >
            <FolderCopy sx={{ mr: 1 }} /> Copy directory path
          </MenuItem>
        </Menu>

        {/* Annotation Editor Popper */}
        <Popper
          open={Boolean(annotationPopperAnchorEl)}
          anchorEl={annotationPopperAnchorEl}
          placement="bottom-start"
          sx={{ zIndex: 1300 }}
        >
          <ClickAwayListener onClickAway={handleClickAway}>
            <Paper
              elevation={8}
              sx={{
                p: 2,
                minWidth: 300,
                maxWidth: 500,
              }}
            >
              <Box>
                <TextField
                  fullWidth
                  multiline
                  rows={3}
                  label="Job Title"
                  value={annotationValue}
                  onChange={(e) => handleAnnotationChange(e.target.value)}
                  onKeyDown={handleKeyDown}
                  placeholder="Enter title for this job..."
                  variant="outlined"
                  size="small"
                  autoFocus
                  helperText="Press Enter to save and close, Escape to cancel, or wait 500ms for auto-save"
                />
              </Box>
            </Paper>
          </ClickAwayListener>
        </Popper>

        <StatusMenu
          job={job}
          anchorEl={statusMenuAnchorEl}
          onClose={() => {
            setStatusMenuAnchorEl(null);
            setJobMenuAnchorEl(null);
          }}
        />
      </>
    )
  );
};

interface StatusMenuProps {
  job: Job;
  anchorEl: null | (EventTarget & Element);
  onClose: () => void;
}

const StatusMenu: React.FC<StatusMenuProps> = ({ job, anchorEl, onClose }) => {
  const api = useApi();

  // Ensure we have a valid project ID for the mutation
  const projectId = job?.project;

  const { mutate: mutateJobs } = api.get_endpoint<Job[]>({
    type: "projects",
    id: projectId,
    endpoint: "jobs",
  });

  const { mutate: mutateJobTree } = api.get_endpoint<any>({
    type: "projects",
    id: projectId,
    endpoint: "job_tree",
  });

  const mutateAllJobs = useCallback(() => {
    console.log(`StatusMenu: mutating jobs for project ${projectId}`);
    mutateJobs();
    mutateJobTree();
  }, [mutateJobs, mutateJobTree, projectId]);

  const setStatus = useCallback(
    async (status: number) => {
      if (!projectId) {
        console.error("StatusMenu: Cannot update status - no project ID");
        return;
      }

      const formData = new FormData();
      formData.append("status", `${status}`);

      console.log(`StatusMenu: Setting job ${job.id} status to ${status}`);
      const result = await api.patch(`jobs/${job.id}/`, formData);
      console.log("StatusMenu: Patch result:", result);

      if (result) {
        mutateAllJobs();
        onClose();
      }
    },
    [job, projectId, mutateAllJobs, onClose, api]
  );

  return (
    <Menu open={Boolean(anchorEl)} anchorEl={anchorEl} onClose={onClose}>
      <MenuItem onClick={() => setStatus(5)}>
        <SmsFailed sx={{ mr: 1 }} /> Failed
      </MenuItem>
      <MenuItem onClick={() => setStatus(6)}>
        <FireExtinguisherRounded sx={{ mr: 1 }} /> Finished
      </MenuItem>
      <MenuItem onClick={() => setStatus(1)}>
        <Pending sx={{ mr: 1 }} /> Pending
      </MenuItem>
    </Menu>
  );
};

export const useJobMenu = () => {
  const context = useContext(JobMenuContext);
  if (!context) {
    throw new Error("useJobMenu must be used within a JobMenuProvider");
  }
  return context;
};
