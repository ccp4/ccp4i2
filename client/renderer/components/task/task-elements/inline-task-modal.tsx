/**
 * InlineTaskModal - Modal dialog for creating and running a helper task inline.
 *
 * Used when a data file widget needs to generate a missing data type
 * (e.g., FreeR flags, ASU contents) without navigating away from the
 * current task. Embeds the helper task's full interface inside a dialog,
 * runs it synchronously, and auto-selects the output file in the parent widget.
 */
import {
  Box,
  Button,
  CircularProgress,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Stack,
  Typography,
} from "@mui/material";
import { PlayArrow as RunIcon } from "@mui/icons-material";
import { useCallback, useEffect, useRef, useState } from "react";
import { useApi } from "../../../api";
import { apiJson } from "../../../api-fetch";
import { useJob, useProject, useProjectFiles } from "../../../utils";
import { Job, File as DjangoFile } from "../../../types/models";
import { TaskProvider } from "../../../providers/task-provider";
import { TaskContainer } from "../task-interfaces/task-container";

/** Configuration for which output parameter to extract from the helper task */
interface OutputConfig {
  /** The param name on the helper task's output (e.g., "FREEROUT", "ASUCONTENTFILE") */
  outputParamName: string;
}

/** Maps helper task names to their output configuration */
const HELPER_TASK_CONFIG: Record<string, OutputConfig> = {
  freerflag: { outputParamName: "FREEROUT" },
  ProvideAsuContents: { outputParamName: "ASUCONTENTFILE" },
};

export interface InlineTaskModalProps {
  /** Whether the modal is open */
  open: boolean;
  /** Called when the modal should close (cancel or completion) */
  onClose: () => void;
  /** The task name to create and run (e.g., "freerflag", "ProvideAsuContents") */
  taskName: string;
  /** The parent job from which the helper task is being launched */
  parentJob: Job;
  /** Called with the output file when the helper task completes successfully */
  onOutputFileReady: (file: DjangoFile) => void;
  /** Human-readable title for the dialog */
  title?: string;
}

type ModalPhase = "creating" | "configuring" | "running" | "succeeded" | "failed";

export const InlineTaskModal: React.FC<InlineTaskModalProps> = ({
  open,
  onClose,
  taskName,
  parentJob,
  onOutputFileReady,
  title,
}) => {
  const api = useApi();
  const { createPeerTask } = useJob(parentJob.id);
  const { mutateJobs } = useProject(parentJob.project);
  const { mutateFiles } = useProjectFiles(parentJob.project);

  const [phase, setPhase] = useState<ModalPhase>("creating");
  const [helperJob, setHelperJob] = useState<Job | null>(null);
  const [errorMessage, setErrorMessage] = useState<string>("");

  // Ref to track whether we've already started creating a job for this modal open.
  // Prevents the cascading re-render loop: createPeerTask calls mutateJobs(),
  // which re-renders hooks, which recreates createPeerTask, which would re-trigger useEffect.
  const isCreatingRef = useRef(false);

  // Stable ref for createPeerTask to avoid dependency on its identity
  const createPeerTaskRef = useRef(createPeerTask);
  createPeerTaskRef.current = createPeerTask;

  // Create the helper job when the modal opens
  useEffect(() => {
    if (!open) {
      // Reset state when closed
      setPhase("creating");
      setHelperJob(null);
      setErrorMessage("");
      isCreatingRef.current = false;
      return;
    }

    // Guard: only create once per modal open
    if (isCreatingRef.current) return;
    isCreatingRef.current = true;

    let cancelled = false;

    const createJob = async () => {
      try {
        setPhase("creating");
        const createdJob = await createPeerTaskRef.current(taskName);
        if (cancelled) return;

        if (!createdJob) {
          setPhase("failed");
          setErrorMessage(`Failed to create ${taskName} task`);
          return;
        }

        setHelperJob(createdJob);
        setPhase("configuring");
      } catch (error) {
        if (cancelled) return;
        setPhase("failed");
        setErrorMessage(`Error creating task: ${error}`);
      }
    };

    createJob();

    return () => {
      cancelled = true;
    };
  }, [open, taskName]);

  // Run the helper task synchronously
  const handleRun = useCallback(async () => {
    if (!helperJob) return;

    setPhase("running");
    setErrorMessage("");

    try {
      // Run synchronously via the run_local endpoint
      const result = await api.post<Job>(`jobs/${helperJob.id}/run_local/`, {
        synchronous: true,
      });

      // Check if the job succeeded (status 6 = Finished)
      if (result.status === 6) {
        // Refresh project files so the new output appears
        await Promise.all([mutateJobs(), mutateFiles()]);

        // Fetch the helper job's files to find the output
        const config = HELPER_TASK_CONFIG[taskName];
        if (!config) {
          setPhase("failed");
          setErrorMessage(`No output configuration for task: ${taskName}`);
          return;
        }

        // Fetch the helper job's output files via GET
        // The files endpoint returns a plain array (not wrapped in {success, data})
        const files: DjangoFile[] = await apiJson<DjangoFile[]>(
          `jobs/${helperJob.id}/files/`
        );

        const outputFile = files.find(
          (f) => f.job_param_name === config.outputParamName
        );

        if (outputFile) {
          setPhase("succeeded");
          onOutputFileReady(outputFile);
          onClose();
        } else {
          setPhase("failed");
          setErrorMessage(
            "Task completed but output file was not found. " +
            "The file may still be available in the project file list."
          );
          // Still refresh so user can manually select
          await Promise.all([mutateJobs(), mutateFiles()]);
        }
      } else {
        setPhase("failed");
        setErrorMessage(
          `Task finished with status: ${getStatusLabel(result.status)}. ` +
          "Check the task log for details."
        );
      }
    } catch (error) {
      setPhase("failed");
      setErrorMessage(`Error running task: ${error}`);
    }
  }, [helperJob, api, taskName, mutateJobs, mutateFiles, onOutputFileReady, onClose]);

  const dialogTitle = title || `Generate: ${taskName}`;

  return (
    <Dialog
      open={open}
      onClose={phase === "running" ? undefined : onClose}
      maxWidth="md"
      fullWidth
      disableEscapeKeyDown={phase === "running"}
    >
      <DialogTitle>{dialogTitle}</DialogTitle>

      <DialogContent dividers sx={{ minHeight: 300 }}>
        {phase === "creating" && (
          <Stack alignItems="center" justifyContent="center" sx={{ py: 4 }}>
            <CircularProgress size={40} />
            <Typography sx={{ mt: 2 }} color="text.secondary">
              Creating {taskName} task...
            </Typography>
          </Stack>
        )}

        {phase === "configuring" && helperJob && (
          <Box sx={{ minHeight: 200 }}>
            <TaskProvider>
              <TaskContainer jobId={helperJob.id} />
            </TaskProvider>
          </Box>
        )}

        {phase === "running" && (
          <Stack alignItems="center" justifyContent="center" sx={{ py: 4 }}>
            <CircularProgress size={60} />
            <Typography sx={{ mt: 2 }} color="text.secondary">
              Running {taskName}...
            </Typography>
            <Typography variant="caption" color="text.secondary" sx={{ mt: 1 }}>
              This should only take a few seconds.
            </Typography>
          </Stack>
        )}

        {phase === "succeeded" && (
          <Stack alignItems="center" justifyContent="center" sx={{ py: 4 }}>
            <Typography color="success.main" variant="h6">
              Done
            </Typography>
          </Stack>
        )}

        {phase === "failed" && (
          <Stack spacing={2} sx={{ py: 2 }}>
            <Typography color="error">
              {errorMessage || "An unexpected error occurred."}
            </Typography>
            {helperJob && (
              <Box sx={{ minHeight: 200 }}>
                <TaskProvider>
                  <TaskContainer jobId={helperJob.id} />
                </TaskProvider>
              </Box>
            )}
          </Stack>
        )}
      </DialogContent>

      <DialogActions>
        {phase === "configuring" && (
          <>
            <Button onClick={onClose}>Cancel</Button>
            <Button
              variant="contained"
              startIcon={<RunIcon />}
              onClick={handleRun}
            >
              Run
            </Button>
          </>
        )}
        {phase === "failed" && (
          <>
            <Button onClick={onClose}>Close</Button>
            {helperJob && (
              <Button
                variant="contained"
                startIcon={<RunIcon />}
                onClick={handleRun}
              >
                Retry
              </Button>
            )}
          </>
        )}
        {phase === "running" && (
          <Button disabled>Running...</Button>
        )}
      </DialogActions>
    </Dialog>
  );
};

/** Convert job status code to a human-readable label */
function getStatusLabel(status: number): string {
  switch (status) {
    case 1: return "Pending";
    case 2: return "Queued";
    case 3: return "Running";
    case 4: return "Failed";
    case 5: return "Unsatisfactory";
    case 6: return "Finished";
    case 7: return "Running remotely";
    default: return `Unknown (${status})`;
  }
}
