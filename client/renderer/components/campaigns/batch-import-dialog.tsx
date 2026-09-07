"use client";

import { useCallback, useReducer, useState } from "react";
import {
  Alert,
  AlertTitle,
  Box,
  Button,
  Chip,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  LinearProgress,
  Paper,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  CloudUpload as UploadIcon,
  PlayArrow as RunIcon,
  Clear as ClearIcon,
} from "@mui/icons-material";
import { useCampaignsApi, useSmilesLookup } from "../../lib/campaigns-api";
import { DropZone } from "../common/drop-zone";
import {
  apiPost,
  apiFetch,
  apiUpload,
  apiUploadWithProgress,
} from "../../api-fetch";
import {
  BatchFileItem,
  BatchFileStatus,
  parseDatasetFilename,
} from "../../types/campaigns";
import { usePopcorn } from "../../providers/popcorn-provider";

const STATUS_META: Record<
  BatchFileStatus,
  { label: string; color: "default" | "primary" | "success" | "error" | "warning" | "info" }
> = {
  idle: { label: "Ready", color: "default" },
  creating_project: { label: "Creating project...", color: "info" },
  creating_job: { label: "Creating job...", color: "info" },
  uploading_coords: { label: "Uploading coords...", color: "info" },
  uploading_reflections: { label: "Uploading data...", color: "info" },
  fetching_smiles: { label: "Fetching SMILES...", color: "info" },
  setting_params: { label: "Setting params...", color: "info" },
  queuing: { label: "Queuing job...", color: "info" },
  done: { label: "Done", color: "success" },
  error: { label: "Error", color: "error" },
};

/** What each step was trying to do, for failure messages. */
const STEP_DESCRIPTION: Record<BatchFileStatus, string> = {
  idle: "starting",
  creating_project: "creating the project",
  creating_job: "creating the job",
  uploading_coords: "uploading the reference coordinates",
  uploading_reflections: "uploading the reflection file",
  fetching_smiles: "looking up the compound SMILES",
  setting_params: "setting the job parameters",
  queuing: "queuing the job",
  done: "finishing",
  error: "processing",
};

function formatBytes(bytes: number): string {
  if (!Number.isFinite(bytes) || bytes < 0) return "?";
  if (bytes < 1024) return `${bytes} B`;
  const units = ["kB", "MB", "GB"];
  let value = bytes / 1024;
  let unit = 0;
  while (value >= 1024 && unit < units.length - 1) {
    value /= 1024;
    unit += 1;
  }
  return `${value.toFixed(value < 10 ? 1 : 0)} ${units[unit]}`;
}

/**
 * Turn whatever went wrong into something a person can act on.
 *
 * A bare transport error ("signal is aborted without reason") tells the
 * user nothing about which of the eight steps failed, and an aborted
 * upload in particular reads as if nothing happened at all.
 */
function describeFailure(step: BatchFileStatus, error: unknown): string {
  const detail =
    error instanceof DOMException && error.name === "AbortError"
      ? "the request was aborted before the server replied"
      : error instanceof Error
        ? error.message
        : "unknown error";
  return `Failed while ${STEP_DESCRIPTION[step]}: ${detail}`;
}

interface BatchImportDialogProps {
  open: boolean;
  onClose: () => void;
  onSuccess?: () => void;
  campaignId: number;
  parentProjectId: number;
  latestCoordsFileId: number;
}

// Reducer for managing file state
type FileAction =
  | { type: "ADD"; file: File }
  | { type: "UPDATE_STATUS"; file: File; status: BatchFileStatus; error?: string }
  | { type: "UPDATE_PROGRESS"; file: File; loaded: number; total: number }
  | { type: "UPDATE_SMILES"; nclId: string; smiles: string }
  | { type: "SET_PROJECT_ID"; file: File; projectId?: number; jobId?: number }
  | { type: "CLEAR" };

function filesReducer(state: BatchFileItem[], action: FileAction): BatchFileItem[] {
  switch (action.type) {
    case "ADD": {
      const parsed = parseDatasetFilename(action.file.name);
      return [
        ...state,
        {
          file: action.file,
          ...parsed,
          status: "idle",
        },
      ];
    }
    case "UPDATE_STATUS": {
      // Entering a step clears the previous step's byte counter; the
      // upload steps repopulate it as data moves.
      return state.map((item) =>
        item.file === action.file
          ? {
              ...item,
              status: action.status,
              error: action.error,
              progress: undefined,
            }
          : item
      );
    }
    case "UPDATE_PROGRESS": {
      return state.map((item) =>
        item.file === action.file
          ? { ...item, progress: { loaded: action.loaded, total: action.total } }
          : item
      );
    }
    case "UPDATE_SMILES": {
      return state.map((item) =>
        item.nclId === action.nclId ? { ...item, smiles: action.smiles } : item
      );
    }
    case "SET_PROJECT_ID": {
      // Recorded as soon as each id exists, not once both do: a job that
      // fails to create must not lose the project that was created for
      // it, or the retry makes a second one.
      return state.map((item) =>
        item.file === action.file
          ? {
              ...item,
              ...(action.projectId !== undefined
                ? { projectId: action.projectId }
                : {}),
              ...(action.jobId !== undefined ? { jobId: action.jobId } : {}),
            }
          : item
      );
    }
    case "CLEAR":
      return [];
    default:
      return state;
  }
}

export function BatchImportDialog({
  open,
  onClose,
  onSuccess,
  campaignId,
  parentProjectId,
  latestCoordsFileId,
}: BatchImportDialogProps) {
  const [files, dispatch] = useReducer(filesReducer, []);
  const [processing, setProcessing] = useState(false);
  const [currentFileIndex, setCurrentFileIndex] = useState(-1);

  const campaignsApi = useCampaignsApi();
  const { setMessage } = usePopcorn();

  // Extract unique regIds for SMILES lookup
  const regIds = files
    .map((f) => (f.nclId ? parseInt(f.nclId) : null))
    .filter((id): id is number => id !== null && id !== 0 && !isNaN(id));
  const uniqueRegIds = [...new Set(regIds)];
  const { smilesMap } = useSmilesLookup(uniqueRegIds);

  // Handle file drop/select
  const handleFilesSelected = useCallback(
    (selectedFiles: FileList | null) => {
      if (!selectedFiles) return;

      Array.from(selectedFiles).forEach((file) => {
        // Only accept MTZ/SCA/CIF files
        const ext = file.name.split(".").pop()?.toLowerCase();
        if (ext === "mtz" || ext === "sca" || ext === "cif") {
          dispatch({ type: "ADD", file });
        }
      });
    },
    []
  );

  // Process a single file
  const processFile = async (item: BatchFileItem): Promise<void> => {
    const { file, nclId } = item;

    // A retry picks up whatever the failed attempt already created.
    // Re-running these steps blind would add a second member project for
    // the same dataset every time Process was clicked, leaving the
    // campaign full of empty projects to unpick by hand.
    let newProjectId = item.projectId;
    let newJobId = item.jobId;

    // The step in flight, so a failure can name what it was doing rather
    // than surfacing a bare transport error.
    let step: BatchFileStatus = "idle";
    const setStatus = (status: BatchFileStatus) => {
      step = status;
      dispatch({ type: "UPDATE_STATUS", file, status });
    };

    try {
      // 1. Create new project
      if (!newProjectId) {
        setStatus("creating_project");
        const projectName = file.name.replace(/\.[^.]+$/, "").replace(/\s/g, "_");
        const projectResponse = await apiPost<{ id: number }>("projects/", {
          name: projectName,
        });
        newProjectId = projectResponse.id;
        dispatch({ type: "SET_PROJECT_ID", file, projectId: newProjectId });

        // Add project to campaign as member
        await campaignsApi.addMember(campaignId, newProjectId, "member");
      }

      // 2. Create SubstituteLigand job
      // The create_task endpoint uses api_success() which wraps the response
      if (!newJobId) {
        setStatus("creating_job");
        const jobResponse = await apiPost<{
          success: boolean;
          data: { new_job: { id: number; uuid: string } };
        }>(
          `projects/${newProjectId}/create_task/`,
          {
            task_name: "SubstituteLigand",
            title: `Process ${file.name}`,
          }
        );
        newJobId = jobResponse.data.new_job.id;
        dispatch({ type: "SET_PROJECT_ID", file, jobId: newJobId });
      }

      // 3. Upload reference coordinates using upload_file_param endpoint
      setStatus("uploading_coords");
      // Get coords file content (id-addressed: files/{id}/download/, not files_by_uuid/{uuid}/download/)
      // Use apiFetch for authenticated download
      const coordsResponse = await apiFetch(
        `files/${latestCoordsFileId}/download/`
      );
      const coordsBlob = await coordsResponse.blob();
      // Extract filename from Content-Disposition header, handling both quoted and unquoted forms
      // e.g., 'filename="XYZOUT.pdb"' or 'filename=XYZOUT.pdb'
      const coordsFileName = coordsResponse.headers.get("Content-Disposition")
        ?.match(/filename="?([^"]+)"?/)?.[1] || "reference.pdb";

      const coordsFormData = new FormData();
      coordsFormData.append("file", coordsBlob, coordsFileName);
      coordsFormData.append("object_path", "SubstituteLigand.inputData.XYZIN");

      // Use apiUpload for authenticated upload
      await apiUpload(`jobs/${newJobId}/upload_file_param/`, coordsFormData);

      // 4. Set SMILES if we have a compound
      if (nclId && nclId !== "00000000") {
        setStatus("fetching_smiles");
        const smiles = smilesMap[parseInt(nclId)];
        if (smiles) {
          setStatus("setting_params");
          await apiPost(`jobs/${newJobId}/set_parameter/`, {
            object_path: "SubstituteLigand.inputData.SMILESIN",
            value: smiles,
          });
        }
      } else {
        // No ligand - set as NONE
        setStatus("setting_params");
        await apiPost(`jobs/${newJobId}/set_parameter/`, {
          object_path: "SubstituteLigand.controlParameters.LIGANDAS",
          value: "NONE",
        });
      }

      // 5. Set pipeline to DIMPLE
      await apiPost(`jobs/${newJobId}/set_parameter/`, {
        object_path: "SubstituteLigand.inputData.PIPELINE",
        value: "DIMPLE",
      });

      // 6. Upload reflection file using upload_file_param endpoint
      setStatus("uploading_reflections");
      const reflFormData = new FormData();
      reflFormData.append("file", file, file.name);
      reflFormData.append("object_path", "SubstituteLigand.inputData.UNMERGEDFILES[0].file");

      // Unmerged reflection files are the big ones - tens to hundreds of
      // megabytes - so this upload reports progress rather than sitting
      // silent for minutes.
      await apiUploadWithProgress(
        `jobs/${newJobId}/upload_file_param/`,
        reflFormData,
        {
          onProgress: ({ loaded, total }) =>
            dispatch({ type: "UPDATE_PROGRESS", file, loaded, total }),
        }
      );

      // 7. Queue job (in Azure mode, this queues via Service Bus)
      setStatus("queuing");
      await apiPost(`jobs/${newJobId}/run/`, {});

      setStatus("done");
    } catch (error) {
      dispatch({
        type: "UPDATE_STATUS",
        file,
        status: "error",
        error: describeFailure(step, error),
      });
      throw error;
    }
  };

  // Process all files sequentially
  const processAllFiles = async () => {
    setProcessing(true);
    let successCount = 0;
    let errorCount = 0;

    for (let i = 0; i < files.length; i++) {
      if (files[i].status === "idle" || files[i].status === "error") {
        setCurrentFileIndex(i);
        try {
          await processFile(files[i]);
          successCount++;
        } catch (error) {
          // Continue with next file even if one fails
          console.error("Error processing file:", error);
          errorCount++;
        }
      }
    }

    setCurrentFileIndex(-1);
    setProcessing(false);

    // Failures used to be reported only by a small chip in the table,
    // with the reason hidden in a tooltip - easy to read as "nothing
    // happened". Say so out loud.
    if (errorCount > 0) {
      setMessage(
        errorCount === 1
          ? "1 dataset failed to import - see the table for the reason"
          : `${errorCount} datasets failed to import - see the table for the reasons`,
        "error"
      );
    }

    // Auto-close on success and show popcorn message
    if (errorCount === 0 && successCount > 0) {
      const message = successCount === 1
        ? "1 dataset queued for processing"
        : `${successCount} datasets queued for processing`;
      setMessage(message, "success");
      dispatch({ type: "CLEAR" });
      onSuccess?.();
      onClose();
    }
  };

  const handleClose = () => {
    if (!processing) {
      dispatch({ type: "CLEAR" });
      onClose();
    }
  };

  const allDone = files.length > 0 && files.every((f) => f.status === "done");
  const failedFiles = files.filter((f) => f.status === "error");
  const hasErrors = failedFiles.length > 0;

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="lg" fullWidth>
      <DialogTitle>Batch Import Datasets</DialogTitle>
      <DialogContent>
        <Stack spacing={2}>
          <Alert severity="info">
            Drop unmerged reflection files (MTZ/SCA/CIF) to create member
            projects. Files should follow the naming convention:
            <br />
            <code>visit_crystal_NCL-XXXXXXXX_processing.mtz</code>
          </Alert>

          {/* Drop zone */}
          <DropZone
            onFilesSelected={handleFilesSelected}
            accept=".mtz,.sca,.cif"
            multiple
            sx={{ p: 4 }}
          >
            <UploadIcon sx={{ fontSize: 48, color: "text.secondary", mb: 1 }} />
            <Typography>
              Drop files here or click to select
            </Typography>
            <Typography variant="caption" color="text.secondary">
              MTZ, SCA, or CIF files
            </Typography>
          </DropZone>

          {/* Files table */}
          {files.length > 0 && (
            <TableContainer component={Paper}>
              <Table size="small">
                <TableHead>
                  <TableRow>
                    <TableCell>Visit</TableCell>
                    <TableCell>Crystal</TableCell>
                    <TableCell>Compound</TableCell>
                    <TableCell>Processing</TableCell>
                    <TableCell align="right">Size</TableCell>
                    <TableCell sx={{ minWidth: 260 }}>Status</TableCell>
                  </TableRow>
                </TableHead>
                <TableBody>
                  {files.map((item, index) => (
                    <TableRow
                      key={index}
                      sx={{
                        bgcolor:
                          currentFileIndex === index
                            ? "action.selected"
                            : undefined,
                      }}
                    >
                      <TableCell>{item.visit || "-"}</TableCell>
                      <TableCell>{item.crystal || "-"}</TableCell>
                      <TableCell>
                        {item.nclId === "00000000" ? (
                          <Chip label="Apo" size="small" variant="outlined" />
                        ) : item.nclId ? (
                          <Tooltip
                            title={
                              smilesMap[parseInt(item.nclId)] || "Loading..."
                            }
                          >
                            <Chip
                              label={`NCL-${item.nclId}`}
                              size="small"
                              color="primary"
                              variant="outlined"
                            />
                          </Tooltip>
                        ) : (
                          <Chip label="Unknown" size="small" color="warning" />
                        )}
                      </TableCell>
                      <TableCell>{item.processing || "-"}</TableCell>
                      <TableCell align="right">
                        {formatBytes(item.file.size)}
                      </TableCell>
                      <TableCell sx={{ minWidth: 260 }}>
                        <StatusChip status={item.status} />
                        {item.progress && (
                          <Box sx={{ mt: 0.5 }}>
                            <LinearProgress
                              variant={
                                item.progress.total > 0
                                  ? "determinate"
                                  : "indeterminate"
                              }
                              value={
                                item.progress.total > 0
                                  ? (item.progress.loaded / item.progress.total) *
                                    100
                                  : undefined
                              }
                            />
                            <Typography variant="caption" color="text.secondary">
                              {formatBytes(item.progress.loaded)}
                              {item.progress.total > 0
                                ? ` of ${formatBytes(item.progress.total)}`
                                : ""}{" "}
                              sent
                            </Typography>
                          </Box>
                        )}
                        {item.status === "error" && item.error && (
                          <Typography
                            variant="caption"
                            color="error"
                            sx={{
                              display: "block",
                              mt: 0.5,
                              whiteSpace: "normal",
                            }}
                          >
                            {item.error}
                          </Typography>
                        )}
                      </TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </TableContainer>
          )}

          {processing && (
            <Box>
              <Typography variant="body2" sx={{ mb: 1 }}>
                Processing file {currentFileIndex + 1} of {files.length}...
              </Typography>
              <LinearProgress
                variant="determinate"
                value={((currentFileIndex + 1) / files.length) * 100}
              />
            </Box>
          )}

          {hasErrors && !processing && (
            <Alert severity="error">
              <AlertTitle>
                {failedFiles.length} of {files.length}{" "}
                {files.length === 1 ? "dataset" : "datasets"} failed
              </AlertTitle>
              <Box component="ul" sx={{ pl: 2, m: 0 }}>
                {failedFiles.map((item, index) => (
                  <li key={index}>
                    <Typography variant="body2" component="span">
                      <strong>{item.file.name}</strong> - {item.error}
                    </Typography>
                  </li>
                ))}
              </Box>
              <Typography variant="body2" sx={{ mt: 1 }}>
                Click Process again to retry the failed datasets. Projects and
                jobs already created are reused, so retrying will not add
                duplicate members to the campaign.
              </Typography>
            </Alert>
          )}

          {allDone && (
            <Alert severity="success">
              All files processed successfully!
            </Alert>
          )}
        </Stack>
      </DialogContent>
      <DialogActions>
        <Button
          onClick={() => dispatch({ type: "CLEAR" })}
          disabled={processing || files.length === 0}
          startIcon={<ClearIcon />}
        >
          Clear
        </Button>
        <Button onClick={handleClose} disabled={processing}>
          {allDone ? "Done" : "Cancel"}
        </Button>
        <Button
          onClick={processAllFiles}
          variant="contained"
          disabled={processing || files.length === 0 || allDone}
          startIcon={<RunIcon />}
        >
          Process
        </Button>
      </DialogActions>
    </Dialog>
  );
}

// Helper component for status display
function StatusChip({ status }: { status: BatchFileStatus }) {
  const { label, color } = STATUS_META[status];
  return <Chip label={label} size="small" color={color} />;
}
