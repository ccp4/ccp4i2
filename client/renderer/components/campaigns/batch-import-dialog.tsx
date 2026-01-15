"use client";

import { useCallback, useReducer, useState } from "react";
import {
  Alert,
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
import { apiPost, apiFetch, apiUpload } from "../../api-fetch";
import {
  BatchFileItem,
  BatchFileStatus,
  parseDatasetFilename,
} from "../../types/campaigns";

interface BatchImportDialogProps {
  open: boolean;
  onClose: () => void;
  campaignId: number;
  parentProjectId: number;
  latestCoordsFileId: number;
}

// Reducer for managing file state
type FileAction =
  | { type: "ADD"; file: File }
  | { type: "UPDATE_STATUS"; file: File; status: BatchFileStatus; error?: string }
  | { type: "UPDATE_SMILES"; nclId: string; smiles: string }
  | { type: "SET_PROJECT_ID"; file: File; projectId: number; jobId: number }
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
      return state.map((item) =>
        item.file === action.file
          ? { ...item, status: action.status, error: action.error }
          : item
      );
    }
    case "UPDATE_SMILES": {
      return state.map((item) =>
        item.nclId === action.nclId ? { ...item, smiles: action.smiles } : item
      );
    }
    case "SET_PROJECT_ID": {
      return state.map((item) =>
        item.file === action.file
          ? { ...item, projectId: action.projectId, jobId: action.jobId }
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
  campaignId,
  parentProjectId,
  latestCoordsFileId,
}: BatchImportDialogProps) {
  const [files, dispatch] = useReducer(filesReducer, []);
  const [processing, setProcessing] = useState(false);
  const [currentFileIndex, setCurrentFileIndex] = useState(-1);

  const campaignsApi = useCampaignsApi();

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

    try {
      // 1. Create new project
      dispatch({ type: "UPDATE_STATUS", file, status: "creating_project" });
      const projectName = file.name.replace(/\.[^.]+$/, "").replace(/\s/g, "_");
      const projectResponse = await apiPost<{ id: number }>("projects/", {
        name: projectName,
        directory: "__default__",
      });
      const newProjectId = projectResponse.id;

      // Add project to campaign as member
      await campaignsApi.addMember(campaignId, newProjectId, "member");

      // 2. Create SubstituteLigand job
      // The create_task endpoint uses api_success() which wraps the response
      dispatch({ type: "UPDATE_STATUS", file, status: "creating_job" });
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
      const newJobId = jobResponse.data.new_job.id;

      dispatch({
        type: "SET_PROJECT_ID",
        file,
        projectId: newProjectId,
        jobId: newJobId,
      });

      // 3. Upload reference coordinates using upload_file_param endpoint
      dispatch({ type: "UPDATE_STATUS", file, status: "uploading_coords" });
      // Get coords file content (use download/ with numeric id, not download_by_uuid/)
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
      coordsFormData.append("objectPath", "SubstituteLigand.inputData.XYZIN");

      // Use apiUpload for authenticated upload
      await apiUpload(`jobs/${newJobId}/upload_file_param/`, coordsFormData);

      // 4. Set SMILES if we have a compound
      if (nclId && nclId !== "00000000") {
        dispatch({ type: "UPDATE_STATUS", file, status: "fetching_smiles" });
        const smiles = smilesMap[parseInt(nclId)];
        if (smiles) {
          dispatch({ type: "UPDATE_STATUS", file, status: "setting_params" });
          await apiPost(`jobs/${newJobId}/set_parameter/`, {
            object_path: "SubstituteLigand.inputData.SMILESIN",
            value: smiles,
          });
        }
      } else {
        // No ligand - set as NONE
        dispatch({ type: "UPDATE_STATUS", file, status: "setting_params" });
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
      dispatch({ type: "UPDATE_STATUS", file, status: "uploading_reflections" });
      const reflFormData = new FormData();
      reflFormData.append("file", file, file.name);
      reflFormData.append("objectPath", "SubstituteLigand.inputData.UNMERGEDFILES[0].file");

      // Use apiUpload for authenticated upload
      await apiUpload(`jobs/${newJobId}/upload_file_param/`, reflFormData);

      // 7. Queue job (in Azure mode, this queues via Service Bus)
      dispatch({ type: "UPDATE_STATUS", file, status: "queuing" });
      await apiPost(`jobs/${newJobId}/run/`, {});

      dispatch({ type: "UPDATE_STATUS", file, status: "done" });
    } catch (error) {
      dispatch({
        type: "UPDATE_STATUS",
        file,
        status: "error",
        error: error instanceof Error ? error.message : "Unknown error",
      });
      throw error;
    }
  };

  // Process all files sequentially
  const processAllFiles = async () => {
    setProcessing(true);

    for (let i = 0; i < files.length; i++) {
      if (files[i].status === "idle" || files[i].status === "error") {
        setCurrentFileIndex(i);
        try {
          await processFile(files[i]);
        } catch (error) {
          // Continue with next file even if one fails
          console.error("Error processing file:", error);
        }
      }
    }

    setCurrentFileIndex(-1);
    setProcessing(false);
  };

  const handleClose = () => {
    if (!processing) {
      dispatch({ type: "CLEAR" });
      onClose();
    }
  };

  const allDone = files.length > 0 && files.every((f) => f.status === "done");
  const hasErrors = files.some((f) => f.status === "error");

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
          <Paper
            sx={{
              p: 4,
              border: "2px dashed",
              borderColor: "divider",
              textAlign: "center",
              cursor: "pointer",
              "&:hover": { borderColor: "primary.main" },
            }}
            onClick={() => {
              const input = document.createElement("input");
              input.type = "file";
              input.multiple = true;
              input.accept = ".mtz,.sca,.cif";
              input.onchange = (e) => {
                handleFilesSelected((e.target as HTMLInputElement).files);
              };
              input.click();
            }}
            onDragOver={(e) => e.preventDefault()}
            onDrop={(e) => {
              e.preventDefault();
              handleFilesSelected(e.dataTransfer.files);
            }}
          >
            <UploadIcon sx={{ fontSize: 48, color: "text.secondary", mb: 1 }} />
            <Typography>
              Drop files here or click to select
            </Typography>
            <Typography variant="caption" color="text.secondary">
              MTZ, SCA, or CIF files
            </Typography>
          </Paper>

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
                    <TableCell>Status</TableCell>
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
                      <TableCell>
                        <StatusChip status={item.status} error={item.error} />
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
            <Alert severity="warning">
              Some files failed to process. You can retry by clicking Process
              again.
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
function StatusChip({
  status,
  error,
}: {
  status: BatchFileStatus;
  error?: string;
}) {
  const statusMap: Record<BatchFileStatus, { label: string; color: "default" | "primary" | "success" | "error" | "warning" | "info" }> = {
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

  const { label, color } = statusMap[status];

  return (
    <Tooltip title={error || ""}>
      <Chip label={label} size="small" color={color} />
    </Tooltip>
  );
}
