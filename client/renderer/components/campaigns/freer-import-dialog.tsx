"use client";

import { useState, useCallback } from "react";
import {
  Box,
  Button,
  CircularProgress,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  LinearProgress,
  Paper,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import { CloudUpload as UploadIcon } from "@mui/icons-material";
import { apiPost, apiUpload } from "../../api-fetch";

interface FreeRImportDialogProps {
  open: boolean;
  onClose: () => void;
  parentProjectId: number;
  parentProjectName: string;
}

export function FreeRImportDialog({
  open,
  onClose,
  parentProjectId,
  parentProjectName,
}: FreeRImportDialogProps) {
  const [fSigFFile, setFSigFFile] = useState<File | null>(null);
  const [freeRFile, setFreeRFile] = useState<File | null>(null);
  const [maxResolution, setMaxResolution] = useState<string>("");
  const [progress, setProgress] = useState(0);
  const [progressMessage, setProgressMessage] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  const handleImport = useCallback(async () => {
    if (!fSigFFile) return;

    setLoading(true);
    setProgress(10);
    setProgressMessage("Creating FreeR import job...");

    try {
      // 1. Create freerflag job
      // The create_task endpoint uses api_success() which wraps the response
      const jobResponse = await apiPost<{
        success: boolean;
        data: { new_job: { id: number; uuid: string } };
      }>(
        `projects/${parentProjectId}/create_task/`,
        {
          task_name: "freerflag",
          title: `Generate FreeR from ${fSigFFile.name}`,
        }
      );
      const newJobId = jobResponse.data.new_job.id;

      // 2. Set parameters
      setProgress(20);
      setProgressMessage("Setting parameters...");

      if (maxResolution && parseFloat(maxResolution) > 0.5) {
        await apiPost(`jobs/${newJobId}/set_parameter/`, {
          object_path: "freerflag.controlParameters.RESMAX",
          value: maxResolution,
        });
      }

      // Set generation mode based on whether we have existing FreeR
      const genMode = freeRFile ? "COMPLETE" : "NEW";
      await apiPost(`jobs/${newJobId}/set_parameter/`, {
        object_path: "freerflag.controlParameters.GEN_MODE",
        value: genMode,
      });

      // 3. Upload F/SigF file using upload_file_param endpoint
      setProgress(30);
      setProgressMessage("Uploading observed data...");

      const fSigFFormData = new FormData();
      fSigFFormData.append("file", fSigFFile, fSigFFile.name);
      fSigFFormData.append("objectPath", "freerflag.inputData.F_SIGF");

      await apiUpload(`jobs/${newJobId}/upload_file_param/`, fSigFFormData);

      // 4. Upload FreeR file if provided
      if (freeRFile) {
        setProgress(50);
        setProgressMessage("Uploading existing FreeR set...");

        const freeRFormData = new FormData();
        freeRFormData.append("file", freeRFile, freeRFile.name);
        freeRFormData.append("objectPath", "freerflag.inputData.FREERFLAG");

        await apiUpload(`jobs/${newJobId}/upload_file_param/`, freeRFormData);
      }

      // 5. Run the job synchronously (freerflag is quick)
      setProgress(70);
      setProgressMessage("Running FreeR generation...");

      // Use run_local with synchronous=true to wait for completion
      await apiPost(`jobs/${newJobId}/run_local/`, { synchronous: true });

      setProgress(100);
      setProgressMessage("Done!");

      // Close after a moment
      setTimeout(() => {
        handleClose();
      }, 1000);
    } catch (error) {
      console.error("Failed to import FreeR:", error);
      setProgressMessage(`Error: ${error instanceof Error ? error.message : "Unknown error"}`);
    } finally {
      setLoading(false);
    }
  }, [fSigFFile, freeRFile, maxResolution, parentProjectId]);

  const handleClose = () => {
    if (!loading) {
      setFSigFFile(null);
      setFreeRFile(null);
      setMaxResolution("");
      setProgress(0);
      setProgressMessage(null);
      onClose();
    }
  };

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="sm" fullWidth>
      <DialogTitle>Import Collective FreeR Set</DialogTitle>
      <DialogContent>
        <Stack spacing={3}>
          {/* Progress indicator */}
          {loading && (
            <Box sx={{ display: "flex", alignItems: "center", gap: 2 }}>
              <CircularProgress variant="determinate" value={progress} size={40} />
              <Typography>{progressMessage}</Typography>
            </Box>
          )}

          {/* F/SigF drop zone */}
          <Box>
            <Typography variant="subtitle2" gutterBottom>
              Exemplar observed data (F/SigF) *
            </Typography>
            <DropZone
              file={fSigFFile}
              onFileSelect={setFSigFFile}
              accept=".mtz"
              disabled={loading}
              placeholder="Drop MTZ file with F/SigF columns"
            />
          </Box>

          {/* FreeR drop zone (optional) */}
          <Box>
            <Typography variant="subtitle2" gutterBottom>
              Starting FreeR set (optional)
            </Typography>
            <DropZone
              file={freeRFile}
              onFileSelect={setFreeRFile}
              accept=".mtz"
              disabled={loading}
              placeholder="Drop MTZ file with existing FreeR flags"
            />
            <Typography variant="caption" color="text.secondary">
              If not provided, a new FreeR set will be generated
            </Typography>
          </Box>

          {/* Max resolution */}
          <TextField
            label="Extend to maximum resolution"
            type="number"
            value={maxResolution}
            onChange={(e) => setMaxResolution(e.target.value)}
            disabled={loading}
            size="small"
            inputProps={{ step: 0.1, min: 0.5 }}
            helperText="Leave empty to use resolution from input file"
          />

          {loading && <LinearProgress variant="determinate" value={progress} />}
        </Stack>
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} disabled={loading}>
          Cancel
        </Button>
        <Button
          onClick={handleImport}
          variant="contained"
          disabled={!fSigFFile || loading}
        >
          Import
        </Button>
      </DialogActions>
    </Dialog>
  );
}

// Reusable drop zone component
function DropZone({
  file,
  onFileSelect,
  accept,
  disabled,
  placeholder,
}: {
  file: File | null;
  onFileSelect: (file: File | null) => void;
  accept: string;
  disabled: boolean;
  placeholder: string;
}) {
  return (
    <Paper
      sx={{
        p: 2,
        border: "2px dashed",
        borderColor: file ? "success.main" : "divider",
        textAlign: "center",
        cursor: disabled ? "default" : "pointer",
        opacity: disabled ? 0.5 : 1,
        "&:hover": disabled ? {} : { borderColor: "primary.main" },
      }}
      onClick={() => {
        if (disabled) return;
        const input = document.createElement("input");
        input.type = "file";
        input.accept = accept;
        input.onchange = (e) => {
          const files = (e.target as HTMLInputElement).files;
          if (files && files.length > 0) {
            onFileSelect(files[0]);
          }
        };
        input.click();
      }}
      onDragOver={(e) => {
        e.preventDefault();
        if (!disabled) e.currentTarget.style.borderColor = "#1976d2";
      }}
      onDragLeave={(e) => {
        e.currentTarget.style.borderColor = file ? "#2e7d32" : "";
      }}
      onDrop={(e) => {
        e.preventDefault();
        if (!disabled && e.dataTransfer.files.length > 0) {
          onFileSelect(e.dataTransfer.files[0]);
        }
      }}
    >
      {file ? (
        <Typography color="success.main">{file.name}</Typography>
      ) : (
        <Typography color="text.secondary">{placeholder}</Typography>
      )}
    </Paper>
  );
}
