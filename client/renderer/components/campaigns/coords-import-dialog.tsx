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
  Typography,
} from "@mui/material";
import { CloudUpload as UploadIcon } from "@mui/icons-material";
import { apiPost, apiUpload } from "../../api-fetch";

interface CoordsImportDialogProps {
  open: boolean;
  onClose: () => void;
  parentProjectId: number;
  parentProjectName: string;
}

export function CoordsImportDialog({
  open,
  onClose,
  parentProjectId,
  parentProjectName,
}: CoordsImportDialogProps) {
  const [coordFile, setCoordFile] = useState<File | null>(null);
  const [progress, setProgress] = useState(0);
  const [progressMessage, setProgressMessage] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  const handleFileSelect = useCallback((files: FileList | null) => {
    if (files && files.length > 0) {
      setCoordFile(files[0]);
    }
  }, []);

  const handleImport = useCallback(async () => {
    if (!coordFile) return;

    setLoading(true);
    setProgress(10);
    setProgressMessage("Creating coordinate import job...");

    try {
      // 1. Create coordinate_selector job
      // The create_task endpoint uses api_success() which wraps the response
      const jobResponse = await apiPost<{
        success: boolean;
        data: { new_job: { id: number; uuid: string } };
      }>(
        `projects/${parentProjectId}/create_task/`,
        {
          task_name: "coordinate_selector",
          title: `Import ${coordFile.name}`,
        }
      );
      const newJobId = jobResponse.data.new_job.id;

      setProgress(30);
      setProgressMessage("Uploading coordinate file...");

      // 2. Upload the coordinate file using upload_file_param endpoint
      // objectPath follows the pattern: taskname.inputData.PARAMNAME
      const formData = new FormData();
      formData.append("file", coordFile, coordFile.name);
      formData.append("objectPath", "coordinate_selector.inputData.XYZIN");

      await apiUpload(`jobs/${newJobId}/upload_file_param/`, formData);

      setProgress(60);
      setProgressMessage("Running coordinate import...");

      // 3. Run the job
      await apiPost(`jobs/${newJobId}/run/`, {});

      setProgress(90);
      setProgressMessage("Job queued, waiting for completion...");

      // 4. Wait a bit for the job to run (coordinate_selector is quick)
      await new Promise((resolve) => setTimeout(resolve, 5000));

      setProgress(100);
      setProgressMessage("Done!");

      // Close after a moment
      setTimeout(() => {
        handleClose();
      }, 1000);
    } catch (error) {
      console.error("Failed to import coordinates:", error);
      setProgressMessage(`Error: ${error instanceof Error ? error.message : "Unknown error"}`);
    } finally {
      setLoading(false);
    }
  }, [coordFile, parentProjectId]);

  const handleClose = () => {
    if (!loading) {
      setCoordFile(null);
      setProgress(0);
      setProgressMessage(null);
      onClose();
    }
  };

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="sm" fullWidth>
      <DialogTitle>Import Starting Coordinates</DialogTitle>
      <DialogContent>
        <Stack spacing={2}>
          {/* Progress indicator */}
          {loading && (
            <Box sx={{ display: "flex", alignItems: "center", gap: 2 }}>
              <CircularProgress variant="determinate" value={progress} size={40} />
              <Typography>{progressMessage}</Typography>
            </Box>
          )}

          {/* Drop zone */}
          <Paper
            sx={{
              p: 4,
              border: "2px dashed",
              borderColor: coordFile ? "success.main" : "divider",
              textAlign: "center",
              cursor: loading ? "default" : "pointer",
              opacity: loading ? 0.5 : 1,
              "&:hover": loading ? {} : { borderColor: "primary.main" },
            }}
            onClick={() => {
              if (loading) return;
              const input = document.createElement("input");
              input.type = "file";
              input.accept = ".pdb,.cif,.ent,.mmcif";
              input.onchange = (e) => {
                handleFileSelect((e.target as HTMLInputElement).files);
              };
              input.click();
            }}
            onDragOver={(e) => {
              e.preventDefault();
              if (!loading) e.currentTarget.style.borderColor = "#1976d2";
            }}
            onDragLeave={(e) => {
              e.currentTarget.style.borderColor = coordFile ? "#2e7d32" : "";
            }}
            onDrop={(e) => {
              e.preventDefault();
              if (!loading) handleFileSelect(e.dataTransfer.files);
            }}
          >
            <UploadIcon sx={{ fontSize: 48, color: "text.secondary", mb: 1 }} />
            {coordFile ? (
              <Typography color="success.main">{coordFile.name}</Typography>
            ) : (
              <>
                <Typography>Drop PDB/CIF file here or click to select</Typography>
                <Typography variant="caption" color="text.secondary">
                  Supported formats: .pdb, .cif, .ent, .mmcif
                </Typography>
              </>
            )}
          </Paper>

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
          disabled={!coordFile || loading}
        >
          Import
        </Button>
      </DialogActions>
    </Dialog>
  );
}
