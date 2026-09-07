"use client";
import React, { useCallback, useRef, useState } from "react";
import {
  Alert,
  Button,
  Container,
  LinearProgress,
  Paper,
  Stack,
  Typography,
} from "@mui/material";
import { Archive, Upload } from "@mui/icons-material";
import { useApi } from "../api";
import { apiUploadWithProgress, UploadProgress } from "../api-fetch";
import { VisuallyHiddenInput } from "./task/task-elements/input-file-upload";
import { useRouter } from "next/navigation";
import { Project } from "../types/models";
import { ImportProjectDirectory } from "./import-project-directory";
import { DropZone } from "./common/drop-zone";

export const ImportProjectContent: React.FC = () => {
  const api = useApi();
  const router = useRouter();
  const [uploading, setUploading] = useState(false);
  const [progress, setProgress] = useState<UploadProgress | null>(null);
  const [error, setError] = useState<string | null>(null);
  const { mutate: mutateProjects } = api.get<Project[]>("projects");

  // Create a ref for the hidden file input
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleFileUpload = useCallback(
    async (selectedFiles: FileList | null) => {
      if (!selectedFiles || selectedFiles.length === 0) return;

      const formData = new FormData();
      for (let i = 0; i < selectedFiles.length; i++) {
        formData.append("files", selectedFiles[i]);
      }

      setUploading(true);
      setProgress(null);
      setError(null);
      try {
        // apiUploadWithProgress, not api.post: the ordinary JSON path puts a
        // 30 s AbortController around the request, and a project zip is
        // routinely far bigger than 30 s of uplink. This one bounds on a
        // stall instead, and can say how far it got.
        const response: any = await apiUploadWithProgress(
          "projects/import_project/",
          formData,
          { onProgress: setProgress }
        );
        if (response?.success === false) {
          setError(response?.error ?? "The project could not be imported");
          return;
        }
        mutateProjects();
        router.push("/ccp4i2");
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err));
      } finally {
        setUploading(false);
        setProgress(null);
      }
    },
    [mutateProjects, router]
  );

  const onChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    void handleFileUpload(event.target.files);
  };

  return (
    <Container
      sx={{
        display: "flex",
        justifyContent: "center",
        // The page supplies a scrolling flex content area beneath the top
        // bar; auto margins centre this card in it and shrink to nothing
        // when the card is taller than the area, so nothing is clipped.
        margin: "auto",
        paddingY: 4,
      }}
    >
      <Stack spacing={3} sx={{ padding: 2, minWidth: "50rem" }}>
        <Stack spacing={0.5}>
          <Typography variant="h4">Import Project(s)</Typography>
          <Typography variant="body2" color="text.secondary">
            Bring in a project from elsewhere — as a zip, or as a folder that
            is already on this machine.
          </Typography>
        </Stack>

        {error && <Alert severity="error">{error}</Alert>}

        <Paper variant="outlined" sx={{ padding: 2 }}>
          <Stack spacing={2}>
            <Stack direction="row" spacing={1} alignItems="center">
              <Archive color="primary" />
              <Typography variant="h6">A project zip</Typography>
            </Stack>
            <Typography variant="body2" color="text.secondary">
              An exported <code>.ccp4_project.zip</code>. Its contents are
              copied into your project store.
            </Typography>
            <Stack spacing={2} direction="row" alignItems="center">
              <DropZone
                onFilesSelected={(files) => void handleFileUpload(files)}
                accept=".zip"
                multiple
                disabled={uploading}
                sx={{ p: 4, flexGrow: 1 }}
              >
                <Typography variant="body1" color="textSecondary">
                  Drag and drop files here, or click here to upload
                </Typography>
              </DropZone>
              <Button
                component="label"
                variant="contained"
                startIcon={<Upload />}
                disabled={uploading}
              >
                <VisuallyHiddenInput
                  ref={fileInputRef}
                  type="file"
                  multiple
                  accept=".zip"
                  onChange={onChange}
                />
              </Button>
            </Stack>

            {uploading && (
              <Stack spacing={0.5}>
                <LinearProgress
                  variant={
                    progress?.fraction != null ? "determinate" : "indeterminate"
                  }
                  value={
                    progress?.fraction != null ? progress.fraction * 100 : undefined
                  }
                />
                <Typography variant="caption" color="text.secondary">
                  {progress?.fraction != null
                    ? `Uploading — ${Math.round(progress.fraction * 100)}%`
                    : "Uploading…"}
                </Typography>
              </Stack>
            )}
          </Stack>
        </Paper>

        {/* Desktop only — renders nothing in a browser, where a folder the
            user picks is not on the server's disk. */}
        <ImportProjectDirectory onImported={() => mutateProjects()} />
      </Stack>
    </Container>
  );
};
