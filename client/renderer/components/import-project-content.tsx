"use client";
import React, { useCallback, useState } from "react";
import {
  Alert,
  Container,
  LinearProgress,
  Paper,
  Stack,
  Typography,
} from "@mui/material";
import { Archive } from "@mui/icons-material";
import { useApi } from "../api";
import { apiUploadWithProgress, UploadProgress } from "../api-fetch";
import { stagingCapability, stageFile } from "../lib/staged-upload";
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

  const handleFileUpload = useCallback(
    async (selectedFiles: FileList | null) => {
      if (!selectedFiles || selectedFiles.length === 0) return;

      // Desktop: if every dropped file has a real local path, import by path
      // (server-side copy) instead of pushing multi-GB zips through the ingress
      // cap -- mirrors utils.ts uploadFileParam, and the server only honours
      // local_path when its gate allows it (desktop, or a staged cloud dir).
      // getPathForFile lives only in the Electron preload; in the browser it is
      // absent, so localPaths stays empty and we upload the bytes as before.
      const files = Array.from(selectedFiles);
      const localPaths = files
        .map((f) => window.electronAPI?.getPathForFile?.(f) || "")
        .filter(Boolean);

      setUploading(true);
      setProgress(null);
      setError(null);
      try {
        const formData = new FormData();
        if (localPaths.length === files.length) {
          // Desktop: import each zip by its local path (server-side copy).
          for (const p of localPaths) formData.append("local_path", p);
        } else {
          // Served deployment: if staging is advertised, deliver each zip in
          // chunks past the ingress/body caps and import by owner-bound handles
          // (import_project takes all-staged or all-body, so stage every file);
          // otherwise upload the bytes as before.
          const cap = await stagingCapability();
          if (cap) {
            const sizes = files.map((f) => f.size);
            const totalBytes = sizes.reduce((a, b) => a + b, 0) || 1;
            let sentBefore = 0;
            for (let i = 0; i < files.length; i++) {
              const handle = await stageFile(files[i], files[i].name, cap, {
                onProgress: (frac) => {
                  const loaded = sentBefore + frac * sizes[i];
                  setProgress({ loaded, total: totalBytes, fraction: loaded / totalBytes });
                },
              });
              sentBefore += sizes[i];
              formData.append("staged_upload", handle);
            }
          } else {
            for (const f of files) formData.append("files", f);
          }
        }

        // apiUploadWithProgress, not api.post: the ordinary JSON path puts a
        // 30 s AbortController around the request, and a project zip is
        // routinely far bigger than 30 s of uplink. This one bounds on a
        // stall instead, and can say how far it got. (With staged handles the
        // body is tiny; the progress above came from staging.)
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
            <DropZone
              onFilesSelected={(files) => void handleFileUpload(files)}
              accept=".zip"
              multiple
              disabled={uploading}
              sx={{ p: 4 }}
            >
              <Typography variant="body1" color="textSecondary">
                Drag and drop files here, or click here to upload
              </Typography>
            </DropZone>

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
