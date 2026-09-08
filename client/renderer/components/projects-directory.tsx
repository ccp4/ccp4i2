"use client";
import {
  Button,
  CircularProgress,
  Stack,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import { Folder } from "@mui/icons-material";
import React, { useCallback, useEffect, useState } from "react";
import { apiGet, apiPatch } from "../api-fetch";

/** Open the native directory picker via Electron IPC; null in the web build. */
async function browseDirectory(title: string): Promise<string | null> {
  const api = typeof window !== "undefined" ? window.electronAPI : undefined;
  if (!api?.invoke) return null;
  try {
    return (await api.invoke("browse-path", { mode: "directory", title })) ?? null;
  } catch {
    return null;
  }
}

/**
 * The default projects directory: where a new project lands unless the user
 * picks somewhere else for that one project. Reads/writes preferences.json
 * via the config API — same desktop/cloud split as Program locations: a
 * cloud deployment reports editable=false (it configures this via the
 * CCP4I2_PROJECTS_DIR environment variable instead) and the panel is
 * read-only there.
 */
export function ProjectsDirectory() {
  const [editable, setEditable] = useState(false);
  const [directory, setDirectory] = useState<string>("");
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);

  const load = useCallback(async () => {
    const resp = await apiGet<any>("config/default-project-parent/");
    const data = resp?.data ?? resp;
    setEditable(Boolean(data?.editable));
    setDirectory(data?.directory ?? "");
  }, []);

  useEffect(() => {
    (async () => {
      setLoading(true);
      try {
        await load();
      } finally {
        setLoading(false);
      }
    })();
  }, [load]);

  const persist = useCallback(
    async (newDirectory: string | null) => {
      setSaving(true);
      try {
        await apiPatch("config/default-project-parent/set/", {
          directory: newDirectory,
        });
        await load();
      } finally {
        setSaving(false);
      }
    },
    [load]
  );

  const handleChange = async () => {
    const picked = await browseDirectory("Select the default projects directory");
    if (picked && picked !== directory) persist(picked);
  };
  const handleReset = () => persist(null);

  return (
    <Stack spacing={2} sx={{ maxWidth: 720, mx: "auto", p: 2 }}>
      <Typography variant="h6">Projects directory</Typography>
      {loading ? (
        <CircularProgress size={20} />
      ) : (
        <Stack direction="row" spacing={2} alignItems="center">
          <TextField
            label="Where a new project goes by default"
            value={directory}
            disabled
            fullWidth
            helperText={
              editable
                ? undefined
                : "Set via the CCP4I2_PROJECTS_DIR environment variable in this deployment."
            }
          />
          {editable && (
            <>
              <Tooltip title="Choose a different default projects directory">
                <Button
                  variant="outlined"
                  startIcon={<Folder />}
                  onClick={handleChange}
                  disabled={saving}
                  sx={{ flexShrink: 0 }}
                >
                  Change
                </Button>
              </Tooltip>
              <Button variant="text" onClick={handleReset} disabled={saving}>
                Reset to default
              </Button>
            </>
          )}
        </Stack>
      )}
    </Stack>
  );
}
