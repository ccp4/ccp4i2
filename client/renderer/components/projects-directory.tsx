"use client";
import {
  Alert,
  Button,
  CircularProgress,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import { FolderOpen } from "@mui/icons-material";
import React, { useCallback, useEffect, useState } from "react";
import {
  DefaultProjectsDir,
  getDefaultProjectsDir,
  setDefaultProjectsDir,
} from "../lib/default-projects-dir";
import { browsePath } from "../utils/browse-path";

/**
 * The default projects directory, and the way back to the built-in one.
 *
 * Read-only in a cloud deployment, which reports editable=false and takes
 * this from the CCP4I2_PROJECTS_DIR environment variable.
 */
export function ProjectsDirectory() {
  const [setting, setSetting] = useState<DefaultProjectsDir | null>(null);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const load = useCallback(
    () => getDefaultProjectsDir().then(setSetting),
    []
  );

  useEffect(() => {
    load().catch((err) => setError(`${err}`));
  }, [load]);

  // Reloaded rather than taken from the request: the server resolves a reset
  // and creates the directory, so only its account of it is reliable.
  const persist = async (directory: string | null) => {
    setSaving(true);
    setError(null);
    try {
      await setDefaultProjectsDir(directory);
    } catch (err) {
      setError(err instanceof Error ? err.message : `${err}`);
    }
    await load().catch(() => undefined);
    setSaving(false);
  };

  const handleChange = async () => {
    const picked = await browsePath({
      mode: "directory",
      title: "Select the default projects directory",
    });
    if (picked && picked !== setting?.directory) persist(picked);
  };

  const isDefault = setting?.directory === setting?.default;

  return (
    <Stack spacing={2} sx={{ maxWidth: 720, mx: "auto", p: 2 }}>
      <Typography variant="h6">Projects directory</Typography>
      {error && <Alert severity="error">{error}</Alert>}
      {!setting ? (
        <CircularProgress size={20} />
      ) : (
        <Stack direction="row" spacing={2} alignItems="flex-start">
          <TextField
            label="Where a new project goes by default"
            value={setting.directory}
            disabled
            fullWidth
            size="small"
            helperText={
              !setting.editable
                ? "Set with the CCP4I2_PROJECTS_DIR environment variable in this deployment."
                : isDefault
                  ? "This is the default."
                  : `Reset restores ${setting.default}`
            }
          />
          {setting.editable && (
            <>
              <Button
                variant="outlined"
                startIcon={<FolderOpen />}
                onClick={handleChange}
                disabled={saving}
                sx={{ flexShrink: 0 }}
              >
                Change
              </Button>
              <Button
                onClick={() => persist(null)}
                disabled={saving || isDefault}
                sx={{ flexShrink: 0 }}
              >
                Reset
              </Button>
            </>
          )}
        </Stack>
      )}
    </Stack>
  );
}
