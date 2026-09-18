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
import { apiGet, apiPatch } from "../api-fetch";
import { browsePath } from "../utils/browse-path";

const ENDPOINT = "config/default-project-parent/";

interface Setting {
  directory: string;
  default: string;
  editable: boolean;
}

const EMPTY: Setting = { directory: "", default: "", editable: false };

/**
 * The default projects directory, and the way back to the built-in one.
 *
 * Read-only in a cloud deployment, which reports editable=false and takes
 * this from the CCP4I2_PROJECTS_DIR environment variable.
 */
export function ProjectsDirectory() {
  const [setting, setSetting] = useState<Setting>(EMPTY);
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const load = useCallback(async () => {
    const resp = await apiGet<any>(ENDPOINT);
    setSetting({ ...EMPTY, ...(resp?.data ?? resp) });
  }, []);

  useEffect(() => {
    load()
      .catch((err) => setError(`${err}`))
      .finally(() => setLoading(false));
  }, [load]);

  // Reloaded rather than taken from the request: the server resolves a reset
  // and creates the directory, so only its account of it is reliable.
  const persist = async (directory: string | null) => {
    setSaving(true);
    setError(null);
    try {
      await apiPatch(ENDPOINT, { directory });
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
    if (picked && picked !== setting.directory) persist(picked);
  };

  const isDefault = Boolean(setting.default) && setting.directory === setting.default;

  return (
    <Stack spacing={2} sx={{ maxWidth: 720, mx: "auto", p: 2 }}>
      <Typography variant="h6">Projects directory</Typography>
      {error && <Alert severity="error">{error}</Alert>}
      {loading ? (
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
