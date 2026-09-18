"use client";
import {
  Alert,
  Button,
  CircularProgress,
  Stack,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import { FolderOpen } from "@mui/icons-material";
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

interface ProjectsDirectoryState {
  directory: string;
  default: string;
  editable: boolean;
}

/**
 * The default projects directory: where a new project lands unless the user
 * picks somewhere else for that one project. Reads/writes preferences.json
 * via the config API — same desktop/cloud split as Program locations: a
 * cloud deployment reports editable=false (it configures this via the
 * CCP4I2_PROJECTS_DIR environment variable instead) and the panel is
 * read-only there.
 *
 * "Reset to default" is the point of the panel. Before it, the only way back
 * to <ccp4i2 home>/projects was to remember the path and navigate to a hidden
 * directory.
 */
export function ProjectsDirectory() {
  const [state, setState] = useState<ProjectsDirectoryState>({
    directory: "",
    default: "",
    editable: false,
  });
  const [draft, setDraft] = useState("");
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const load = useCallback(async () => {
    const resp = await apiGet<any>("config/default-project-parent/");
    const data = resp?.data ?? resp;
    setState({
      directory: data?.directory ?? "",
      default: data?.default ?? "",
      editable: Boolean(data?.editable),
    });
    setDraft(data?.directory ?? "");
  }, []);

  useEffect(() => {
    (async () => {
      setLoading(true);
      try {
        await load();
      } catch (err) {
        setError(`Could not read the projects directory: ${err}`);
      } finally {
        setLoading(false);
      }
    })();
  }, [load]);

  const persist = useCallback(
    async (directory: string | null) => {
      setSaving(true);
      setError(null);
      try {
        // Reload rather than trusting the request body: the server creates
        // the directory and resolves a reset, so what it reports back is the
        // only account of what actually happened.
        await apiPatch("config/default-project-parent/set/", { directory });
        await load();
      } catch (err) {
        // The server sends a sentence explaining the refusal; the wrapper
        // puts it in the Error's message.
        setError(err instanceof Error ? err.message : `${err}`);
        await load().catch(() => undefined);
      } finally {
        setSaving(false);
      }
    },
    [load]
  );

  const handleBrowse = async () => {
    const picked = await browseDirectory("Select the default projects directory");
    if (picked && picked !== state.directory) persist(picked);
  };

  const handleBlur = () => {
    const typed = draft.trim();
    if (typed && typed !== state.directory) persist(typed);
    else setDraft(state.directory);
  };

  const isDefault = Boolean(state.default) && state.directory === state.default;

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
            value={draft}
            onChange={(event) => setDraft(event.target.value)}
            onBlur={handleBlur}
            disabled={!state.editable || saving}
            fullWidth
            size="small"
            helperText={
              !state.editable
                ? "Set with the CCP4I2_PROJECTS_DIR environment variable in this deployment."
                : isDefault
                  ? "This is the default."
                  : `Reset restores ${state.default}`
            }
          />
          {state.editable && (
            <>
              <Tooltip title="Choose a different default projects directory">
                <span>
                  <Button
                    variant="outlined"
                    startIcon={<FolderOpen />}
                    onClick={handleBrowse}
                    disabled={saving}
                    sx={{ flexShrink: 0 }}
                  >
                    Change
                  </Button>
                </span>
              </Tooltip>
              <Tooltip
                title={
                  isDefault
                    ? "Already the default"
                    : `Restore ${state.default}`
                }
              >
                <span>
                  <Button
                    variant="text"
                    onClick={() => persist(null)}
                    disabled={saving || isDefault}
                    sx={{ flexShrink: 0 }}
                  >
                    Reset
                  </Button>
                </span>
              </Tooltip>
            </>
          )}
        </Stack>
      )}
    </Stack>
  );
}
