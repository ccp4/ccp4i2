"use client";
import {
  Box,
  Button,
  Chip,
  CircularProgress,
  IconButton,
  Paper,
  Stack,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  CheckCircle,
  ErrorOutline,
  Add,
  Delete,
  Refresh,
  FolderOpen,
} from "@mui/icons-material";
import React, { useCallback, useEffect, useState } from "react";
import { apiGet, apiPatch } from "../api-fetch";

interface ProgramStatus {
  name: string;
  path: string | null;
  source: string;
  tasks?: string[];
}

// Explicit path fields shown as their own inputs (the rest go via exePaths).
// `browse` picks the native-picker mode: "file" for a single executable,
// "directory" for a suite install dir.
const EXPLICIT_FIELDS: {
  key: string;
  label: string;
  help: string;
  browse: "file" | "directory";
}[] = [
  { key: "COOT_EXECUTABLE", label: "Coot executable", help: "Full path to the coot binary", browse: "file" },
  { key: "CCP4MG_EXECUTABLE", label: "CCP4mg executable", help: "Full path to the ccp4mg binary", browse: "file" },
  { key: "SHELXDIR", label: "SHELX directory", help: "Directory containing shelxc/d/e/l", browse: "directory" },
  { key: "DIALSDIR", label: "DIALS directory", help: "Directory containing dials binaries", browse: "directory" },
  { key: "BUSTERDIR", label: "BUSTER directory", help: "Directory containing the BUSTER refine binary", browse: "directory" },
];

/** Open the native picker via Electron IPC; returns null in the web build. */
async function browsePath(
  mode: "file" | "directory",
  title?: string
): Promise<string | null> {
  const api = typeof window !== "undefined" ? window.electronAPI : undefined;
  if (!api?.invoke) return null;
  try {
    return (await api.invoke("browse-path", { mode, title })) ?? null;
  } catch {
    return null;
  }
}

/** True in the Electron desktop build (native picker available). */
function canBrowse(): boolean {
  return (
    typeof window !== "undefined" && Boolean(window.electronAPI?.invoke)
  );
}

const SOURCE_LABEL: Record<string, string> = {
  executable_pref: "explicit path",
  suite_dir: "suite directory",
  exe_paths: "extra search path",
  path: "PATH",
  missing: "not found",
};

/**
 * Program locations preferences: let a user point CCP4i2 at programs installed
 * off PATH. Reads/writes the desktop preferences.json via the config API; in a
 * cloud deployment the API reports editable=false and the panel is read-only
 * (values come from environment variables there).
 */
export function ProgramLocations() {
  const [editable, setEditable] = useState<boolean>(false);
  const [prefs, setPrefs] = useState<Record<string, any>>({});
  const [exePaths, setExePaths] = useState<string[]>([]);
  const [statuses, setStatuses] = useState<ProgramStatus[]>([]);
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [newPath, setNewPath] = useState("");

  const loadPrefs = useCallback(async () => {
    const resp = await apiGet<any>("config/program-preferences/");
    const data = resp?.data ?? resp;
    setEditable(Boolean(data?.editable));
    const p = data?.preferences ?? {};
    setPrefs(p);
    setExePaths(Array.isArray(p.exePaths) ? p.exePaths : []);
  }, []);

  const loadStatuses = useCallback(async () => {
    const resp = await apiGet<any>("config/discover-programs/");
    const data = resp?.data ?? resp;
    setStatuses(data?.programs ?? []);
  }, []);

  useEffect(() => {
    (async () => {
      setLoading(true);
      try {
        await Promise.all([loadPrefs(), loadStatuses()]);
      } finally {
        setLoading(false);
      }
    })();
  }, [loadPrefs, loadStatuses]);

  const persist = useCallback(
    async (updates: Record<string, any>) => {
      setSaving(true);
      try {
        await apiPatch("config/program-preferences/set/", updates);
        // Re-probe after a change so the status chips reflect it.
        await Promise.all([loadPrefs(), loadStatuses()]);
      } finally {
        setSaving(false);
      }
    },
    [loadPrefs, loadStatuses]
  );

  const handleExplicitBlur = (key: string, value: string) => {
    if ((prefs[key] ?? "") !== value) persist({ [key]: value });
  };
  const handleBrowseExplicit = async (
    key: string,
    mode: "file" | "directory",
    label: string
  ) => {
    const picked = await browsePath(mode, `Select ${label}`);
    if (picked && picked !== (prefs[key] ?? "")) {
      setPrefs((p) => ({ ...p, [key]: picked }));
      persist({ [key]: picked });
    }
  };
  const addPath = (p: string) => {
    const trimmed = p.trim();
    if (!trimmed || exePaths.includes(trimmed)) return;
    const next = [...exePaths, trimmed];
    setExePaths(next);
    setNewPath("");
    persist({ exePaths: next });
  };
  const handleAddPath = () => addPath(newPath);
  const handleBrowseAddPath = async () => {
    const picked = await browsePath("directory", "Add executable search directory");
    if (picked) addPath(picked);
  };
  const handleRemovePath = (p: string) => {
    const next = exePaths.filter((x) => x !== p);
    setExePaths(next);
    persist({ exePaths: next });
  };

  if (loading) {
    return (
      <Box sx={{ display: "flex", justifyContent: "center", py: 4 }}>
        <CircularProgress />
      </Box>
    );
  }

  return (
    <Stack spacing={3} sx={{ maxWidth: 720, mx: "auto", p: 2 }}>
      <Box>
        <Typography variant="h6">Program locations</Typography>
        <Typography variant="body2" color="text.secondary">
          Point CCP4i2 at programs installed outside your PATH. Resolution order:
          explicit path → suite directory → extra search paths → PATH.
        </Typography>
        {!editable && (
          <Typography variant="body2" color="warning.main" sx={{ mt: 1 }}>
            Read-only: in a server deployment these are set via environment
            variables.
          </Typography>
        )}
      </Box>

      {/* Explicit per-program fields */}
      <Stack spacing={2}>
        {EXPLICIT_FIELDS.map((f) => (
          <Stack key={f.key} direction="row" spacing={1} alignItems="flex-start">
            <TextField
              label={f.label}
              helperText={f.help}
              value={prefs[f.key] ?? ""}
              disabled={!editable || saving}
              size="small"
              fullWidth
              onChange={(e) =>
                setPrefs((p) => ({ ...p, [f.key]: e.target.value }))
              }
              onBlur={(e) => handleExplicitBlur(f.key, e.target.value.trim())}
            />
            {editable && canBrowse() && (
              <Button
                startIcon={<FolderOpen />}
                onClick={() => handleBrowseExplicit(f.key, f.browse, f.label)}
                disabled={saving}
                sx={{ mt: 0.25, whiteSpace: "nowrap" }}
              >
                Browse
              </Button>
            )}
          </Stack>
        ))}
      </Stack>

      {/* Extra search paths (exePaths / legacy EXEPATHLIST) */}
      <Box>
        <Typography variant="subtitle2" sx={{ mb: 1 }}>
          Extra executable search directories
        </Typography>
        <Stack spacing={1}>
          {exePaths.map((p) => (
            <Stack key={p} direction="row" spacing={1} alignItems="center">
              <Typography variant="body2" sx={{ flexGrow: 1, wordBreak: "break-all" }}>
                {p}
              </Typography>
              {editable && (
                <IconButton size="small" onClick={() => handleRemovePath(p)}>
                  <Delete fontSize="small" />
                </IconButton>
              )}
            </Stack>
          ))}
          {exePaths.length === 0 && (
            <Typography variant="body2" color="text.secondary">
              None.
            </Typography>
          )}
        </Stack>
        {editable && (
          <Stack direction="row" spacing={1} sx={{ mt: 1 }}>
            <TextField
              size="small"
              placeholder="/path/to/bin"
              value={newPath}
              onChange={(e) => setNewPath(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && handleAddPath()}
              fullWidth
            />
            {canBrowse() && (
              <Button
                startIcon={<FolderOpen />}
                onClick={handleBrowseAddPath}
                disabled={saving}
                sx={{ whiteSpace: "nowrap" }}
              >
                Browse
              </Button>
            )}
            <Button
              startIcon={<Add />}
              onClick={handleAddPath}
              disabled={!newPath.trim() || saving}
            >
              Add
            </Button>
          </Stack>
        )}
      </Box>

      {/* Live discovery status */}
      <Box>
        <Stack direction="row" alignItems="center" spacing={1} sx={{ mb: 1 }}>
          <Typography variant="subtitle2">Task programs</Typography>
          {statuses.length > 0 && (
            <Typography variant="caption" color="text.secondary">
              {statuses.filter((s) => s.path != null).length}/{statuses.length} found
            </Typography>
          )}
          <Tooltip title="Re-probe">
            <IconButton size="small" onClick={loadStatuses} disabled={saving}>
              <Refresh fontSize="small" />
            </IconButton>
          </Tooltip>
          {saving && <CircularProgress size={16} />}
        </Stack>
        <Typography variant="caption" color="text.secondary" sx={{ mb: 1, display: "block" }}>
          The external programs each task runs (from the task registry). Hover for
          the resolved path and which tasks use it.
        </Typography>
        <Paper variant="outlined" sx={{ p: 1.5 }}>
          <Box
            sx={{
              display: "grid",
              gridTemplateColumns: "repeat(auto-fill, minmax(200px, 1fr))",
              gap: 1,
            }}
          >
            {[...statuses]
              .sort((a, b) => {
                // Missing first (they're what the user needs to act on), then name.
                const am = a.path == null ? 0 : 1;
                const bm = b.path == null ? 0 : 1;
                return am - bm || a.name.localeCompare(b.name);
              })
              .map((s) => {
                const found = s.path != null;
                const usedBy = s.tasks && s.tasks.length
                  ? `\nused by: ${s.tasks.join(", ")}`
                  : "";
                return (
                  <Tooltip
                    key={s.name}
                    title={
                      (found
                        ? `${s.path} (${SOURCE_LABEL[s.source] ?? s.source})`
                        : "not found on PATH or in your program locations") + usedBy
                    }
                  >
                    <Chip
                      icon={found ? <CheckCircle /> : <ErrorOutline />}
                      color={found ? "success" : "warning"}
                      variant={found ? "filled" : "outlined"}
                      label={s.name}
                      size="small"
                    />
                  </Tooltip>
                );
              })}
          </Box>
        </Paper>
      </Box>
    </Stack>
  );
}
