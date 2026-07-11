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
} from "@mui/icons-material";
import React, { useCallback, useEffect, useState } from "react";
import { apiGet, apiPatch } from "../api-fetch";

interface ProgramStatus {
  name: string;
  path: string | null;
  source: string;
}

// Explicit path fields shown as their own inputs (the rest go via exePaths).
const EXPLICIT_FIELDS: { key: string; label: string; help: string }[] = [
  { key: "COOT_EXECUTABLE", label: "Coot executable", help: "Full path to the coot binary" },
  { key: "CCP4MG_EXECUTABLE", label: "CCP4mg executable", help: "Full path to the ccp4mg binary" },
  { key: "SHELXDIR", label: "SHELX directory", help: "Directory containing shelxc/d/e/l" },
  { key: "DIALSDIR", label: "DIALS directory", help: "Directory containing dials binaries" },
  { key: "BUSTERDIR", label: "BUSTER directory", help: "Directory containing the BUSTER refine binary" },
];

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
  const handleAddPath = () => {
    const p = newPath.trim();
    if (!p || exePaths.includes(p)) return;
    const next = [...exePaths, p];
    setExePaths(next);
    setNewPath("");
    persist({ exePaths: next });
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
          <TextField
            key={f.key}
            label={f.label}
            helperText={f.help}
            defaultValue={prefs[f.key] ?? ""}
            disabled={!editable || saving}
            size="small"
            fullWidth
            onBlur={(e) => handleExplicitBlur(f.key, e.target.value.trim())}
          />
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
          <Typography variant="subtitle2">Discovered programs</Typography>
          <Tooltip title="Re-probe">
            <IconButton size="small" onClick={loadStatuses} disabled={saving}>
              <Refresh fontSize="small" />
            </IconButton>
          </Tooltip>
          {saving && <CircularProgress size={16} />}
        </Stack>
        <Paper variant="outlined" sx={{ p: 1.5 }}>
          <Box
            sx={{
              display: "grid",
              gridTemplateColumns: "repeat(auto-fill, minmax(220px, 1fr))",
              gap: 1,
            }}
          >
            {statuses.map((s) => {
              const found = s.path != null;
              return (
                <Tooltip
                  key={s.name}
                  title={found ? `${s.path} (${SOURCE_LABEL[s.source] ?? s.source})` : "not found"}
                >
                  <Chip
                    icon={found ? <CheckCircle /> : <ErrorOutline />}
                    color={found ? "success" : "default"}
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
