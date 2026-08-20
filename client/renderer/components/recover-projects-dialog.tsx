"use client";
import React, { useCallback, useEffect, useState } from "react";
import {
  Alert,
  AlertTitle,
  Box,
  Button,
  Checkbox,
  Chip,
  CircularProgress,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Divider,
  FormControlLabel,
  LinearProgress,
  List,
  ListItem,
  ListItemText,
  Stack,
  Typography,
} from "@mui/material";
import {
  CheckCircle,
  FolderOpen,
  Restore,
  WarningAmber,
} from "@mui/icons-material";
import { useApi } from "../api";
import { apiGet } from "../api-fetch";

/** One project directory that could be restored. */
interface Candidate {
  directory: string;
  project_name: string | null;
  project_uuid: string | null;
  recorded_directory: string | null;
  jobs: number;
  files: number;
  relocated: boolean;
  restored: boolean;
  skipped_reason: string | null;
  paths_rewritten: number;
  missing_job_directories: string[];
  warnings: string[];
}

interface RestorableResponse {
  source: { kind: string; path: string; entries?: number };
  candidates: Candidate[];
  restorable: number;
}

interface RestoreResult {
  restored: Candidate[];
  skipped: Candidate[];
  total: number;
  dry_run: boolean;
}

async function browseDirectory(message: string): Promise<string | null> {
  const api = typeof window !== "undefined" ? window.electronAPI : undefined;
  if (!api?.invoke) return null;
  try {
    return (
      (await api.invoke("browse-path", {
        mode: "directory",
        title: message,
        message,
        buttonLabel: "Look here",
      })) ?? null
    );
  } catch {
    return null;
  }
}

function canBrowse(): boolean {
  return typeof window !== "undefined" && Boolean(window.electronAPI?.invoke);
}

interface RecoverProjectsDialogProps {
  open: boolean;
  onClose: () => void;
  onRecovered?: () => void;
}

/**
 * Rebuild the database from the snapshots CCP4i2 leaves in project directories.
 *
 * For a database that has been lost or corrupted. Each project directory holds
 * a `DATABASE.db.xml` describing everything the database knew about it, and a
 * registry outside the database records where the projects are — so recovery
 * normally needs nothing but a confirmation. If the registry has gone too,
 * pointing this at the folder the projects live in finds them directly.
 *
 * Nothing is moved and nothing already in the database is touched unless the
 * user explicitly asks: restoring over a live project would lose work rather
 * than recover it.
 */
export function RecoverProjectsDialog({
  open,
  onClose,
  onRecovered,
}: RecoverProjectsDialogProps) {
  const api = useApi();

  const [scanPath, setScanPath] = useState<string | null>(null);
  const [found, setFound] = useState<RestorableResponse | null>(null);
  const [result, setResult] = useState<RestoreResult | null>(null);
  const [replace, setReplace] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [looking, setLooking] = useState(false);
  const [restoring, setRestoring] = useState(false);

  const look = useCallback(
    async (path: string | null) => {
      setLooking(true);
      setError(null);
      setFound(null);
      try {
        const query = path ? `?scan=${encodeURIComponent(path)}` : "";
        // A one-off read, not a subscription: SWR would cache a view of the
        // filesystem that goes stale the moment the user looks elsewhere.
        const payload: any = await apiGet(`projects/restorable/${query}`);
        if (payload?.success === false) {
          setError(payload?.error ?? "Could not look for projects");
          return;
        }
        setFound(payload?.data ?? payload);
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err));
      } finally {
        setLooking(false);
      }
    },
    []
  );

  useEffect(() => {
    if (!open) return;
    setScanPath(null);
    setFound(null);
    setResult(null);
    setError(null);
    setReplace(false);
    void look(null);
  }, [open, look]);

  const handleBrowse = useCallback(async () => {
    const chosen = await browseDirectory(
      "Choose the folder your CCP4i2 projects are in"
    );
    if (!chosen) return;
    setScanPath(chosen);
    setResult(null);
    void look(chosen);
  }, [look]);

  const restore = useCallback(async () => {
    setRestoring(true);
    setError(null);
    try {
      const response: any = await api.post("projects/restore/", {
        source: scanPath ? "scan" : "registry",
        ...(scanPath ? { path: scanPath } : {}),
        replace,
      });
      if (response?.success === false) {
        setError(response?.error ?? "Could not restore the projects");
        return;
      }
      setResult(response?.data ?? response);
      onRecovered?.();
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setRestoring(false);
    }
  }, [api, scanPath, replace, onRecovered]);

  const busy = looking || restoring;
  const candidates = found?.candidates ?? [];
  const usable = candidates.filter((c) => c.project_uuid !== null);
  const unusable = candidates.filter((c) => c.project_uuid === null);

  return (
    <Dialog
      open={open}
      onClose={busy ? undefined : onClose}
      maxWidth="md"
      fullWidth
    >
      <DialogTitle>
        <Stack direction="row" spacing={1} alignItems="center">
          <Restore color="primary" />
          <span>Recover projects from disk</span>
        </Stack>
      </DialogTitle>

      <DialogContent dividers>
        <Typography variant="body2" color="text.secondary">
          CCP4i2 keeps a record of each project inside its own folder, so a lost
          or damaged database can be rebuilt from what is still on disk. Nothing
          is moved, and projects already in the database are left alone.
        </Typography>

        <Stack
          direction="row"
          spacing={1}
          alignItems="center"
          sx={{ mt: 2, flexWrap: "wrap" }}
        >
          <Typography variant="body2" sx={{ flexGrow: 1 }}>
            {scanPath ? (
              <>
                Looking in{" "}
                <Box component="span" sx={{ fontFamily: "monospace" }}>
                  {scanPath}
                </Box>
              </>
            ) : (
              "Using CCP4i2's own record of where your projects are"
            )}
          </Typography>
          {canBrowse() && (
            <Button
              size="small"
              variant="outlined"
              startIcon={<FolderOpen />}
              onClick={handleBrowse}
              disabled={busy}
            >
              Look elsewhere
            </Button>
          )}
        </Stack>

        {looking && <LinearProgress sx={{ mt: 2 }} />}

        {error && (
          <Alert severity="error" sx={{ mt: 2 }}>
            {error}
          </Alert>
        )}

        {found && !result && (
          <Box sx={{ mt: 2 }}>
            {usable.length === 0 ? (
              <Alert severity="warning">
                <AlertTitle>Nothing found to recover</AlertTitle>
                {scanPath
                  ? "No project folders with a recovery record were found there."
                  : "CCP4i2 has no record of where your projects are. If you " +
                    "know the folder they live in, use Look elsewhere."}
              </Alert>
            ) : (
              <Alert severity="info">
                <AlertTitle>
                  {usable.length} project{usable.length === 1 ? "" : "s"} can be
                  recovered
                </AlertTitle>
                Each will be rebuilt exactly as it was, under the job numbers
                its folders already use.
              </Alert>
            )}

            <List dense sx={{ maxHeight: 260, overflow: "auto", mt: 1 }}>
              {usable.map((candidate) => (
                <ListItem key={candidate.directory} disableGutters>
                  <CheckCircle
                    fontSize="small"
                    color="success"
                    sx={{ mr: 1, flexShrink: 0 }}
                  />
                  <ListItemText
                    primary={
                      <Stack direction="row" spacing={1} alignItems="center">
                        <span>{candidate.project_name}</span>
                        <Chip
                          size="small"
                          label={`${candidate.jobs} job${
                            candidate.jobs === 1 ? "" : "s"
                          }`}
                        />
                        {candidate.relocated && (
                          <Chip size="small" color="info" label="moved" />
                        )}
                      </Stack>
                    }
                    secondary={candidate.directory}
                    secondaryTypographyProps={{
                      sx: { fontFamily: "monospace", wordBreak: "break-all" },
                    }}
                  />
                </ListItem>
              ))}
              {unusable.map((candidate) => (
                <ListItem key={candidate.directory} disableGutters>
                  <WarningAmber
                    fontSize="small"
                    color="warning"
                    sx={{ mr: 1, flexShrink: 0 }}
                  />
                  <ListItemText
                    primary={candidate.directory}
                    secondary={candidate.skipped_reason}
                    primaryTypographyProps={{
                      sx: { fontFamily: "monospace", wordBreak: "break-all" },
                    }}
                  />
                </ListItem>
              ))}
            </List>

            {usable.length > 0 && (
              <FormControlLabel
                sx={{ mt: 1 }}
                control={
                  <Checkbox
                    checked={replace}
                    disabled={busy}
                    onChange={(event) => setReplace(event.target.checked)}
                  />
                }
                label={
                  <Typography variant="body2">
                    Also rebuild projects already in the database, discarding
                    what it currently holds for them
                  </Typography>
                }
              />
            )}
          </Box>
        )}

        {result && (
          <Box sx={{ mt: 2 }}>
            <Alert severity={result.restored.length > 0 ? "success" : "warning"}>
              <AlertTitle>
                Recovered {result.restored.length} project
                {result.restored.length === 1 ? "" : "s"}
              </AlertTitle>
              {result.restored.some((c) => c.paths_rewritten > 0) && (
                <Typography variant="body2">
                  Some had moved since their record was written; the paths
                  inside their files were updated to match.
                </Typography>
              )}
            </Alert>

            {result.restored.some(
              (c) => c.missing_job_directories.length > 0
            ) && (
              <>
                <Divider sx={{ my: 2 }} />
                <Alert severity="warning">
                  <AlertTitle>Some job folders are not on disk</AlertTitle>
                  <Typography variant="body2" gutterBottom>
                    These jobs are recorded but their folders are missing, so
                    their results cannot be opened.
                  </Typography>
                  {result.restored
                    .filter((c) => c.missing_job_directories.length > 0)
                    .map((c) => (
                      <Typography variant="body2" key={c.directory}>
                        {c.project_name}: job{" "}
                        {c.missing_job_directories.slice(0, 8).join(", ")}
                      </Typography>
                    ))}
                </Alert>
              </>
            )}

            {result.skipped.filter((c) => c.skipped_reason).length > 0 && (
              <List dense sx={{ maxHeight: 200, overflow: "auto", mt: 1 }}>
                {result.skipped
                  .filter((c) => c.skipped_reason)
                  .map((c) => (
                    <ListItem key={c.directory} disableGutters>
                      <ListItemText
                        primary={c.project_name ?? c.directory}
                        secondary={c.skipped_reason}
                      />
                    </ListItem>
                  ))}
              </List>
            )}
          </Box>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose} disabled={busy}>
          {result ? "Close" : "Cancel"}
        </Button>
        <Button
          variant="contained"
          onClick={restore}
          disabled={busy || usable.length === 0 || Boolean(result)}
          startIcon={restoring ? <CircularProgress size={16} /> : <Restore />}
        >
          Recover
        </Button>
      </DialogActions>
    </Dialog>
  );
}

export default RecoverProjectsDialog;
