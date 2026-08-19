"use client";
import React, { useCallback, useEffect, useMemo, useState } from "react";
import {
  Alert,
  AlertTitle,
  Box,
  Button,
  Chip,
  CircularProgress,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Divider,
  LinearProgress,
  List,
  ListItem,
  ListItemText,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import { CheckCircle, FolderOpen, LinkOff } from "@mui/icons-material";
import { useApi } from "../api";

interface BrokenProject {
  id: number;
  uuid: string;
  name: string;
  directory: string;
}

interface MappedProject {
  id: number;
  name: string;
  old_directory: string;
  new_directory: string;
  reason?: string;
  rewritten?: number;
  error?: string;
}

interface RootRebaseSummary {
  old_root: string;
  new_root: string;
  matched: MappedProject[];
  missing: MappedProject[];
  unaffected: { id: number; name: string }[];
  relocated?: MappedProject[];
  failed?: MappedProject[];
  preference_updated?: boolean;
}

/** Open the native directory picker; returns null in the web build. */
async function browseDirectory(message: string): Promise<string | null> {
  const api = typeof window !== "undefined" ? window.electronAPI : undefined;
  if (!api?.invoke) return null;
  try {
    return (
      (await api.invoke("browse-path", {
        mode: "directory",
        title: message,
        message,
        buttonLabel: "Use this folder",
      })) ?? null
    );
  } catch {
    return null;
  }
}

function canBrowse(): boolean {
  return typeof window !== "undefined" && Boolean(window.electronAPI?.invoke);
}

/** Everything but the last path component. */
function dirName(path: string): string {
  const trimmed = path.replace(/[\\/]+$/, "");
  const index = Math.max(trimmed.lastIndexOf("/"), trimmed.lastIndexOf("\\"));
  return index > 0 ? trimmed.slice(0, index) : "";
}

/**
 * The deepest folder all these paths share.
 *
 * When a drive is renamed every project moves together, so their common parent
 * is almost always the root that changed — which saves the user working out
 * what to type.
 */
function commonParent(paths: string[]): string {
  if (paths.length === 0) return "";
  const split = (p: string) => dirName(p).split(/[\\/]/);
  let prefix = split(paths[0]);
  for (const path of paths.slice(1)) {
    const parts = split(path);
    let i = 0;
    while (i < prefix.length && i < parts.length && prefix[i] === parts[i]) i++;
    prefix = prefix.slice(0, i);
  }
  const joined = prefix.join("/");
  return joined === "" && paths[0].startsWith("/") ? "/" : joined;
}

interface ReconnectProjectsDialogProps {
  open: boolean;
  onClose: () => void;
  brokenProjects: BrokenProject[];
  /** Called after anything changes so the caller can refresh. */
  onReconnected?: () => void;
}

/**
 * Reconnect projects whose recorded directory is no longer on disk.
 *
 * A renamed drive or a projects folder moved outside CCP4i2 invalidates every
 * project at once, so this works on the storage root rather than project by
 * project: say where the projects used to be and where they are now, and each
 * one is re-pointed and its internal paths rebased. Nothing is moved.
 */
export function ReconnectProjectsDialog({
  open,
  onClose,
  brokenProjects,
  onReconnected,
}: ReconnectProjectsDialogProps) {
  const api = useApi();

  const [oldRoot, setOldRoot] = useState("");
  const [newRoot, setNewRoot] = useState("");
  const [plan, setPlan] = useState<RootRebaseSummary | null>(null);
  const [result, setResult] = useState<RootRebaseSummary | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [planning, setPlanning] = useState(false);
  const [applying, setApplying] = useState(false);
  const [locating, setLocating] = useState<number | null>(null);

  const suggestedOldRoot = useMemo(
    () => commonParent(brokenProjects.map((p) => p.directory)),
    [brokenProjects]
  );

  useEffect(() => {
    if (!open) return;
    setOldRoot(suggestedOldRoot);
    setNewRoot("");
    setPlan(null);
    setResult(null);
    setError(null);
  }, [open, suggestedOldRoot]);

  const runDryRun = useCallback(async () => {
    if (!oldRoot || !newRoot) return;
    setPlanning(true);
    setError(null);
    setPlan(null);
    try {
      const response: any = await api.post("projects/rebase_root/", {
        old_root: oldRoot,
        new_root: newRoot,
        dry_run: true,
      });
      if (response?.success === false) {
        setError(response?.error ?? "Could not work out the mapping");
        return;
      }
      setPlan(response?.data ?? response);
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setPlanning(false);
    }
  }, [oldRoot, newRoot]);

  // Same reasoning as the move dialog: the check is read-only, so there is no
  // point making the user ask for it.
  useEffect(() => {
    if (!open || !oldRoot || !newRoot || oldRoot === newRoot || result) return;
    const timer = setTimeout(() => {
      void runDryRun();
    }, 700);
    return () => clearTimeout(timer);
  }, [open, oldRoot, newRoot, result, runDryRun]);

  const handleBrowse = useCallback(async () => {
    const chosen = await browseDirectory(
      "Choose the folder the projects are in now"
    );
    if (!chosen) return;
    setNewRoot(chosen);
    setPlan(null);
    setResult(null);
    setError(null);
  }, []);

  const apply = useCallback(async () => {
    if (!oldRoot || !newRoot) return;
    setApplying(true);
    setError(null);
    try {
      const response: any = await api.post("projects/rebase_root/", {
        old_root: oldRoot,
        new_root: newRoot,
      });
      if (response?.success === false) {
        setError(response?.error ?? "Could not reconnect the projects");
        return;
      }
      setResult(response?.data ?? response);
      setPlan(null);
      onReconnected?.();
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setApplying(false);
    }
  }, [oldRoot, newRoot, onReconnected]);

  /** Point one stubborn project at a folder the user picks by hand. */
  const locate = useCallback(
    async (project: BrokenProject) => {
      const chosen = await browseDirectory(
        `Find the folder for "${project.name}"`
      );
      if (!chosen) return;
      setLocating(project.id);
      setError(null);
      try {
        const response: any = await api.post(
          `projects/${project.id}/relocate/`,
          { directory: chosen }
        );
        if (response?.success === false) {
          setError(response?.error ?? `Could not reconnect ${project.name}`);
          return;
        }
        onReconnected?.();
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err));
      } finally {
        setLocating(null);
      }
    },
    [onReconnected]
  );

  const busy = planning || applying || locating !== null;
  const stillBroken = result
    ? [...(result.missing ?? []), ...(result.failed ?? [])]
    : [];

  return (
    <Dialog
      open={open}
      onClose={busy ? undefined : onClose}
      maxWidth="md"
      fullWidth
    >
      <DialogTitle>
        <Stack direction="row" spacing={1} alignItems="center">
          <LinkOff color="warning" />
          <span>Reconnect projects</span>
        </Stack>
      </DialogTitle>

      <DialogContent dividers>
        <Typography variant="body2" color="text.secondary" gutterBottom>
          {brokenProjects.length} project
          {brokenProjects.length === 1 ? "" : "s"} cannot be found where CCP4i2
          recorded {brokenProjects.length === 1 ? "it" : "them"}. If a drive was
          renamed or the projects folder was moved, say where they are now and
          they will be reconnected. Nothing is moved.
        </Typography>

        <TextField
          fullWidth
          size="small"
          label="They used to be in"
          helperText="Taken from the projects that are missing - edit if it is not right."
          value={oldRoot}
          disabled={busy}
          onChange={(event) => {
            setOldRoot(event.target.value);
            setPlan(null);
            setResult(null);
          }}
          sx={{ mt: 2 }}
        />

        <Stack direction="row" spacing={1} alignItems="flex-start" sx={{ mt: 2 }}>
          <TextField
            fullWidth
            size="small"
            label="They are now in"
            value={newRoot}
            disabled={busy}
            onChange={(event) => {
              setNewRoot(event.target.value);
              setPlan(null);
              setResult(null);
            }}
          />
          {canBrowse() && (
            <Button
              variant="outlined"
              startIcon={<FolderOpen />}
              onClick={handleBrowse}
              disabled={busy}
              sx={{ whiteSpace: "nowrap", mt: 0.25 }}
            >
              Browse
            </Button>
          )}
        </Stack>

        {planning && <LinearProgress sx={{ mt: 2 }} />}

        {error && (
          <Alert severity="error" sx={{ mt: 2 }}>
            {error}
          </Alert>
        )}

        {plan && (
          <Box sx={{ mt: 2 }}>
            <Alert severity={plan.matched.length > 0 ? "info" : "warning"}>
              <AlertTitle>
                {plan.matched.length} of {plan.matched.length + plan.missing.length}{" "}
                found in the new location
              </AlertTitle>
              {plan.matched.length > 0 && (
                <Typography variant="body2">
                  Each will be re-pointed and the paths inside its files
                  rebased. Projects elsewhere are left alone.
                </Typography>
              )}
              {plan.matched.length === 0 && (
                <Typography variant="body2">
                  Nothing was found under that folder. Check it is the folder
                  that directly contains the project folders.
                </Typography>
              )}
            </Alert>

            <List dense sx={{ maxHeight: 240, overflow: "auto", mt: 1 }}>
              {plan.matched.map((entry) => (
                <ListItem key={entry.id} disableGutters>
                  <CheckCircle
                    fontSize="small"
                    color="success"
                    sx={{ mr: 1, flexShrink: 0 }}
                  />
                  <ListItemText
                    primary={entry.name}
                    secondary={entry.new_directory}
                    secondaryTypographyProps={{
                      sx: { fontFamily: "monospace", wordBreak: "break-all" },
                    }}
                  />
                </ListItem>
              ))}
              {plan.missing.map((entry) => (
                <ListItem key={entry.id} disableGutters>
                  <LinkOff
                    fontSize="small"
                    color="warning"
                    sx={{ mr: 1, flexShrink: 0 }}
                  />
                  <ListItemText
                    primary={entry.name}
                    secondary={`${entry.reason} at ${entry.new_directory}`}
                    secondaryTypographyProps={{
                      sx: { fontFamily: "monospace", wordBreak: "break-all" },
                    }}
                  />
                </ListItem>
              ))}
            </List>
          </Box>
        )}

        {result && (
          <Alert severity="success" sx={{ mt: 2 }}>
            <AlertTitle>
              Reconnected {result.relocated?.length ?? 0} project
              {(result.relocated?.length ?? 0) === 1 ? "" : "s"}
            </AlertTitle>
            {result.preference_updated && (
              <Typography variant="body2">
                Your default projects folder was updated to match, so new
                projects will be created in the right place.
              </Typography>
            )}
          </Alert>
        )}

        {stillBroken.length > 0 && (
          <>
            <Divider sx={{ my: 2 }} />
            <Alert severity="warning">
              <AlertTitle>Still not found</AlertTitle>
              <Typography variant="body2" gutterBottom>
                These did not move with the rest. Point each one at its folder:
              </Typography>
              <Stack spacing={1} sx={{ mt: 1 }}>
                {stillBroken.map((entry) => {
                  const project = brokenProjects.find((p) => p.id === entry.id);
                  return (
                    <Stack
                      key={entry.id}
                      direction="row"
                      spacing={1}
                      alignItems="center"
                      flexWrap="wrap"
                    >
                      <Box component="span" sx={{ flexGrow: 1 }}>
                        {entry.name}
                      </Box>
                      <Chip
                        size="small"
                        label={entry.error ?? entry.reason ?? "not found"}
                      />
                      {canBrowse() && project && (
                        <Button
                          size="small"
                          variant="outlined"
                          disabled={busy}
                          onClick={() => locate(project)}
                          startIcon={
                            locating === entry.id ? (
                              <CircularProgress size={14} />
                            ) : (
                              <FolderOpen />
                            )
                          }
                        >
                          Locate
                        </Button>
                      )}
                    </Stack>
                  );
                })}
              </Stack>
            </Alert>
          </>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose} disabled={busy}>
          {result ? "Close" : "Cancel"}
        </Button>
        <Button
          onClick={runDryRun}
          disabled={busy || !oldRoot || !newRoot || Boolean(result)}
        >
          Re-check
        </Button>
        <Button
          variant="contained"
          onClick={apply}
          disabled={
            busy || !plan || plan.matched.length === 0 || Boolean(result)
          }
          startIcon={applying ? <CircularProgress size={16} /> : undefined}
        >
          Reconnect
        </Button>
      </DialogActions>
    </Dialog>
  );
}

export default ReconnectProjectsDialog;
