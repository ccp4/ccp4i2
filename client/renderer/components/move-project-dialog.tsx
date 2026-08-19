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
import { DriveFileMove, FolderOpen } from "@mui/icons-material";
import { useApi } from "../api";
import { Project } from "../types/models";

/** One file the rebase would rewrite. */
interface Rewrite {
  path: string;
  occurrences: number;
}

/** Server summary of a planned or completed move / repair. */
interface MoveSummary {
  source: string;
  destination: string;
  rewrites: Rewrite[];
  total_occurrences: number;
  caches_to_delete: string[];
  rewritten?: number;
  cross_device?: boolean;
  same_filesystem?: boolean;
  stale_roots?: Record<string, number>;
  skipped?: {
    binary: number;
    provenance: number;
    too_large: string[];
  };
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

/** Join a parent directory and a folder name using the parent's separator. */
function joinPath(parent: string, name: string): string {
  const separator = parent.includes("\\") && !parent.includes("/") ? "\\" : "/";
  return `${parent.replace(/[\\/]+$/, "")}${separator}${name}`;
}

/** Last path component, whichever separator the platform uses. */
function baseName(path: string): string {
  const parts = path.replace(/[\\/]+$/, "").split(/[\\/]/);
  return parts[parts.length - 1] ?? path;
}

/** Everything but the last path component. */
function dirName(path: string): string {
  const trimmed = path.replace(/[\\/]+$/, "");
  const index = Math.max(trimmed.lastIndexOf("/"), trimmed.lastIndexOf("\\"));
  return index > 0 ? trimmed.slice(0, index) : "";
}

interface MoveProjectDialogProps {
  project: Project | null;
  open: boolean;
  onClose: () => void;
  /** Called after a successful move or repair so the caller can refresh. */
  onMoved?: () => void;
}

/**
 * Move a project to a new directory, or repair paths left behind by an earlier
 * move.
 *
 * The dialog always runs a dry run first: the user sees exactly which files
 * would be rewritten before anything is touched.
 */
export function MoveProjectDialog({
  project,
  open,
  onClose,
  onMoved,
}: MoveProjectDialogProps) {
  const api = useApi();

  // Kept as two fields rather than one path so there is no ambiguity about
  // whether the project's own folder name is included: the user picks the
  // folder to move *into*, and names the folder that will be created in it.
  const [parentDir, setParentDir] = useState("");
  const [folderName, setFolderName] = useState("");
  const [plan, setPlan] = useState<MoveSummary | null>(null);
  const [result, setResult] = useState<MoveSummary | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [planning, setPlanning] = useState(false);
  const [moving, setMoving] = useState(false);
  const [repairingRoot, setRepairingRoot] = useState<string | null>(null);
  /** Set from the last move/repair result; otherwise the fetched value wins. */
  const [staleRootsOverride, setStaleRootsOverride] = useState<Record<
    string,
    number
  > | null>(null);

  // The projects list uses a lightweight serializer that omits `directory`, so
  // fetch the full record to learn where the project actually lives.
  const { data: detail } = api.get<Project>(
    open && project ? `projects/${project.id}/` : null
  );
  const currentDirectory = detail?.directory ?? project?.directory ?? "";

  const destination =
    parentDir && folderName ? joinPath(parentDir, folderName) : "";
  const unchanged = Boolean(destination) && destination === currentDirectory;

  // Scanning a project for leftover roots costs a full tree walk, so only ask
  // for it while the dialog is actually open.
  const { data: staleRootsResponse, mutate: refreshStaleRoots } = api.get<any>(
    open && project ? `projects/${project.id}/stale_roots/` : null
  );

  const staleRoots: Record<string, number> =
    staleRootsOverride ??
    (staleRootsResponse?.data ?? staleRootsResponse)?.stale_roots ??
    {};

  // Reset whenever the dialog is opened for a (possibly different) project.
  useEffect(() => {
    if (!open) return;
    // Seed from where the project is now: a move usually keeps the folder name
    // and changes the parent, so only one field needs touching.
    setParentDir(dirName(currentDirectory));
    setFolderName(baseName(currentDirectory));
    setPlan(null);
    setResult(null);
    setError(null);
    setStaleRootsOverride(null);
  }, [open, project?.id, currentDirectory]);

  const handleBrowse = useCallback(async () => {
    if (!project) return;
    const parent = await browseDirectory(
      `Choose the folder to move "${project.name}" into. ` +
        "A new folder for the project will be created inside it."
    );
    if (!parent) return;
    setParentDir(parent);
    setPlan(null);
    setResult(null);
    setError(null);
  }, [project]);

  const runDryRun = useCallback(async () => {
    if (!project || !destination) return;
    setPlanning(true);
    setError(null);
    setPlan(null);
    try {
      const response: any = await api.post(`projects/${project.id}/move/`, {
        directory: destination,
        dry_run: true,
      });
      if (response?.success === false) {
        setError(response?.error ?? "Could not plan the move");
        return;
      }
      setPlan(response?.data ?? response);
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setPlanning(false);
    }
  }, [project, destination]);

  // Run the preview automatically once the destination stops changing. Gating
  // the Move button on a manual Preview click just made it look broken, and the
  // dry run is read-only, so there is no reason to make the user ask for it.
  useEffect(() => {
    if (!open || !destination || unchanged || result || moving) return;
    const timer = setTimeout(() => {
      void runDryRun();
    }, 700);
    return () => clearTimeout(timer);
  }, [open, destination, unchanged, result, moving, runDryRun]);

  const runMove = useCallback(async () => {
    if (!project || !destination) return;
    setMoving(true);
    setError(null);
    try {
      const response: any = await api.post(`projects/${project.id}/move/`, {
        directory: destination,
      });
      if (response?.success === false) {
        setError(response?.error ?? "Could not move the project");
        return;
      }
      const summary: MoveSummary = response?.data ?? response;
      setResult(summary);
      setPlan(null);
      setStaleRootsOverride(summary.stale_roots ?? {});
      refreshStaleRoots();
      onMoved?.();
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setMoving(false);
    }
  }, [project, destination, onMoved, refreshStaleRoots]);

  const repairRoot = useCallback(
    async (oldRoot: string) => {
      if (!project) return;
      setRepairingRoot(oldRoot);
      setError(null);
      try {
        const response: any = await api.post(
          `projects/${project.id}/repair_paths/`,
          { old_directory: oldRoot }
        );
        if (response?.success === false) {
          setError(response?.error ?? "Could not repair paths");
          return;
        }
        const summary: MoveSummary = response?.data ?? response;
        setResult(summary);
        setStaleRootsOverride(summary.stale_roots ?? {});
        refreshStaleRoots();
        onMoved?.();
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err));
      } finally {
        setRepairingRoot(null);
      }
    },
    [project, onMoved, refreshStaleRoots]
  );

  const busy = planning || moving || repairingRoot !== null;
  const staleRootEntries = useMemo(
    () => Object.entries(staleRoots),
    [staleRoots]
  );

  if (!project) return null;

  return (
    <Dialog open={open} onClose={busy ? undefined : onClose} maxWidth="md" fullWidth>
      <DialogTitle>
        <Stack direction="row" spacing={1} alignItems="center">
          <DriveFileMove color="primary" />
          <span>Move &quot;{project.name}&quot;</span>
        </Stack>
      </DialogTitle>

      <DialogContent dividers>
        <Typography variant="body2" color="text.secondary" gutterBottom>
          Currently at
        </Typography>
        <Typography
          variant="body2"
          sx={{ fontFamily: "monospace", wordBreak: "break-all", mb: 2 }}
        >
          {currentDirectory}
        </Typography>

        <Stack direction="row" spacing={1} alignItems="flex-start">
          <TextField
            fullWidth
            size="small"
            label="Move into this folder"
            helperText="An existing folder. The project gets its own new folder inside it."
            value={parentDir}
            disabled={busy}
            onChange={(event) => {
              setParentDir(event.target.value);
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

        <TextField
          fullWidth
          size="small"
          label="Project folder name"
          helperText="Created for you - it must not already exist."
          value={folderName}
          disabled={busy}
          onChange={(event) => {
            setFolderName(event.target.value);
            setPlan(null);
            setResult(null);
          }}
          sx={{ mt: 2 }}
        />

        {destination && (
          <Box sx={{ mt: 2 }}>
            <Typography variant="body2" color="text.secondary">
              The project will end up at
            </Typography>
            <Typography
              variant="body2"
              sx={{ fontFamily: "monospace", wordBreak: "break-all" }}
            >
              {destination}
            </Typography>
            {unchanged && (
              <Typography variant="body2" color="warning.main" sx={{ mt: 1 }}>
                That is where the project already is.
              </Typography>
            )}
          </Box>
        )}

        {planning && <LinearProgress sx={{ mt: 2 }} />}

        {error && (
          <Alert severity="error" sx={{ mt: 2 }}>
            {error}
          </Alert>
        )}

        {plan && (
          <Box sx={{ mt: 2 }}>
            <Alert severity="info">
              <AlertTitle>
                {plan.rewrites.length} file
                {plan.rewrites.length === 1 ? "" : "s"} will be updated
              </AlertTitle>
              <Typography variant="body2">
                {plan.total_occurrences} path reference
                {plan.total_occurrences === 1 ? "" : "s"} will be rewritten, and{" "}
                {plan.caches_to_delete.length} cached report
                {plan.caches_to_delete.length === 1 ? "" : "s"} deleted (these
                regenerate on demand).
              </Typography>
              {plan.same_filesystem === false && (
                <Typography variant="body2" sx={{ mt: 1 }}>
                  The destination is on a different filesystem, so the project
                  will be copied and then the original removed. This needs
                  enough free space for a second copy, and will take longer
                  than a move within one disk.
                </Typography>
              )}
              <Typography variant="body2" sx={{ mt: 1 }}>
                Logs and command scripts are left untouched: they record what
                actually ran.
              </Typography>
            </Alert>

            <List dense sx={{ maxHeight: 220, overflow: "auto", mt: 1 }}>
              {plan.rewrites.map((rewrite) => (
                <ListItem key={rewrite.path} disableGutters>
                  <ListItemText
                    primaryTypographyProps={{
                      variant: "body2",
                      sx: { fontFamily: "monospace", wordBreak: "break-all" },
                    }}
                    primary={rewrite.path}
                    secondary={`${rewrite.occurrences} reference${
                      rewrite.occurrences === 1 ? "" : "s"
                    }`}
                  />
                </ListItem>
              ))}
            </List>
          </Box>
        )}

        {result && (
          <Alert severity="success" sx={{ mt: 2 }}>
            <AlertTitle>Done</AlertTitle>
            <Typography variant="body2">
              {result.rewritten ?? 0} file
              {(result.rewritten ?? 0) === 1 ? "" : "s"} updated. The project is
              now at{" "}
              <Box component="span" sx={{ fontFamily: "monospace" }}>
                {result.destination}
              </Box>
              .
            </Typography>
          </Alert>
        )}

        {staleRootEntries.length > 0 && (
          <>
            <Divider sx={{ my: 2 }} />
            <Alert severity="warning">
              <AlertTitle>This project also refers to other locations</AlertTitle>
              <Typography variant="body2" gutterBottom>
                Paths from an earlier move, an import from another machine, or a
                drive that has since been renamed. Rewrite them to point at the
                project&apos;s current directory:
              </Typography>
              <Stack spacing={1} sx={{ mt: 1 }}>
                {staleRootEntries.map(([root, count]) => (
                  <Stack
                    key={root}
                    direction="row"
                    spacing={1}
                    alignItems="center"
                    flexWrap="wrap"
                  >
                    <Box
                      component="span"
                      sx={{
                        fontFamily: "monospace",
                        fontSize: "0.8rem",
                        wordBreak: "break-all",
                        flexGrow: 1,
                      }}
                    >
                      {root}
                    </Box>
                    <Chip size="small" label={`${count} refs`} />
                    <Button
                      size="small"
                      variant="outlined"
                      disabled={busy}
                      onClick={() => repairRoot(root)}
                      startIcon={
                        repairingRoot === root ? (
                          <CircularProgress size={14} />
                        ) : undefined
                      }
                    >
                      Repair
                    </Button>
                  </Stack>
                ))}
              </Stack>
            </Alert>
          </>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose} disabled={busy}>
          {result ? "Close" : "Cancel"}
        </Button>
        {/* The preview runs itself; this is for retrying after a bad path. */}
        <Button
          onClick={runDryRun}
          disabled={busy || !destination || unchanged || Boolean(result)}
        >
          Re-check
        </Button>
        <Button
          variant="contained"
          onClick={runMove}
          disabled={busy || !destination || unchanged || Boolean(result)}
          startIcon={moving ? <CircularProgress size={16} /> : <DriveFileMove />}
        >
          Move project
        </Button>
      </DialogActions>
    </Dialog>
  );
}

export default MoveProjectDialog;
