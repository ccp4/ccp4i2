"use client";
import React, { useCallback, useState } from "react";
import {
  Alert,
  Box,
  Button,
  Checkbox,
  Chip,
  LinearProgress,
  List,
  ListItem,
  ListItemText,
  Paper,
  Stack,
  Typography,
} from "@mui/material";
import { DriveFolderUpload, FolderOpen } from "@mui/icons-material";
import { useApi } from "../api";
import { apiGet } from "../api-fetch";

/** One project directory found by the scan, plus what a dry run made of it. */
interface Candidate {
  directory: string;
  project_name: string | null;
  project_uuid: string | null;
  recorded_directory: string | null;
  jobs: number;
  files: number;
  relocated?: boolean;
  /** Filled in by the dry run: why this one cannot be imported. */
  skipped_reason?: string | null;
}

function canBrowse(): boolean {
  return typeof window !== "undefined" && Boolean(window.electronAPI?.invoke);
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
        buttonLabel: "Use this folder",
      })) ?? null
    );
  } catch {
    return null;
  }
}

interface ImportProjectDirectoryProps {
  /** Called after at least one project has been imported. */
  onImported?: () => void;
}

/**
 * Import a project that is already unpacked on disk, in place.
 *
 * The other half of the import page takes a zip and copies its contents into
 * the project store. This takes a folder and adopts it where it lies: nothing
 * is copied and nothing is renumbered, so the job directories keep the names
 * the paths inside their own `params.xml` files already point at.
 *
 * Desktop only, and the constraint is real rather than cosmetic — the folder
 * has to be on the same machine as the server for adopting it to mean
 * anything. The server refuses these routes off the desktop as well.
 */
export const ImportProjectDirectory: React.FC<ImportProjectDirectoryProps> = ({
  onImported,
}) => {
  const api = useApi();

  const [path, setPath] = useState<string | null>(null);
  const [candidates, setCandidates] = useState<Candidate[] | null>(null);
  const [chosen, setChosen] = useState<Set<string>>(new Set());
  const [error, setError] = useState<string | null>(null);
  const [looking, setLooking] = useState(false);
  const [importing, setImporting] = useState(false);
  const [done, setDone] = useState<string | null>(null);

  /**
   * Ask what is in a folder, then dry-run each candidate.
   *
   * The scan only reads each `DATABASE.db.xml`; it is the dry run that knows
   * whether the name or the UUID is already taken. Doing both here means the
   * clash is on screen before the user commits to anything, rather than
   * arriving as a skip in the results.
   */
  const look = useCallback(
    async (directory: string) => {
      setLooking(true);
      setError(null);
      setCandidates(null);
      setDone(null);
      try {
        const payload: any = await apiGet(
          `projects/restorable/?scan=${encodeURIComponent(directory)}`
        );
        if (payload?.success === false) {
          setError(payload?.error ?? "Could not look in that folder");
          return;
        }
        const found: Candidate[] = (payload?.data ?? payload)?.candidates ?? [];
        const usable = found.filter((one) => one.project_uuid !== null);

        const checked = await Promise.all(
          usable.map(async (one) => {
            try {
              const response: any = await api.post("projects/restore/", {
                source: "directory",
                path: one.directory,
                dry_run: true,
              });
              const result = response?.data ?? response;
              const report = (result?.skipped ?? []).find(
                (entry: Candidate) => entry.directory === one.directory
              );
              return { ...one, skipped_reason: report?.skipped_reason ?? null };
            } catch {
              // A dry run that will not run is not a reason to hide the
              // project: let the user try, and report the real failure then.
              return { ...one, skipped_reason: null };
            }
          })
        );

        setCandidates(checked);
        setChosen(
          new Set(
            checked.filter((one) => !one.skipped_reason).map((one) => one.directory)
          )
        );
        if (found.length === 0) {
          setError(
            "No CCP4i2 project in that folder. Pick either a project folder " +
              "itself or the folder your projects sit in."
          );
        }
      } catch (err) {
        setError(err instanceof Error ? err.message : String(err));
      } finally {
        setLooking(false);
      }
    },
    [api]
  );

  const handleBrowse = useCallback(async () => {
    const picked = await browseDirectory("Choose a CCP4i2 project folder");
    if (!picked) return;
    setPath(picked);
    void look(picked);
  }, [look]);

  const handleImport = useCallback(async () => {
    if (!candidates) return;
    setImporting(true);
    setError(null);
    const imported: string[] = [];
    const failed: string[] = [];
    try {
      for (const one of candidates) {
        if (!chosen.has(one.directory)) continue;
        const response: any = await api.post("projects/restore/", {
          source: "directory",
          path: one.directory,
        });
        const result = response?.data ?? response;
        if (response?.success === false) {
          failed.push(`${one.project_name}: ${response?.error}`);
        } else if ((result?.restored ?? []).length > 0) {
          imported.push(one.project_name ?? one.directory);
        } else {
          const report = (result?.skipped ?? [])[0];
          failed.push(
            `${one.project_name}: ${report?.skipped_reason ?? "not imported"}`
          );
        }
      }
      if (failed.length > 0) setError(failed.join("; "));
      if (imported.length > 0) {
        setDone(
          `Imported ${imported.join(", ")}. ${
            imported.length === 1 ? "It stays" : "They stay"
          } where ${imported.length === 1 ? "it is" : "they are"} on disk.`
        );
        setCandidates(null);
        onImported?.();
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setImporting(false);
    }
  }, [api, candidates, chosen, onImported]);

  if (!canBrowse()) return null;

  const busy = looking || importing;
  const importable = (candidates ?? []).filter((one) => !one.skipped_reason);

  return (
    <Paper variant="outlined" sx={{ padding: 2 }}>
      <Stack spacing={2}>
        <Stack direction="row" spacing={1} alignItems="center">
          <DriveFolderUpload color="primary" />
          <Typography variant="h6">…or a project folder</Typography>
        </Stack>
        <Typography variant="body2" color="text.secondary">
          If the project is already unpacked on this machine, import the folder
          directly — there is no need to zip it up first. It is adopted where it
          lies: nothing is copied and nothing is renumbered.
        </Typography>

        <Stack direction="row" spacing={2} alignItems="center">
          <Button
            variant="outlined"
            startIcon={<FolderOpen />}
            onClick={handleBrowse}
            disabled={busy}
          >
            Choose folder
          </Button>
          {path && (
            <Typography
              variant="body2"
              sx={{ fontFamily: "monospace", wordBreak: "break-all" }}
            >
              {path}
            </Typography>
          )}
        </Stack>

        {busy && <LinearProgress />}
        {error && <Alert severity="error">{error}</Alert>}
        {done && <Alert severity="success">{done}</Alert>}

        {candidates && candidates.length > 0 && (
          <>
            <List dense>
              {candidates.map((one) => (
                <ListItem
                  key={one.directory}
                  secondaryAction={
                    one.skipped_reason ? (
                      <Chip size="small" color="warning" label="Cannot import" />
                    ) : (
                      <Checkbox
                        edge="end"
                        checked={chosen.has(one.directory)}
                        disabled={busy}
                        onChange={(event) => {
                          const next = new Set(chosen);
                          if (event.target.checked) next.add(one.directory);
                          else next.delete(one.directory);
                          setChosen(next);
                        }}
                      />
                    )
                  }
                >
                  <ListItemText
                    primary={one.project_name ?? one.directory}
                    secondary={
                      <>
                        <Box component="span">
                          {one.jobs} job{one.jobs === 1 ? "" : "s"}, {one.files}{" "}
                          file{one.files === 1 ? "" : "s"}
                        </Box>
                        {one.relocated && one.recorded_directory && (
                          <Box component="span" sx={{ display: "block" }}>
                            Last recorded at {one.recorded_directory}; paths will
                            be repointed here.
                          </Box>
                        )}
                        {one.skipped_reason && (
                          <Box
                            component="span"
                            sx={{ display: "block", color: "warning.main" }}
                          >
                            {one.skipped_reason}
                          </Box>
                        )}
                      </>
                    }
                  />
                </ListItem>
              ))}
            </List>
            <Button
              variant="contained"
              disabled={busy || chosen.size === 0}
              onClick={handleImport}
            >
              Import {chosen.size > 1 ? `${chosen.size} projects` : "project"} in
              place
            </Button>
            {importable.length === 0 && (
              <Typography variant="body2" color="text.secondary">
                Nothing here can be imported: a project is identified by its name
                and its UUID, and both must be free.
              </Typography>
            )}
          </>
        )}
      </Stack>
    </Paper>
  );
};

export default ImportProjectDirectory;
