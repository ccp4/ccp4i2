"use client";
import React, { useCallback } from "react";
import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Alert,
  Button,
  CircularProgress,
  Dialog,
  DialogContent,
  DialogTitle,
  Paper,
  Stack,
  Switch,
  Typography,
} from "@mui/material";
import {
  Cancel,
  Check,
  RadioButtonUnchecked,
  ExpandMore,
  RocketLaunch,
} from "@mui/icons-material";
import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import { useCCP4i2Window } from "../app-context";
import { usePopcorn } from "../providers/popcorn-provider";
import { ThemeToggle } from "./theme-toggle";
import { useTheme } from "../theme/theme-provider";

/**
 * The "get ready" body of the launch screen. Owns the Electron setup/launch
 * lifecycle (IPC to the main process) but presents it as one job: get the user
 * into their work.
 *
 * - When everything is ready, the setup detail collapses and a single "Enter
 *   CCP4i2" action launches the backend; there is no installer jargon.
 * - When something needs attention (missing CCP4, requirements not installed),
 *   the setup section opens itself and shows exactly what to fix.
 * - No teams / auth / accounts UI: this is a single-user desktop tool.
 */
// Grace period (seconds) before auto-entering CCP4i2 once setup is complete.
const AUTO_LAUNCH_SECONDS = 5;

export const ConfigContent: React.FC = () => {
  const { setTheme } = useTheme();
  const [config, setConfig] = useState<any | null>(null);
  const router = useRouter();
  const { devMode, setDevMode } = useCCP4i2Window();
  const [existingFiles, setExistingFiles] = useState<any | null>(null);
  const [requirementsExist, setRequirementsExist] = useState<boolean>(false);
  const [serverVersion, setServerVersion] = useState<string | null>(null);
  const { setMessage } = usePopcorn();
  const [launching, setLaunching] = useState(false);
  const [setupOpen, setSetupOpen] = useState(false);
  // Auto-launch countdown: when setup is complete, count down then enter CCP4i2
  // automatically, unless the user intercepts. `countdown` is null when idle.
  const [countdown, setCountdown] = useState<number | null>(null);
  // Set when the user cancels this launch's countdown, so it doesn't re-arm.
  const [autoLaunchCancelled, setAutoLaunchCancelled] = useState(false);
  const [installProgress, setInstallProgress] = useState<{
    isInstalling: boolean;
    output: string[];
    status: "idle" | "started" | "installing" | "completed" | "failed";
  }>({
    isInstalling: false,
    output: [],
    status: "idle",
  });

  useEffect(() => {
    if (typeof window !== "undefined" && window.electronAPI) {
      window.electronAPI.sendMessage("get-config");
      window.electronAPI.onMessage("message-from-main", messageHandler);
      return () => {
        window.electronAPI.removeMessageListener(
          "message-from-main",
          messageHandler
        );
      };
    }
  }, []);

  const messageHandler = useCallback(
    (event: any, data: any) => {
      if (data.message === "get-config") {
        setConfig(data.config);
        setDevMode(data.config.devMode);
        setTheme(data.config.theme || "dark");
      } else if (data.message === "start-uvicorn") {
        router.push("/ccp4i2");
      } else if (data.message === "check-file-exists") {
        if (config) {
          if (data.path === config.CCP4I2_PROJECTS_DIR) {
            setExistingFiles((prev: any) => ({
              ...prev,
              CCP4I2_PROJECTS_DIR: data.exists,
            }));
          }
          if (data.path === config.CCP4Dir) {
            setExistingFiles((prev: any) => ({ ...prev, CCP4Dir: data.exists }));
          }
          if (data.path === config.venv_python) {
            setExistingFiles((prev: any) => ({
              ...prev,
              venv_python: data.exists,
            }));
          }
        }
      } else if (data.message === "requirements-exist") {
        setRequirementsExist(true);
        setServerVersion(data.version || null);
      } else if (data.message === "requirements-missing") {
        setRequirementsExist(false);
        setServerVersion(null);
        setMessage(data.error || "Requirements are missing");
      } else if (data.message === "install-requirements-progress") {
        setInstallProgress((prev) => {
          const newOutput = data.output
            ? [...prev.output, data.output]
            : prev.output;
          return {
            isInstalling:
              data.status === "started" || data.status === "installing",
            output: newOutput,
            status: data.status,
          };
        });
        if (data.status === "completed") {
          setMessage("Requirements installed successfully");
          setRequirementsExist(true);
        } else if (data.status === "failed") {
          setMessage("Requirements installation failed");
        }
      }
    },
    [config, setMessage]
  );

  useEffect(() => {
    if (config && typeof window !== "undefined" && window.electronAPI) {
      window.electronAPI.removeMessageListener(
        "message-from-main",
        messageHandler
      );
      window.electronAPI.onMessage("message-from-main", messageHandler);
      window.electronAPI.sendMessage("check-file-exists", {
        path: config.CCP4I2_PROJECTS_DIR,
      });
      window.electronAPI.sendMessage("check-file-exists", {
        path: config.CCP4Dir,
      });
      window.electronAPI.sendMessage("check-file-exists", {
        path: config.venv_python,
      });
      window.electronAPI.sendMessage("check-requirements", {
        path: config.venv_python,
      });
      return () => {
        if (typeof window !== "undefined" && window.electronAPI) {
          window.electronAPI.removeMessageListener(
            "message-from-main",
            messageHandler
          );
        }
      };
    }
  }, [config]);

  const onLocateCcp4 = () => window.electronAPI?.sendMessage("locate-ccp4");
  const onSelectProjectsDir = () =>
    window.electronAPI?.sendMessage("locate-ccp4i2-project-directory");

  const onEnter = () => {
    if (config && window.electronAPI) {
      setLaunching(true);
      window.electronAPI.sendMessage("start-uvicorn", {
        ...config,
        CCP4Dir: config.CCP4Dir.path,
      });
    } else {
      // Web mode: the server is already running server-side.
      router.push("/ccp4i2");
    }
  };

  // Intercept the countdown: stay on this screen and turn auto-launch off for
  // good (the Setup switch can turn it back on). Enter now stays available.
  const onCancelAutoLaunch = () => {
    setAutoLaunchCancelled(true);
    setCountdown(null);
    window.electronAPI?.sendMessage("set-auto-launch", { enabled: false });
  };

  // Setup switch: persist the preference and, when re-enabling, allow the
  // countdown to re-arm this session.
  const onToggleAutoLaunch = (ev: React.ChangeEvent<HTMLInputElement>) => {
    const enabled = ev.target.checked;
    if (enabled) setAutoLaunchCancelled(false);
    window.electronAPI?.sendMessage("set-auto-launch", { enabled });
  };

  const onInstallRequirements = () => {
    if (window.electronAPI) {
      // Starting an install intercepts any running countdown for this session
      // (without persisting auto-launch off) so we never enter mid-pip-install.
      setAutoLaunchCancelled(true);
      setCountdown(null);
      setInstallProgress({ isInstalling: true, output: [], status: "started" });
      window.electronAPI.sendMessage("install-requirements", {
        ...config,
        CCP4Dir: config.CCP4Dir.path,
      });
    }
  };

  const onToggleDevMode = (ev: React.ChangeEvent<HTMLInputElement>) => {
    window.electronAPI?.sendMessage("toggle-dev-mode", {});
    ev.preventDefault();
    ev.stopPropagation();
  };

  const onCloseProgressDialog = () =>
    setInstallProgress({ isInstalling: false, status: "idle", output: [] });

  const hasElectron =
    typeof window !== "undefined" && !!window.electronAPI;

  // The exact backend version this app build is pinned to (from the main
  // process). Shown so the user can see what's EXPECTED, not just what's
  // installed — a mismatch (e.g. installed 3.1.0a1 under a 3.1.0a3 app) is then
  // obvious rather than silent.
  const requiredVersion: string | undefined = config?.requiredServerVersion;
  // Unpacked/dev build (npm run start): the backend must be an EDITABLE install
  // of server/, not the pinned wheel — so the version reference doesn't apply.
  // This mirrors the main-process probe's isDev gate (not the devMode toggle).
  const serverIsDev = !!config?.isDev;
  const norm = (v?: string | null) =>
    (v || "").trim().toLowerCase().replace(/[-_]/g, "");
  const versionMismatch =
    !serverIsDev &&
    !!requirementsExist &&
    !!requiredVersion &&
    !!serverVersion &&
    norm(serverVersion) !== norm(requiredVersion);

  // Ready to enter when the venv exists AND the backend actually imports. Dev
  // mode relaxes only the *version match* (handled in the requirements probe —
  // it never version-gates in dev), NOT the existence of the backend: with
  // nothing installed the probe fails and Django can't even start (e.g.
  // ModuleNotFoundError: corsheaders), so "nothing installed" is never ready,
  // dev or not. requirementsExist already encodes "importable, version-relaxed
  // in dev", so gate on it directly.
  const isReady =
    config && hasElectron && existingFiles?.venv_python && requirementsExist;

  // What still needs doing, in the user's words.
  const blockers: string[] = [];
  if (config && hasElectron) {
    if (!existingFiles?.CCP4Dir) blockers.push("Locate your CCP4 installation");
    if (!existingFiles?.venv_python) blockers.push("A Python environment is missing");
    if (!requirementsExist)
      blockers.push(
        serverIsDev
          ? "Install the CCP4i2 backend as an editable install of your server directory"
          : "Install the CCP4i2 backend into your CCP4 environment"
      );
  }

  // Open the setup section automatically when something needs attention.
  useEffect(() => {
    if (config && hasElectron && !isReady) setSetupOpen(true);
  }, [config, hasElectron, isReady]);

  // Auto-launch preference (persisted in the electron store; on by default).
  const autoLaunchPref = config?.autoLaunch !== false;
  // The countdown is "armed" only when setup is complete, we're not already
  // launching or installing, the preference is on, and the user hasn't
  // cancelled this round.
  const autoLaunchArmed =
    !!isReady &&
    hasElectron &&
    !launching &&
    !installProgress.isInstalling &&
    autoLaunchPref &&
    !autoLaunchCancelled;

  // Run the countdown while armed; disarming (cancel, launch, config change)
  // clears it. Uses functional updates so the interval closure stays valid.
  useEffect(() => {
    if (!autoLaunchArmed) {
      setCountdown(null);
      return;
    }
    setCountdown(AUTO_LAUNCH_SECONDS);
    const id = setInterval(() => {
      setCountdown((c) => (c === null ? null : c <= 1 ? 0 : c - 1));
    }, 1000);
    return () => clearInterval(id);
  }, [autoLaunchArmed]);

  // Enter automatically the instant the countdown reaches zero.
  useEffect(() => {
    if (countdown === 0 && autoLaunchArmed) onEnter();
    // onEnter flips `launching`, which disarms and stops any re-fire.
  }, [countdown, autoLaunchArmed]);

  return (
    <Stack spacing={2}>
      {/* Primary action — enter, or tell the user what's blocking it. */}
      {isReady || !config ? (
        <Paper
          variant="outlined"
          sx={{ p: 3, textAlign: "center", borderColor: "primary.main" }}
        >
          <Stack spacing={1.5} alignItems="center">
            <Typography variant="body2" color="text.secondary">
              {serverVersion ? `Ready · ccp4i2 ${serverVersion}` : "Ready to go"}
            </Typography>
            {requiredVersion && (
              <Typography variant="caption" color="text.secondary">
                This app expects ccp4i2 {requiredVersion}
              </Typography>
            )}
            {countdown !== null && (
              <Typography variant="body2" fontWeight={600}>
                Entering CCP4i2 automatically in {countdown}s…
              </Typography>
            )}
            {countdown !== null ? (
              <Stack direction="row" spacing={1.5}>
                <Button
                  variant="contained"
                  size="large"
                  startIcon={
                    launching ? (
                      <CircularProgress size={18} color="inherit" />
                    ) : (
                      <RocketLaunch />
                    )
                  }
                  onClick={onEnter}
                  disabled={launching}
                  sx={{ minWidth: 160 }}
                >
                  {launching ? "Starting…" : "Enter now"}
                </Button>
                <Button
                  variant="outlined"
                  size="large"
                  startIcon={<Cancel />}
                  onClick={onCancelAutoLaunch}
                  disabled={launching}
                >
                  Cancel
                </Button>
              </Stack>
            ) : (
              <Button
                variant="contained"
                size="large"
                startIcon={
                  launching ? (
                    <CircularProgress size={18} color="inherit" />
                  ) : (
                    <RocketLaunch />
                  )
                }
                onClick={onEnter}
                disabled={launching}
                sx={{ minWidth: 220 }}
              >
                {launching ? "Starting…" : "Enter CCP4i2"}
              </Button>
            )}
          </Stack>
        </Paper>
      ) : (
        <Alert severity="info" icon={false}>
          <Typography fontWeight={600} sx={{ mb: 0.5 }}>
            A little setup first
          </Typography>
          <Stack component="ul" sx={{ m: 0, pl: 2.5 }} spacing={0.25}>
            {blockers.map((b) => (
              <li key={b}>
                <Typography variant="body2">{b}</Typography>
              </li>
            ))}
          </Stack>
        </Alert>
      )}

      {/* Setup detail — collapsed when all is well, open when it isn't. */}
      {config && (
        <Accordion
          expanded={setupOpen}
          onChange={(_, v) => setSetupOpen(v)}
          variant="outlined"
          disableGutters
          sx={{ "&:before": { display: "none" }, borderRadius: 1 }}
        >
          <AccordionSummary expandIcon={<ExpandMore />}>
            <Typography variant="subtitle2" fontWeight={600}>
              Setup
            </Typography>
          </AccordionSummary>
          <AccordionDetails>
            <Stack spacing={1}>
              <SetupRow
                ok={existingFiles?.CCP4Dir}
                label="CCP4"
                value={config.CCP4Dir}
                action={hasElectron ? { label: "Change", onClick: onLocateCcp4 } : undefined}
              />
              <SetupRow
                ok={existingFiles?.venv_python}
                label="Python"
                value={config.venv_python || "Not found"}
              />
              <SetupRow
                ok={existingFiles?.CCP4I2_PROJECTS_DIR}
                // Not yet created is the normal first-run state, not a fault:
                // the server creates this directory when it starts. A red cross
                // here told first-time users something was wrong and that they
                // had to go and fix it.
                pendingCreation
                label="Projects"
                value={
                  existingFiles?.CCP4I2_PROJECTS_DIR
                    ? config.CCP4I2_PROJECTS_DIR
                    : `${config.CCP4I2_PROJECTS_DIR} — will be created`
                }
                action={
                  hasElectron
                    ? { label: "Change", onClick: onSelectProjectsDir }
                    : undefined
                }
              />
              <SetupRow
                // A version mismatch is not "OK" — this app runs only its exact
                // partner backend, so surface it as needing action.
                ok={requirementsExist && !versionMismatch}
                label="ccp4i2 backend"
                value={
                  !requirementsExist
                    ? serverIsDev
                      ? "Not installed — needs an editable install of server/"
                      : requiredVersion
                        ? `Not installed — needs ${requiredVersion}`
                        : "Not installed"
                    : versionMismatch
                      ? `Installed ${serverVersion}, but this app needs ${requiredVersion}`
                      : serverVersion
                        ? serverIsDev
                          ? `ccp4i2 ${serverVersion} (editable)`
                          : `ccp4i2 ${serverVersion}`
                        : "Installed"
                }
                action={
                  hasElectron
                    ? {
                        // Name what the button will install so the action is
                        // unambiguous: the exact version in production
                        // (e.g. "Install 3.1.0a3"), an editable install in dev.
                        label: !requirementsExist || versionMismatch
                          ? serverIsDev
                            ? "Install (editable)"
                            : requiredVersion
                              ? `Install ${requiredVersion}`
                              : "Install"
                          : "Reinstall",
                        onClick: onInstallRequirements,
                        variant:
                          !requirementsExist || versionMismatch
                            ? "contained"
                            : "text",
                        disabled: !existingFiles?.venv_python,
                      }
                    : undefined
                }
              />

              <Stack
                direction="row"
                spacing={1}
                alignItems="center"
                justifyContent="space-between"
                sx={{ pt: 1 }}
              >
                <Typography variant="body2">
                  Launch automatically when ready
                </Typography>
                <Switch
                  size="small"
                  checked={autoLaunchPref}
                  onChange={onToggleAutoLaunch}
                  disabled={!hasElectron}
                />
              </Stack>

              <Stack
                direction="row"
                spacing={2}
                alignItems="center"
                justifyContent="space-between"
              >
                <Stack direction="row" spacing={1} alignItems="center">
                  <Typography variant="body2">Theme</Typography>
                  <ThemeToggle />
                </Stack>
                <Stack direction="row" spacing={1} alignItems="center">
                  <Typography variant="body2">Developer mode</Typography>
                  <Switch
                    size="small"
                    checked={devMode}
                    onChange={onToggleDevMode}
                    disabled={!hasElectron}
                  />
                </Stack>
                <Typography variant="caption" color="text.secondary">
                  Ports {config.NEXT_PORT} / {config.UVICORN_PORT}
                </Typography>
              </Stack>
            </Stack>
          </AccordionDetails>
        </Accordion>
      )}

      {/* Web mode (no electron config): just a theme control. */}
      {!config && (
        <Stack direction="row" spacing={1} alignItems="center" justifyContent="center">
          <Typography variant="body2">Theme</Typography>
          <ThemeToggle />
        </Stack>
      )}

      {/* Install progress. */}
      <Dialog
        open={
          installProgress.isInstalling ||
          installProgress.status === "completed" ||
          installProgress.status === "failed"
        }
        onClose={onCloseProgressDialog}
        maxWidth="md"
        fullWidth
      >
        <DialogTitle>
          {installProgress.status === "completed"
            ? "Requirements installed"
            : installProgress.status === "failed"
              ? "Installation failed"
              : "Installing requirements"}
          {installProgress.status === "installing" && (
            <CircularProgress size={20} sx={{ ml: 2 }} />
          )}
        </DialogTitle>
        <DialogContent>
          <Paper
            elevation={0}
            sx={{
              p: 2,
              bgcolor: "grey.900",
              color: "grey.100",
              maxHeight: 400,
              overflowY: "auto",
              fontFamily: "monospace",
              fontSize: "0.875rem",
            }}
          >
            {installProgress.output.length === 0 ? (
              <Typography>Initializing…</Typography>
            ) : (
              installProgress.output.map((line, i) => (
                <Typography
                  key={i}
                  component="pre"
                  sx={{ m: 0, whiteSpace: "pre-wrap", wordBreak: "break-word" }}
                >
                  {line}
                </Typography>
              ))
            )}
          </Paper>
          {(installProgress.status === "completed" ||
            installProgress.status === "failed") && (
            <Button
              onClick={onCloseProgressDialog}
              variant="contained"
              sx={{ mt: 2 }}
              fullWidth
            >
              Close
            </Button>
          )}
        </DialogContent>
      </Dialog>
    </Stack>
  );
};

const SetupRow: React.FC<{
  ok: boolean | undefined;
  label: string;
  value: string;
  /** Absent is normal — the item is created on demand, so show it as pending
   *  rather than as a fault the user has to go and fix. */
  pendingCreation?: boolean;
  action?: {
    label: string;
    onClick: () => void;
    variant?: "text" | "contained" | "outlined";
    disabled?: boolean;
  };
}> = ({ ok, label, value, action, pendingCreation }) => (
  <Stack direction="row" alignItems="center" spacing={1}>
    {ok ? (
      <Check color="success" fontSize="small" />
    ) : pendingCreation ? (
      <RadioButtonUnchecked color="disabled" fontSize="small" />
    ) : (
      <Cancel color="error" fontSize="small" />
    )}
    <Typography variant="body2" sx={{ minWidth: 96, fontWeight: 500 }}>
      {label}
    </Typography>
    <Typography
      variant="body2"
      color="text.secondary"
      sx={{
        flex: 1,
        fontFamily: "monospace",
        fontSize: "0.75rem",
        overflow: "hidden",
        textOverflow: "ellipsis",
        whiteSpace: "nowrap",
      }}
    >
      {value}
    </Typography>
    {action && (
      <Button
        size="small"
        variant={action.variant}
        onClick={action.onClick}
        disabled={action.disabled}
      >
        {action.label}
      </Button>
    )}
  </Stack>
);
