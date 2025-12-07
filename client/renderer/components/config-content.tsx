"use client";
import React, { useCallback } from "react";
import {
  Button,
  Dialog,
  DialogContent,
  DialogTitle,
  Stack,
  Switch,
  Typography,
  Paper,
  CircularProgress,
} from "@mui/material";
import { green } from "@mui/material/colors";
import { useApi } from "../api";
import { Cancel, Check, Folder } from "@mui/icons-material";
import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import { useCCP4i2Window } from "../app-context";
import { usePopcorn } from "../providers/popcorn-provider";
import { ThemeToggle } from "./theme-toggle";
import { useTheme } from "../theme/theme-provider";

export const ConfigContent: React.FC = () => {
  const api = useApi();
  const { setTheme } = useTheme();
  const [config, setConfig] = useState<any | null>(null);
  const router = useRouter();
  const { devMode, setDevMode } = useCCP4i2Window();
  const [existingFiles, setExistingFiles] = useState<any | null>(null);
  const [requirementsExist, setRequirementsExist] = useState<boolean>(false);
  const { setMessage } = usePopcorn();
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
    // Send a message to the main process to get the config
    if (typeof window !== "undefined" && window.electronAPI) {
      window.electronAPI.sendMessage("get-config");
      // Listen for messages from the main process

      window.electronAPI.onMessage("message-from-main", messageHandler);

      return () => {
        window.electronAPI.removeMessageListener(
          "message-from-main",
          messageHandler
        );
      };
    } else {
      console.log("Electron API is not available - running in web mode");
    }
  }, []);

  const messageHandler = useCallback(
    (event: any, data: any) => {
      if (data.message === "get-config") {
        setConfig(data.config);
        setDevMode(data.config.devMode);
        setTheme(data.config.theme || "dark");
      } else if (data.message === "start-uvicorn") {
        router.push("/");
      } else if (data.message === "check-file-exists") {
        if (config) {
          if (data.path === config.CCP4I2_PROJECTS_DIR) {
            setExistingFiles((prevState: any) => ({
              ...prevState,
              CCP4I2_PROJECTS_DIR: data.exists,
            }));
          }
          if (data.path === config.CCP4Dir) {
            setExistingFiles((prevState: any) => ({
              ...prevState,
              CCP4Dir: data.exists,
            }));
          }
          if (data.path === config.venv_python) {
            setExistingFiles((prevState: any) => ({
              ...prevState,
              venv_python: data.exists,
            }));
          }
        }
      } else if (data.message === "requirements-exist") {
        setRequirementsExist(true);
      } else if (data.message === "requirements-missing") {
        setRequirementsExist(false);
        setMessage(data.error || "Requirements are missing");
      }
      // Add new handler for installation progress
      else if (data.message === "install-requirements-progress") {
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

        // Show success message when completed
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
    if (config) {
      if (typeof window !== "undefined" && window.electronAPI) {
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
      }
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

  const onLaunchBrowser = async () => {
    if (typeof window !== "undefined" && window?.electronAPI) {
      console.log("Gonna send locate-ccp4");
      console.log(window.electronAPI);
      window.electronAPI.sendMessage("locate-ccp4");
    } else {
      console.error("Electron API is not available");
    }
  };

  const onSelectProjectsDir = async () => {
    if (typeof window !== "undefined" && window.electronAPI) {
      console.log("Gonna send locate-ccp4");
      console.log(window.electronAPI);
      window.electronAPI.sendMessage("locate-ccp4i2-project-directory");
    } else {
      console.error("Electron API is not available");
    }
  };

  const onStartUvicorn = async () => {
    if (typeof window !== "undefined" && window.electronAPI) {
      window.electronAPI.sendMessage("start-uvicorn", {
        ...config,
        CCP4Dir: config.CCP4Dir.path,
      });
    } else {
      console.error("Electron API is not available");
    }
  };

  const onInstallRequirements = async () => {
    if (typeof window !== "undefined" && window.electronAPI) {
      // Reset progress state
      setInstallProgress({
        isInstalling: true,
        output: [],
        status: "started",
      });

      window.electronAPI.sendMessage("install-requirements", {
        ...config,
        CCP4Dir: config.CCP4Dir.path,
      });
    } else {
      console.error("Electron API is not available");
    }
  };

  const onToggleDevMode = async (
    ev: React.ChangeEvent<HTMLInputElement>
  ): Promise<void> => {
    if (typeof window !== "undefined" && window?.electronAPI) {
      window.electronAPI.sendMessage("toggle-dev-mode", {});
    }
    ev.preventDefault();
    ev.stopPropagation();
  };

  const onCloseProgressDialog = () => {
    setInstallProgress({
      isInstalling: false,
      status: "idle",
      output: [],
    });
  };

  // Check if ready to launch
  const isReadyToLaunch =
    config &&
    typeof window !== "undefined" &&
    window.electronAPI &&
    existingFiles?.venv_python &&
    (devMode || requirementsExist);

  return (
    <Stack spacing={1.5}>
      {/* Launch Button - Prominent at top */}
      <Paper
        variant="outlined"
        sx={{
          p: 2,
          textAlign: "center",
          bgcolor: isReadyToLaunch ? "success.dark" : undefined,
          borderColor: isReadyToLaunch ? "success.main" : undefined,
        }}
      >
        <Stack direction="row" spacing={2} alignItems="center" justifyContent="center">
          <Button
            variant="contained"
            size="large"
            startIcon={<Folder />}
            onClick={config ? onStartUvicorn : () => router.push("/")}
            disabled={config ? !isReadyToLaunch : false}
            sx={{
              minWidth: 200,
              bgcolor: isReadyToLaunch || !config ? green[500] : undefined,
              "&:hover": {
                bgcolor: isReadyToLaunch || !config ? green[600] : undefined,
              },
            }}
          >
            Launch CCP4i2
          </Button>
          {config && !isReadyToLaunch && (
            <Typography variant="caption" color="error">
              Configure paths and install requirements first
            </Typography>
          )}
        </Stack>
      </Paper>

      {config && (
        <>
          {/* File Paths - Compact table format */}
          <Paper variant="outlined" sx={{ p: 1.5 }}>
            <Typography variant="subtitle2" color="primary" fontWeight={600} sx={{ mb: 1 }}>
              Configuration
            </Typography>
            <Stack spacing={1}>
              {/* CCP4 Directory */}
              <Stack direction="row" alignItems="center" spacing={1}>
                {existingFiles?.CCP4Dir ? (
                  <Check color="success" fontSize="small" />
                ) : (
                  <Cancel color="error" fontSize="small" />
                )}
                <Typography variant="body2" sx={{ minWidth: 100, fontWeight: 500 }}>
                  CCP4 Dir
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
                  {config.CCP4Dir}
                </Typography>
                <Button
                  size="small"
                  onClick={onLaunchBrowser}
                  disabled={!(typeof window !== "undefined" && window.electronAPI)}
                >
                  Change
                </Button>
              </Stack>

              {/* Python */}
              <Stack direction="row" alignItems="center" spacing={1}>
                {existingFiles?.venv_python ? (
                  <Check color="success" fontSize="small" />
                ) : (
                  <Cancel color="error" fontSize="small" />
                )}
                <Typography variant="body2" sx={{ minWidth: 100, fontWeight: 500 }}>
                  Python
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
                  {config.venv_python || "Not found"}
                </Typography>
              </Stack>

              {/* Projects Directory */}
              <Stack direction="row" alignItems="center" spacing={1}>
                {existingFiles?.CCP4I2_PROJECTS_DIR ? (
                  <Check color="success" fontSize="small" />
                ) : (
                  <Cancel color="error" fontSize="small" />
                )}
                <Typography variant="body2" sx={{ minWidth: 100, fontWeight: 500 }}>
                  Projects
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
                  {config.CCP4I2_PROJECTS_DIR}
                </Typography>
                <Button
                  size="small"
                  onClick={onSelectProjectsDir}
                  disabled={!(typeof window !== "undefined" && window.electronAPI)}
                >
                  Change
                </Button>
              </Stack>

              {/* Requirements */}
              <Stack direction="row" alignItems="center" spacing={1}>
                {requirementsExist ? (
                  <Check color="success" fontSize="small" />
                ) : (
                  <Cancel color="error" fontSize="small" />
                )}
                <Typography variant="body2" sx={{ minWidth: 100, fontWeight: 500 }}>
                  Requirements
                </Typography>
                <Typography variant="body2" color="text.secondary" sx={{ flex: 1 }}>
                  {requirementsExist ? "Installed" : "Not installed"}
                </Typography>
                <Button
                  size="small"
                  variant={requirementsExist ? "text" : "contained"}
                  onClick={onInstallRequirements}
                  disabled={
                    !(typeof window !== "undefined" && window.electronAPI) ||
                    !existingFiles?.venv_python
                  }
                >
                  {requirementsExist ? "Reinstall" : "Install"}
                </Button>
              </Stack>
            </Stack>
          </Paper>

          {/* Settings Row - Theme, Dev Mode, Ports */}
          <Stack direction="row" spacing={1.5}>
            {/* Theme */}
            <Paper variant="outlined" sx={{ p: 1.5, flex: 1 }}>
              <Stack direction="row" alignItems="center" justifyContent="space-between">
                <Typography variant="body2" fontWeight={500}>
                  Theme
                </Typography>
                <ThemeToggle />
              </Stack>
            </Paper>

            {/* Dev Mode */}
            <Paper variant="outlined" sx={{ p: 1.5, flex: 1 }}>
              <Stack direction="row" alignItems="center" justifyContent="space-between">
                <Typography variant="body2" fontWeight={500}>
                  Dev Mode
                </Typography>
                <Switch
                  size="small"
                  checked={devMode}
                  onChange={onToggleDevMode}
                  disabled={!(typeof window !== "undefined" && window.electronAPI)}
                />
              </Stack>
            </Paper>

            {/* Ports */}
            <Paper variant="outlined" sx={{ p: 1.5, flex: 1 }}>
              <Stack direction="row" alignItems="center" justifyContent="space-between">
                <Typography variant="body2" fontWeight={500}>
                  Ports
                </Typography>
                <Typography variant="body2" color="text.secondary">
                  {config.NEXT_PORT} / {config.UVICORN_PORT}
                </Typography>
              </Stack>
            </Paper>
          </Stack>
        </>
      )}

      {/* Minimal config when no electron config */}
      {!config && (
        <Paper variant="outlined" sx={{ p: 1.5 }}>
          <Stack direction="row" alignItems="center" justifyContent="space-between">
            <Typography variant="body2" fontWeight={500}>
              Theme
            </Typography>
            <ThemeToggle />
          </Stack>
        </Paper>
      )}

      {/* Installation Progress Dialog */}
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
            ? "Installation Complete"
            : installProgress.status === "failed"
              ? "Installation Failed"
              : "Installing Requirements"}
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
              <Typography>Initializing installation...</Typography>
            ) : (
              installProgress.output.map((line, index) => (
                <Typography
                  key={index}
                  component="pre"
                  sx={{
                    m: 0,
                    whiteSpace: "pre-wrap",
                    wordBreak: "break-word",
                  }}
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
