"use client";
import React, { useCallback } from "react";
import {
  Box,
  Button,
  Dialog,
  DialogContent,
  DialogTitle,
  FormControlLabel,
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
  const { setTheme, mode, customColors } = useTheme();
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

  return (
    <Stack spacing={3}>
      {/* Appearance Section - Available regardless of config */}
      <Stack spacing={2}>
        <Typography variant="h6" color="primary" fontWeight={600}>
          Appearance
        </Typography>

        <Paper variant="outlined" sx={{ p: 2 }}>
          <Stack
            direction="row"
            alignItems="center"
            justifyContent="space-between"
          >
            <Box>
              <Typography variant="subtitle1" fontWeight={500}>
                Theme
              </Typography>
              <Typography variant="body2" color="text.secondary">
                Choose between light and dark theme
              </Typography>
            </Box>
            <ThemeToggle />
          </Stack>
        </Paper>
      </Stack>

      {/* Launch Section - Available regardless of config */}
      <Stack spacing={2}>
        <Typography variant="h6" color="primary" fontWeight={600}>
          Launch Application
        </Typography>

        <Paper variant="outlined" sx={{ p: 3, textAlign: "center" }}>
          <Stack spacing={2} alignItems="center">
            <Typography variant="body1" color="text.secondary">
              {config
                ? typeof window !== "undefined" && window.electronAPI
                  ? existingFiles?.venv_python && (devMode || requirementsExist)
                    ? "Ready to start CCP4i2 with current configuration"
                    : "Configure paths and install requirements to launch CCP4i2"
                  : "CCP4i2 configuration - running in web mode"
                : "Access CCP4i2 web interface"}
            </Typography>
            <Button
              variant="contained"
              size="large"
              startIcon={<Folder />}
              onClick={config ? onStartUvicorn : () => router.push("/")}
              disabled={
                config
                  ? !(typeof window !== "undefined" && window.electronAPI) ||
                    !existingFiles?.venv_python ||
                    (!devMode && !requirementsExist)
                  : false
              }
              sx={{
                minWidth: 300,
                bgcolor: config
                  ? typeof window !== "undefined" &&
                    window.electronAPI &&
                    existingFiles?.venv_python &&
                    (devMode || requirementsExist)
                    ? green[500]
                    : undefined
                  : green[500], // Always green for web launch
                "&:hover": {
                  bgcolor: config
                    ? typeof window !== "undefined" &&
                      window.electronAPI &&
                      existingFiles?.venv_python &&
                      (devMode || requirementsExist)
                      ? green[600]
                      : undefined
                    : green[600], // Always green hover for web launch
                },
                "&.Mui-disabled": {
                  bgcolor: "action.disabledBackground",
                  color: "action.disabled",
                },
              }}
            >
              {config
                ? typeof window !== "undefined" && window.electronAPI
                  ? "Launch CCP4i2"
                  : "Web Mode - Launch Unavailable"
                : "Launch CCP4i2"}
            </Button>
            {config &&
              (!existingFiles?.venv_python ||
                (!devMode && !requirementsExist) ||
                !(typeof window !== "undefined" && window.electronAPI)) && (
                <Typography variant="caption" color="error">
                  {typeof window !== "undefined" && window.electronAPI
                    ? "Please configure paths and install requirements first"
                    : "CCP4i2 launch requires Electron environment"}
                </Typography>
              )}
          </Stack>
        </Paper>
      </Stack>

      {config && (
        <>
          {/* File Paths Section */}
          <Stack spacing={2}>
            <Typography variant="h6" color="primary" fontWeight={600}>
              File Paths
            </Typography>

            {/* CCP4Dir */}
            <Paper variant="outlined" sx={{ p: 2 }}>
              <Stack spacing={2}>
                <Box display="flex" alignItems="center" gap={1}>
                  <Typography variant="subtitle1" fontWeight={500}>
                    CCP4 Installation Directory
                  </Typography>
                  {existingFiles?.CCP4Dir ? (
                    <Check color="success" />
                  ) : (
                    <Cancel color="error" />
                  )}
                </Box>
                <Typography
                  variant="body2"
                  color="text.secondary"
                  sx={{
                    fontFamily: "monospace",
                    bgcolor: mode === "dark" ? "grey.800" : "grey.50",
                    p: 1,
                    borderRadius: 1,
                  }}
                >
                  {config.CCP4Dir}
                </Typography>
                <Button
                  variant="outlined"
                  startIcon={<Folder />}
                  onClick={onLaunchBrowser}
                  disabled={
                    !(typeof window !== "undefined" && window.electronAPI)
                  }
                  sx={{ alignSelf: "flex-start" }}
                >
                  Select Directory...
                </Button>
              </Stack>
            </Paper>

            {/* Python Interpreter (ccp4-python or venv) */}
            <Paper variant="outlined" sx={{ p: 2 }}>
              <Stack spacing={2}>
                <Box display="flex" alignItems="center" gap={1}>
                  <Typography variant="subtitle1" fontWeight={500}>
                    Python Interpreter
                  </Typography>
                  {existingFiles?.venv_python ? (
                    <Check color="success" />
                  ) : (
                    <Cancel color="error" />
                  )}
                </Box>
                <Typography
                  variant="body2"
                  color="text.secondary"
                  sx={{
                    fontFamily: "monospace",
                    bgcolor: mode === "dark" ? "grey.800" : "grey.50",
                    p: 1,
                    borderRadius: 1,
                  }}
                >
                  {config.venv_python || "Not found - check CCP4 installation"}
                </Typography>
              </Stack>
            </Paper>

            {/* Projects Directory */}
            <Paper variant="outlined" sx={{ p: 2 }}>
              <Stack spacing={2}>
                <Box display="flex" alignItems="center" gap={1}>
                  <Typography variant="subtitle1" fontWeight={500}>
                    CCP4i2 Projects Directory
                  </Typography>
                  {existingFiles?.CCP4I2_PROJECTS_DIR ? (
                    <Check color="success" />
                  ) : (
                    <Cancel color="error" />
                  )}
                </Box>
                <Typography
                  variant="body2"
                  color="text.secondary"
                  sx={{
                    fontFamily: "monospace",
                    bgcolor: mode === "dark" ? "grey.800" : "grey.50",
                    p: 1,
                    borderRadius: 1,
                  }}
                >
                  {config.CCP4I2_PROJECTS_DIR}
                </Typography>
                <Button
                  variant="outlined"
                  startIcon={<Folder />}
                  onClick={onSelectProjectsDir}
                  disabled={
                    !(typeof window !== "undefined" && window.electronAPI)
                  }
                  sx={{ alignSelf: "flex-start" }}
                >
                  Select Directory...
                </Button>
              </Stack>
            </Paper>
          </Stack>

          {/* Server Configuration Section */}
          <Stack spacing={2}>
            <Typography variant="h6" color="primary" fontWeight={600}>
              Server Configuration
            </Typography>

            <Box display="flex" gap={2} flexWrap="wrap">
              <Paper variant="outlined" sx={{ p: 2, flex: 1, minWidth: 200 }}>
                <Stack spacing={1}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Next.js Port
                  </Typography>
                  <Typography variant="h4" color="primary" fontWeight={600}>
                    {config.NEXT_PORT}
                  </Typography>
                </Stack>
              </Paper>

              <Paper variant="outlined" sx={{ p: 2, flex: 1, minWidth: 200 }}>
                <Stack spacing={1}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Uvicorn Port
                  </Typography>
                  <Typography variant="h4" color="primary" fontWeight={600}>
                    {config.UVICORN_PORT}
                  </Typography>
                </Stack>
              </Paper>
            </Box>
          </Stack>

          {/* System Setup Section */}
          <Stack spacing={2}>
            <Typography variant="h6" color="primary" fontWeight={600}>
              System Setup
            </Typography>

            {/* Requirements */}
            <Paper variant="outlined" sx={{ p: 2 }}>
              <Stack
                direction="row"
                alignItems="center"
                justifyContent="space-between"
              >
                <Box display="flex" alignItems="center" gap={2}>
                  <Box>
                    <Typography variant="subtitle1" fontWeight={500}>
                      Python Requirements
                    </Typography>
                    <Typography variant="body2" color="text.secondary">
                      Required packages for CCP4i2 operation
                    </Typography>
                  </Box>
                  {requirementsExist ? (
                    <Check color="success" />
                  ) : (
                    <Cancel color="error" />
                  )}
                </Box>
                <Button
                  variant="contained"
                  onClick={onInstallRequirements}
                  disabled={
                    !(typeof window !== "undefined" && window.electronAPI) ||
                    !existingFiles?.venv_python
                  }
                >
                  {requirementsExist ? "Reinstall" : "Install"}
                </Button>
              </Stack>
            </Paper>

            {/* Dev Mode */}
            <Paper variant="outlined" sx={{ p: 2 }}>
              <Stack
                direction="row"
                alignItems="center"
                justifyContent="space-between"
              >
                <Box>
                  <Typography variant="subtitle1" fontWeight={500}>
                    Development Mode
                  </Typography>
                  <Typography variant="body2" color="text.secondary">
                    Enable additional debugging features
                  </Typography>
                </Box>
                <FormControlLabel
                  control={
                    <Switch
                      checked={devMode}
                      onChange={onToggleDevMode}
                      name="devModeToggle"
                      color="primary"
                      disabled={
                        !(typeof window !== "undefined" && window.electronAPI)
                      }
                    />
                  }
                  label={devMode ? "Enabled" : "Disabled"}
                />
              </Stack>
            </Paper>
          </Stack>
        </>
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
