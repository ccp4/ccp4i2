import { BrowserWindow } from "electron";
import { startDjangoServer } from "./ccp4i2-django-server";
import { platform } from "node:os";
import Store from "electron-store";
import { dialog } from "electron";
import path from "node:path";
import fs from "node:fs";
import { spawn, ChildProcessWithoutNullStreams } from "node:child_process";
import { fileURLToPath } from "node:url";
import { StoreSchema } from "../types/store";
import { getProjectRoot } from "./ccp4i2-master";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

/**
 * Finds the Python executable in the project's virtual environment.
 * Checks multiple possible venv locations in order of preference.
 *
 * @param projectRoot - The path to the ccp4i2 project root.
 * @returns The path to the Python executable, or null if not found.
 */
function findVenvPython(projectRoot: string): string | null {
  const isWindows = platform() === "win32";
  const pythonBin = isWindows ? "python.exe" : "python";
  const binDir = isWindows ? "Scripts" : "bin";

  // Check possible venv locations in order of preference
  const venvPaths = [
    path.join(projectRoot, ".venv", binDir, pythonBin),
    path.join(projectRoot, ".venv.py311", binDir, pythonBin),
    path.join(projectRoot, ".venv.py39", binDir, pythonBin),
  ];

  for (const venvPath of venvPaths) {
    if (fs.existsSync(venvPath)) {
      return venvPath;
    }
  }

  return null;
}

/**
 * Finds the Python executable, preferring ccp4-python from the CCP4 installation.
 * Falls back to project's virtual environment if ccp4-python is not available.
 *
 * @param CCP4Dir - The path to the CCP4 installation directory.
 * @param projectRoot - The path to the project root (fallback for .venv).
 * @returns The path to the Python executable, or null if not found.
 */
function findPython(CCP4Dir: string, projectRoot: string): string | null {
  const isWindows = platform() === "win32";

  // Prefer ccp4-python from CCP4 installation
  const ccp4PythonBin = isWindows ? "ccp4-python.bat" : "ccp4-python";
  const ccp4PythonPath = path.join(CCP4Dir, "bin", ccp4PythonBin);
  if (fs.existsSync(ccp4PythonPath)) {
    return ccp4PythonPath;
  }

  // Fallback to project's virtual environment
  return findVenvPython(projectRoot);
}

/**
 * Installs IPC handlers for communication between the Electron main process and renderer processes.
 *
 * @param ipcMain - The Electron IpcMain instance used to listen for IPC events.
 * @param getMainWindow - A function that returns the main BrowserWindow instance or null if unavailable.
 * @param store - An instance of the Electron Store used to persist application state.
 * @param djangoServerPort - The port number for the Django server.
 * @param nextServerPort - The port number for the Next.js server.
 * @param isDev - A boolean indicating whether the application is running in development mode.
 * @param setDjangoServer - A function to set the Django server process instance.
 *
 * @remarks
 * This function sets up various IPC event listeners to handle tasks such as:
 * - Locating the CCP4 directory.
 * - Checking if a file exists.
 * - Selecting a directory for CCP4I2 projects.
 * - Starting the Uvicorn (Django) server.
 * - Retrieving the current configuration.
 * - Toggling developer mode.
 * - Adjusting zoom levels (zoom in, zoom out, reset).
 *
 * Each IPC event listener performs specific actions and replies to the renderer process with the results.
 */
export const installIpcHandlers = (
  ipcMain: Electron.IpcMain,
  getMainWindow: () => BrowserWindow | null,
  store: Store<StoreSchema>,
  djangoServerPort: number,
  nextServerPort: number,
  isDev: boolean,
  setDjangoServer: (server: ChildProcessWithoutNullStreams) => void
) => {
  const getConfigResponse = () => {
    const projectRoot = getProjectRoot(); // Always computed, not from store
    const CCP4Dir = store.get("CCP4Dir") || "";
    // Prefer ccp4-python from CCP4 installation, fall back to venv
    const python_path = findPython(CCP4Dir, projectRoot);
    const config: any = store.store; // Retrieve all values from the electron-store
    // Override projectRoot with computed value (not user-configurable)
    config.projectRoot = projectRoot;
    // Use venv_python key for backwards compatibility with the UI
    config.venv_python = python_path;
    config.UVICORN_PORT = djangoServerPort;
    config.NEXT_PORT = nextServerPort;
    return {
      message: "get-config",
      status: "Success",
      config,
    };
  };

  const getCwdResponse = () => {
    return {
      message: "get-cwd",
      status: "Success",
      cwd: isDev
        ? path.join(process.cwd(), "..", "server")
        : "", // In packaged mode, ccp4i2 is pip-installed — no server directory on disk
    };
  };

  // IPC communication to trigger file dialog to locate a valid CCP4 directory
  ipcMain.on("locate-ccp4", (event, _data) => {
    const mainWindow: BrowserWindow | null = getMainWindow();
    if (!mainWindow) {
      console.error("Main window is not available");
      return;
    }

    dialog
      .showOpenDialog(mainWindow, {
        properties: ["openDirectory"],
      })
      .then((result) => {
        if (!result.canceled) {
          console.log("Selected CCP4 directory:", result.filePaths);
          store.set("CCP4Dir", result.filePaths[0]);
          event.reply("message-from-main", getConfigResponse());
          process.env.CCP4 = result.filePaths[0];
        }
      });
  });

  // IPC communication to trigger file dialog to locate the ccp4i2 project root
  ipcMain.on("locate-project-root", (event, _data) => {
    const mainWindow: BrowserWindow | null = getMainWindow();
    if (!mainWindow) {
      console.error("Main window is not available");
      return;
    }

    dialog
      .showOpenDialog(mainWindow, {
        properties: ["openDirectory"],
        title: "Select CCP4i2 Project Directory",
        message: "Select the root directory of the ccp4i2 project (containing .venv)",
      })
      .then((result) => {
        if (!result.canceled) {
          const selectedPath = result.filePaths[0];
          // Validate that this looks like a valid project root
          const venvPython = findVenvPython(selectedPath);
          if (venvPython) {
            console.log("Selected project root:", selectedPath);
            store.set("projectRoot", selectedPath);
            event.reply("message-from-main", getConfigResponse());
          } else {
            event.reply("message-from-main", {
              message: "locate-project-root",
              status: "Error",
              error: `No Python virtual environment found in ${selectedPath}. Please select a directory containing .venv`,
            });
          }
        }
      });
  });

  // IPC communication to trigger file dialog to select parent directory for new projects
  // and set the CCP4I2_PROJECTS_DIR in the store
  ipcMain.on("check-file-exists", (event, data) => {
    event.reply("message-from-main", {
      message: "check-file-exists",
      path: data.path,
      exists: fs.existsSync(data.path),
    });
  });

  // IPC communication to trigger file dialog to select parent directory for new projects
  // and set the CCP4I2_PROJECTS_DIR in the store
  ipcMain.on("locate-ccp4i2-project-directory", (event, data) => {
    const mainWindow: BrowserWindow | null = getMainWindow();
    if (!mainWindow) {
      console.error("Main window is not available");
      return;
    }
    dialog
      .showOpenDialog(mainWindow, {
        properties: ["openDirectory", "createDirectory"],
      })
      .then((result) => {
        if (!result.canceled) {
          console.log("Selected directory:", result.filePaths);
          store.set("CCP4I2_PROJECTS_DIR", result.filePaths[0]);
          event.reply("message-from-main", getConfigResponse());
        }
      });
  });

  // IPC communication to trigger file dialog to start the uvicorn (django) server
  ipcMain.on("start-uvicorn", async (event, _data) => {
    if (!djangoServerPort) return;
    if (!nextServerPort) return;
    const djangoServer = await startDjangoServer(
      store.get("CCP4Dir"),
      getProjectRoot(), // Always computed, not from store
      djangoServerPort,
      nextServerPort,
      isDev,
      store.get("CCP4I2_PROJECTS_DIR")
    );
    setDjangoServer(djangoServer);
    event.reply("message-from-main", {
      message: "start-uvicorn",
      status: "Success",
    });
  });

  // IPC communication to prompt reply with current config response
  ipcMain.on("get-config", (event, _data) => {
    event.reply("message-from-main", getConfigResponse());
  });

  // IPC communication to prompt reply with current config response
  ipcMain.on("get-cwd", (event, _data) => {
    event.reply("message-from-main", getCwdResponse());
  });

  // IPC communication to trigger file dialog to toggle the state of the developer mode
  ipcMain.on("toggle-dev-mode", (event, data) => {
    store.set("devMode", !store.get("devMode"));
    event.reply("message-from-main", getConfigResponse());
  });

  // IPC communication to set theme mode
  ipcMain.on("set-theme-mode", (event, data) => {
    if (data.theme !== "light" && data.theme !== "dark") {
      console.error("Invalid theme mode:", data.theme);
      return;
    }
    store.set("theme", data.theme);
    event.reply("message-from-main", {
      message: "set-theme-mode",
      status: "Success",
      theme: data.theme,
    });
  });

  // Find-in-page: relay to focused window's webContents
  ipcMain.on("find-in-page", (event, data) => {
    const win = BrowserWindow.getFocusedWindow();
    if (!win) return;
    const { text, forward = true, findNext = false } = data;
    if (!text) return;
    win.webContents.findInPage(text, { forward, findNext });
  });

  ipcMain.on("stop-find-in-page", (event, _data) => {
    const win = BrowserWindow.getFocusedWindow();
    if (!win) return;
    win.webContents.stopFindInPage("clearSelection");
  });

  // Save file to user-chosen location via native save dialog.
  // Used by "Save to..." context menu items for file export.
  ipcMain.on("save-file-as", async (event, data) => {
    const { filePath, fileName } = data;
    if (!filePath || !fs.existsSync(filePath)) {
      console.error("save-file-as: source not found:", filePath);
      return;
    }
    const win = BrowserWindow.getFocusedWindow();
    if (!win) return;

    const result = await dialog.showSaveDialog(win, {
      defaultPath: fileName,
    });
    if (!result.canceled && result.filePath) {
      try {
        fs.copyFileSync(filePath, result.filePath);
      } catch (err) {
        console.error("save-file-as: copy failed:", err);
      }
    }
  });

  ipcMain.on("zoom-in", (event, _data) => {
    BrowserWindow.getAllWindows().forEach((win) => {
      const current = win.webContents.getZoomLevel();
      win.webContents.setZoomLevel(current + 1);
      store.set("zoomLevel", current + 1);
    });
  });

  ipcMain.on("zoom-out", (event, _data) => {
    BrowserWindow.getAllWindows().forEach((win) => {
      const current = win.webContents.getZoomLevel();
      win.webContents.setZoomLevel(current - 1);
      store.set("zoomLevel", current - 1);
    });
  });

  ipcMain.on("zoom-reset", (event, _data) => {
    BrowserWindow.getAllWindows().forEach((win) => {
      win.webContents.setZoomLevel(0);
      store.set("zoomLevel", 0);
    });
  });

  ipcMain.on("check-requirements", (event, _data) => {
    const projectRoot = store.get("projectRoot") || "";
    const CCP4Dir = store.get("CCP4Dir") || "";
    const pythonPath = findPython(CCP4Dir, projectRoot);


    // Validate that the executable exists before spawning
    if (!pythonPath) {
      event.reply("message-from-main", {
        message: "requirements-missing",
        error: `Python not found. Please configure CCP4 installation or project virtual environment.`,
      });
      return;
    }

    let errorOutput = "";

    // Add error handling for spawn
    try {
      const child = spawn(pythonPath, ["-c", "import rest_framework"], {
        stdio: ["ignore", "ignore", "pipe"],
        // Add shell option for Windows compatibility
        shell: process.platform === "win32",
      });

      child.stderr?.on("data", (data: Buffer) => {
        errorOutput += data.toString();
      });

      child.on("exit", (code: number) => {
        if (code === 0) {
          event.reply("message-from-main", { message: "requirements-exist" });
        } else {
          event.reply("message-from-main", {
            message: "requirements-missing",
            error: errorOutput.trim() || `Process exited with code ${code}`,
          });
        }
      });

      child.on("error", (error: Error) => {
        console.error("Spawn error:", error);
        event.reply("message-from-main", {
          message: "requirements-missing",
          error: `Failed to execute: ${error.message}`,
        });
      });
    } catch (error) {
      console.error("Failed to spawn process:", error);
      event.reply("message-from-main", {
        message: "requirements-missing",
        error: `Spawn failed: ${error instanceof Error ? error.message : String(error)}`,
      });
    }
  });

  ipcMain.on("install-requirements", (event, _config) => {
    // In packaged mode, dependencies are managed by the CCP4 installation (pip install)
    if (!isDev) {
      event.sender.send("message-from-main", {
        message: "install-requirements-progress",
        status: "failed",
        output: "Dependencies are managed by the CCP4 installation. Please reinstall or update CCP4.",
      });
      return;
    }

    const projectRoot = store.get("projectRoot") || "";
    const CCP4Dir = store.get("CCP4Dir") || "";
    const pythonPath = findPython(CCP4Dir, projectRoot);

    if (!pythonPath) {
      event.sender.send("message-from-main", {
        message: "install-requirements-progress",
        status: "failed",
        output: `Python not found. Please configure CCP4 installation or project virtual environment.`,
      });
      return;
    }

    // In dev mode, install from pyproject.toml
    const serverPath = path.join(process.cwd(), "..", "server");
    const requirementsPath = path.join(serverPath, "requirements.txt");

    // Spawn pip install process
    const pipProcess = spawn(pythonPath, [
      "-m",
      "pip",
      "install",
      "-r",
      requirementsPath,
      "--verbose", // For more detailed output
    ]);

    // Send start message
    event.sender.send("message-from-main", {
      message: "install-requirements-progress",
      status: "started",
    });

    // Capture stdout
    pipProcess.stdout.on("data", (data) => {
      const output = data.toString();
      event.sender.send("message-from-main", {
        message: "install-requirements-progress",
        status: "installing",
        output: output,
      });
    });

    // Capture stderr (pip sends progress info here too)
    pipProcess.stderr.on("data", (data) => {
      const output = data.toString();
      event.sender.send("message-from-main", {
        message: "install-requirements-progress",
        status: "installing",
        output: output,
      });
    });

    // Handle completion
    pipProcess.on("close", (code) => {
      if (code === 0) {
        event.sender.send("message-from-main", {
          message: "install-requirements-progress",
          status: "completed",
          output: "All requirements installed successfully",
        });
        // Recheck requirements
        event.sender.send("message-from-main", {
          message: "requirements-exist",
        });
      } else {
        event.sender.send("message-from-main", {
          message: "install-requirements-progress",
          status: "failed",
          output: `Installation failed with code ${code}`,
        });
      }
    });

    // Handle errors
    pipProcess.on("error", (error) => {
      event.sender.send("message-from-main", {
        message: "install-requirements-progress",
        status: "failed",
        output: `Error: ${error.message}`,
      });
    });
  });
};
