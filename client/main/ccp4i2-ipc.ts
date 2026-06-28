import { BrowserWindow } from "electron";
import { startDjangoServer } from "./ccp4i2-django-server";
import os, { platform } from "node:os";
import Store from "electron-store";
import { dialog } from "electron";
import path from "node:path";
import fs from "node:fs";
import { spawn, ChildProcessWithoutNullStreams } from "node:child_process";
import { fileURLToPath } from "node:url";
import { StoreSchema } from "../types/store";
import { getProjectRoot } from "./ccp4i2-master";
import { loadPreferences, updatePreferences, sqliteUrl } from "./ccp4i2-preferences";
import {
  CCP4I2_SERVER_VERSION_FLOOR,
  meetsServerVersionFloor,
} from "./ccp4i2-server-version";

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
    // Shared bootstrap keys live in ~/.ccp4i2/preferences.json (the file the
    // server and the i2/i2run CLI also read). Prefer it, fall back to the store.
    const filePrefs = loadPreferences();
    const CCP4Dir = filePrefs.ccp4Dir || store.get("CCP4Dir") || "";
    // Prefer ccp4-python from CCP4 installation, fall back to venv
    const python_path = findPython(CCP4Dir, projectRoot);
    const config: any = store.store; // Retrieve all values from the electron-store
    // Override projectRoot with computed value (not user-configurable)
    config.projectRoot = projectRoot;
    // Use venv_python key for backwards compatibility with the UI
    config.venv_python = python_path;
    config.UVICORN_PORT = djangoServerPort;
    config.NEXT_PORT = nextServerPort;
    // Overlay the shared keys so the GUI reflects what the server/CLI will use
    // (including values set via the file or a future `i2 preferences set`).
    config.CCP4Dir = CCP4Dir;
    if (filePrefs.projectsDir) config.CCP4I2_PROJECTS_DIR = filePrefs.projectsDir;
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
          // Write the shared key so the server and the i2/i2run CLI see it too.
          updatePreferences({ ccp4Dir: result.filePaths[0] });
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
          const projectsDir = result.filePaths[0];
          store.set("CCP4I2_PROJECTS_DIR", projectsDir);
          // Write the shared keys so the server and the i2/i2run CLI agree —
          // including the database location, so all three open the SAME
          // db.sqlite3 inside the projects dir (not a divergent default).
          updatePreferences({
            projectsDir,
            database: sqliteUrl(path.join(projectsDir, "db.sqlite3")),
          });
          event.reply("message-from-main", getConfigResponse());
        }
      });
  });

  // IPC communication to trigger file dialog to start the uvicorn (django) server
  ipcMain.on("start-uvicorn", async (event, _data) => {
    if (!djangoServerPort) return;
    if (!nextServerPort) return;
    // Use the canonical shared values (file first, store fallback) so the server
    // is launched with the same CCP4 / projects dir the GUI shows and the CLI reads.
    const filePrefs = loadPreferences();
    const djangoServer = await startDjangoServer(
      filePrefs.ccp4Dir || store.get("CCP4Dir"),
      getProjectRoot(), // Always computed, not from store
      djangoServerPort,
      nextServerPort,
      isDev,
      filePrefs.projectsDir || store.get("CCP4I2_PROJECTS_DIR")
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

  // Probe whether this python can actually run the django backend, and reply on
  // "message-from-main" with requirements-exist (carrying the version) or
  // requirements-missing (carrying a real diagnostic). Shared by the
  // check-requirements handler and the post-install recheck so the readiness
  // verdict — including the packaged version-floor gate — is computed in exactly
  // one place. `send` is event.reply / event.sender.send (same target).
  const runRequirementsProbe = (send: (payload: any) => void) => {
    const projectRoot = store.get("projectRoot") || "";
    const CCP4Dir = store.get("CCP4Dir") || "";
    const pythonPath = findPython(CCP4Dir, projectRoot);

    // Validate that the executable exists before spawning
    if (!pythonPath) {
      send({
        message: "requirements-missing",
        error: `Python not found. Please configure CCP4 installation or project virtual environment.`,
      });
      return;
    }

    // Probe the *actual* server entrypoint rather than a transitive dependency.
    // `import rest_framework` only proves some DRF is installed somewhere; it
    // says nothing about whether ccp4i2 is present, importable, or which
    // version. Instead we construct the same ASGI application uvicorn loads at
    // launch (ccp4i2.config.asgi:application) — if that succeeds the server
    // will boot — and emit ccp4i2.__version__ so the UI can show / gate on it.
    // Run with the same DJANGO_SETTINGS_MODULE and, in dev, the same server/
    // cwd as startDjangoServer so the probe can't diverge from the real launch.
    //
    // The probe is written to a temp .py file and run as `python <file>` rather
    // than `python -c "<multi-statement string>"`. On Windows the multi-line
    // string (quotes, JSON braces, semicolons) would be mangled by cmd.exe if
    // spawned with shell:true, and ccp4-python is a .bat that historically
    // wanted a shell. A plain file path is a single clean argv that needs no
    // shell on any platform — matching the (working, shell-less) install spawn.
    const probeSource = [
      "import json, ccp4i2",
      "from ccp4i2.config.asgi import application",
      'print(json.dumps({"version": ccp4i2.__version__}))',
      "",
    ].join("\n");

    let probePath: string;
    try {
      probePath = path.join(os.tmpdir(), `ccp4i2-probe-${process.pid}.py`);
      fs.writeFileSync(probePath, probeSource);
    } catch (error) {
      send({
        message: "requirements-missing",
        error: `Could not write probe file: ${
          error instanceof Error ? error.message : String(error)
        }`,
      });
      return;
    }

    const cleanupProbe = () => {
      try {
        fs.unlinkSync(probePath);
      } catch {
        /* best-effort */
      }
    };

    const serverCwd = isDev
      ? path.join(process.cwd(), "..", "server")
      : undefined;

    let stdoutBuf = "";
    let errorOutput = "";

    // Add error handling for spawn
    try {
      const child = spawn(pythonPath, [probePath], {
        stdio: ["ignore", "pipe", "pipe"],
        env: {
          ...process.env,
          DJANGO_SETTINGS_MODULE: "ccp4i2.config.settings",
          MPLBACKEND: "Agg",
        },
        ...(serverCwd && { cwd: serverCwd }),
      });

      child.stdout?.on("data", (data: Buffer) => {
        stdoutBuf += data.toString();
      });

      child.stderr?.on("data", (data: Buffer) => {
        errorOutput += data.toString();
      });

      child.on("exit", (code: number) => {
        cleanupProbe();
        if (code === 0) {
          let version: string | undefined;
          try {
            version = JSON.parse(stdoutBuf.trim()).version;
          } catch {
            // ASGI app built but version line unparseable — still "ready".
          }
          // Packaged app: an importable-but-stale backend is "not ready" so the
          // UI can offer an upgrade. Dev mode (unpacked) uses the local tree and
          // is never floor-gated — whatever the checkout provides is what runs.
          if (!isDev && !meetsServerVersionFloor(version)) {
            send({
              message: "requirements-missing",
              version,
              error:
                `Installed ccp4i2 ${version || "(unknown)"} is older than the ` +
                `required ${CCP4I2_SERVER_VERSION_FLOOR}. Click Install to upgrade.`,
            });
            return;
          }
          send({
            message: "requirements-exist",
            version,
          });
        } else {
          send({
            message: "requirements-missing",
            error: errorOutput.trim() || `Process exited with code ${code}`,
          });
        }
      });

      child.on("error", (error: Error) => {
        cleanupProbe();
        console.error("Spawn error:", error);
        send({
          message: "requirements-missing",
          error: `Failed to execute: ${error.message}`,
        });
      });
    } catch (error) {
      cleanupProbe();
      console.error("Failed to spawn process:", error);
      send({
        message: "requirements-missing",
        error: `Spawn failed: ${error instanceof Error ? error.message : String(error)}`,
      });
    }
  };

  ipcMain.on("check-requirements", (event, _data) => {
    runRequirementsProbe((payload) => event.reply("message-from-main", payload));
  });

  ipcMain.on("install-requirements", (event, _config) => {
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

    // Mode is keyed strictly on app.isPackaged (isDev):
    //   • Packaged → install the django backend from PyPI as a version floor.
    //     `--upgrade "ccp4i2>=<floor>"` brings a stale/absent backend up to at
    //     least the tested minimum while still letting users pick up newer
    //     compatible releases without an app rebuild.
    //   • Dev (unpacked) → editable install of the local checkout so the
    //     developer's tree is what runs. `-e ./server` also pulls the same
    //     declared dependencies (it reads pyproject.toml), so it strictly
    //     supersedes the old `-r requirements.txt` install.
    const pipArgs = isDev
      ? [
          "-m",
          "pip",
          "install",
          "-e",
          path.join(process.cwd(), "..", "server"),
          "--verbose",
        ]
      : [
          "-m",
          "pip",
          "install",
          "--upgrade",
          `ccp4i2>=${CCP4I2_SERVER_VERSION_FLOOR}`,
          "--verbose",
        ];

    // Spawn pip install process
    const pipProcess = spawn(pythonPath, pipArgs);

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

    // Handle completion.
    //
    // Pip's exit code is NOT treated as the authority on success. Hand-rolled
    // CCP4 python environments can carry distributions with corrupt *.dist-info
    // metadata (e.g. meson/scons with a missing `Version:` field); modern pip
    // iterates all installed dists to print its post-install summary and
    // crashes with BadMetadata AFTER the requested package is already fully
    // installed — exiting non-zero on a successful install. The real authority
    // is the probe: can we import and build ccp4i2.config.asgi (and, packaged,
    // meet the version floor)? So on ANY exit we run the probe and let its
    // verdict decide. Pip's own output is surfaced as the diagnostic only if
    // the probe also says the backend isn't usable.
    pipProcess.on("close", (code) => {
      runRequirementsProbe((payload) => {
        const probeSucceeded = payload.message === "requirements-exist";
        if (probeSucceeded) {
          // Backend imports — install effectively succeeded, even if pip exited
          // non-zero on a broken-metadata summary crash.
          event.sender.send("message-from-main", {
            message: "install-requirements-progress",
            status: "completed",
            output:
              code === 0
                ? "ccp4i2 installed successfully"
                : `ccp4i2 installed successfully (pip exited ${code}, likely ` +
                  `a benign metadata-summary error in the CCP4 environment; ` +
                  `the backend imports correctly).`,
          });
        } else {
          // Backend still not usable — a genuine failure. Prefer pip's exit
          // code context, fall back to the probe's diagnostic.
          event.sender.send("message-from-main", {
            message: "install-requirements-progress",
            status: "failed",
            output:
              `Installation did not produce a usable ccp4i2 backend ` +
              `(pip exit ${code}). ${payload.error || ""}`.trim(),
          });
        }
        // Forward the probe verdict itself so the config page updates its
        // requirements/version state (requirements-exist | requirements-missing).
        event.sender.send("message-from-main", payload);
      });
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
