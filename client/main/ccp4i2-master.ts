import { app, BrowserWindow, ipcMain, session } from "electron";
import path from "path";
import fs from "fs";
import Store from "electron-store";
import detectPort from "detect-port";
import { fileURLToPath } from "node:url";
import { startNextServer } from "./ccp4i2-next-server";
import { installIpcHandlers } from "./ccp4i2-ipc";
import { Server } from "node:http";
import { installWillDownloadHandler } from "./ccp4i2-session";
import { StoreSchema } from "../types/store";
import { createWindow } from "./ccp4i2-create-window";
import { addNewWindowMenuItem } from "./ccp4i2-menu";
import { setupZoomLevel } from "./ccp4i2-zoom";
import os from "os";

const isDev = !app.isPackaged; // ✅ Works in compiled builds

// Change the current working directory to the Resources folder
if (!isDev) {
  const asarDir = app.getAppPath();
  const resourcesDir = path.resolve(asarDir, ".."); // Resolve the parent directory of the app.asar file
  // Change the working directory to the Resources folder where .next is
  process.chdir(resourcesDir);
}

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const getProjectsDir = () => {
  const homeDir = os.homedir();
  const isWindows = process.platform === "win32";
  return isWindows
    ? path.join(homeDir, "ccp4x", "CCP4X_PROJECTS")
    : path.join(homeDir, ".ccp4x", "CCP4X_PROJECTS");
};

// Compute projectRoot based on dev vs packaged mode
// This is NOT user-configurable - it's always derived from the app location
// In dev: parent of client/ directory (the ccp4i2 repo root where qticons/, svgicons/ exist)
// In packaged: resources/ccp4i2 (bundled modules and icons)
export const getProjectRoot = () => {
  if (isDev) {
    // In dev mode, client is at ccp4i2/client, so parent is projectRoot
    return path.resolve(process.cwd(), "..");
  }
  // In packaged mode, bundled ccp4i2 is in resources
  return path.join((process as any).resourcesPath as string, "ccp4i2");
};

// Get a sensible default for CCP4Dir by checking common locations
const getDefaultCCP4Dir = () => {
  const isWindows = process.platform === "win32";
  const isMac = process.platform === "darwin";

  // Check common CCP4 installation locations in order of preference
  const possiblePaths: string[] = [];

  if (isDev) {
    // In dev mode, check sibling directory (../ccp4-* patterns)
    const parentDir = path.resolve(process.cwd(), "../..");
    try {
      const siblings = fs.readdirSync(parentDir);
      // Sort descending to prefer newer versions (ccp4-20251105 > ccp4-9)
      const ccp4Dirs = siblings
        .filter((name) => name.startsWith("ccp4"))
        .sort()
        .reverse();
      for (const dir of ccp4Dirs) {
        possiblePaths.push(path.join(parentDir, dir));
      }
    } catch {
      // Ignore errors reading directory
    }
  }

  // Standard installation locations
  if (isMac) {
    possiblePaths.push("/Applications/ccp4-9");
    possiblePaths.push("/Applications/ccp4-8");
  } else if (isWindows) {
    possiblePaths.push("C:\\CCP4\\ccp4-9");
    possiblePaths.push("C:\\CCP4\\ccp4-8");
  } else {
    // Linux
    possiblePaths.push("/opt/ccp4");
    possiblePaths.push("/usr/local/ccp4");
  }

  // Return first path that exists with ccp4-python
  for (const ccp4Path of possiblePaths) {
    const pythonBin = isWindows ? "ccp4-python.bat" : "ccp4-python";
    const pythonPath = path.join(ccp4Path, "bin", pythonBin);
    if (fs.existsSync(pythonPath)) {
      return ccp4Path;
    }
  }

  // Fallback to first standard location (user will need to configure)
  return isMac ? "/Applications/ccp4-9x" : isWindows ? "C:\\CCP4\\ccp4-9" : "/opt/ccp4";
};

export const store = new Store<StoreSchema>({
  defaults: {
    CCP4Dir: getDefaultCCP4Dir(),
    projectRoot: getProjectRoot(), // Computed, not user-configurable
    devMode: false,
    zoomLevel: -2,
    CCP4I2_PROJECTS_DIR: getProjectsDir(),
    theme: "dark",
  },
});

// Note: projectRoot is always computed via getProjectRoot(), not read from store
// The store default is just for schema compatibility

let mainWindow: BrowserWindow | null = null;
let nextServerPort: number | null = null;
let nextServer: Server | null = null;
let djangoServerPort: number | null = null;
let djangoServer: any | null = null;

const setDjangoServer = (server) => {
  if (djangoServer) {
    djangoServer.kill();
  }
  djangoServer = server;
};

const getMainWindow = () => {
  if (mainWindow) {
    return mainWindow;
  } else {
    console.error("getMainWindow: Main window is not available");
    return null;
  }
};

app
  .whenReady()
  .then(async () => {
    nextServerPort = await detectPort(3000);
    djangoServerPort = await detectPort(nextServerPort + 1);
    installIpcHandlers(
      ipcMain,
      getMainWindow,
      store,
      djangoServerPort,
      nextServerPort,
      isDev,
      setDjangoServer
    );
    installWillDownloadHandler(session.defaultSession);
    addNewWindowMenuItem(nextServerPort, djangoServerPort);
    setupZoomLevel(store);
    // Set both API_BASE_URL (runtime, for server-side proxy routes) and
    // NEXT_PUBLIC_API_BASE_URL (for any client-side code that needs it)
    process.env.API_BASE_URL = `http://localhost:${djangoServerPort}`;
    process.env.NEXT_PUBLIC_API_BASE_URL = `http://localhost:${djangoServerPort}`;
    nextServer = await startNextServer(isDev, nextServerPort, djangoServerPort);
  })
  .then(async () => {
    mainWindow = await createWindow(
      `http://localhost:${nextServerPort}/config`,
      store
    );
    console.log({
      CCP4Dir: store.get("CCP4Dir"),
      djangoServerPort,
      nextServerPort,
    });
  });

app.on("window-all-closed", () => {
  if (process.platform !== "darwin") {
    app.quit();
  }
});

app.on("before-quit", () => {
  nextServer?.close();
  djangoServer?.kill();
});

app.on("activate", async () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    mainWindow = await createWindow(
      `http://localhost:${nextServerPort}`,
      store
    );
  }
});
