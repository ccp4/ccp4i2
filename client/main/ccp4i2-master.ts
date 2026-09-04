import { app, BrowserWindow, ipcMain, Menu, session } from "electron";
import path from "path";
import fs from "fs";
import Store from "electron-store";
import detectPort from "detect-port";
import { fileURLToPath } from "node:url";
// Side-effect import: generates the per-launch local-session token and
// installs CCP4I2_LOCAL_SESSION_TOKEN + CCP4I2_LOCAL_USER_EMAIL on
// process.env, so they propagate to both the Django child process and
// the Electron preload. Must be imported before any BrowserWindow is
// created or the Django server is spawned.
import "./ccp4i2-local-session";
import { startNextServer } from "./ccp4i2-next-server";
import { installIpcHandlers } from "./ccp4i2-ipc";
import { Server } from "node:http";
import { installWillDownloadHandler } from "./ccp4i2-session";
import { StoreSchema } from "../types/store";
import { ccp4i2Home, defaultProjectsDir } from "./ccp4i2-preferences";
import { createWindow } from "./ccp4i2-create-window";
import { setupZoomLevel } from "./ccp4i2-zoom";
import { assessPython, listCcp4Dirs } from "./ccp4i2-python-suitability";

const isDev = !app.isPackaged; // ✅ Works in compiled builds

// The AppImage intentionally does NOT force --no-sandbox. appendSwitch("no-sandbox")
// runs too late to be reliable — Chromium's zygote reads the sandbox flag before
// this JS executes — so it never dependably took effect in a packaged AppImage,
// yet it robbed capable users of a real Chromium sandbox. So the AppImage now
// behaves as it did in a11: on permissive kernels it sandboxes normally; on
// locked-down kernels (Ubuntu 24.04+) either enable unprivileged user namespaces
//   sudo sysctl -w kernel.apparmor_restrict_unprivileged_userns=0
// or install the .deb, whose postinst setuids chrome-sandbox + ships an AppArmor
// profile — the proper, sandbox-preserving fix. (Restored at Paul Bond's request.)

// Change the current working directory to the Resources folder
if (!isDev) {
  const asarDir = app.getAppPath();
  const resourcesDir = path.resolve(asarDir, ".."); // Resolve the parent directory of the app.asar file
  // Change the working directory to the Resources folder where .next is
  process.chdir(resourcesDir);
}

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Where new projects go. Delegated to the shared resolver rather than computed
// here: this function used to return ~/.ccp4i2/CCP4X_PROJECTS while the Django
// settings said ~/.ccp4i2-django/CCP4X_PROJECTS, so the desktop app and the
// server it launched disagreed about the project store — which is why testers
// reported seeing several directories where they expected one.
const getProjectsDir = () => defaultProjectsDir();

// Compute projectRoot based on dev vs packaged mode
// This is NOT user-configurable - it's always derived from the app location
// In dev: parent of client/ directory (the ccp4i2 repo root)
// In packaged: empty string (ccp4i2 is pip-installed in ccp4-python, no local project root)
export const getProjectRoot = () => {
  if (isDev) {
    // In dev mode, client is at ccp4i2/client, so parent is projectRoot
    return path.resolve(process.cwd(), "..");
  }
  // In packaged mode, ccp4i2 is pip-installed — no bundled resources
  return "";
};

// Get a sensible default for CCP4Dir by checking common locations

const getDefaultCCP4Dir = () => {
  const isWindows = process.platform === "win32";
  const isMac = process.platform === "darwin";

  // Check common CCP4 installation locations in order of preference
  const possiblePaths: string[] = [];

  // Highest priority: an already-sourced CCP4 environment. When the app is
  // launched from a shell that has sourced ccp4.setup-sh, $CCP4 points at the
  // install root (which holds bin/ccp4-python) and ccp4-python is on PATH.
  // Honour that so a user who already configured CCP4 isn't asked to re-select
  // a directory by hand (the common Linux tarball-in-a-custom-location case).
  const pythonBinName = isWindows ? "ccp4-python.bat" : "ccp4-python";
  if (process.env.CCP4) {
    possiblePaths.push(process.env.CCP4);
  }
  for (const dir of (process.env.PATH || "").split(path.delimiter)) {
    // A ccp4-python on PATH lives in <root>/bin, so the root is its parent.
    if (dir && fs.existsSync(path.join(dir, pythonBinName))) {
      possiblePaths.push(path.dirname(dir));
    }
  }

  if (isDev) {
    // In dev mode, check sibling directory (../ccp4-* patterns)
    // Newest first by number: a plain string sort put ccp4-9 ahead of
    // ccp4-20260702, since "9" > "2".
    const parentDir = path.resolve(process.cwd(), "../..");
    for (const dir of listCcp4Dirs(parentDir)) {
      possiblePaths.push(path.join(parentDir, dir));
    }
  }

  // Standard installation roots. Every ccp4-* directory under each root is a
  // candidate, newest first — "ccp4-20260702" before "ccp4-9" — so the
  // numbering scheme decides the order, not a hard-coded name. The old list
  // named /Applications/ccp4-9 outright, which on 2026-09-04 put a CCP4 9 in
  // front of a user who had just cleared their settings for a demo: its
  // Python 3.9 cannot host this backend and its site-packages is root-owned,
  // so "Install" could only fail.
  const roots = isMac
    ? ["/Applications"]
    : isWindows
      ? ["C:\\CCP4"]
      : ["/opt", "/usr/local"];
  for (const root of roots) {
    for (const name of listCcp4Dirs(root)) possiblePaths.push(path.join(root, name));
  }
  if (!isMac && !isWindows) {
    possiblePaths.push("/opt/ccp4");
    possiblePaths.push("/usr/local/ccp4");
  }

  // Return the first path that exists with a ccp4-python this backend can run
  // on. A CCP4 whose Python is below the floor is skipped, never offered: the
  // setup page asking for a location is better than confidently naming one
  // that cannot work.
  const seen = new Set<string>();
  for (const ccp4Path of possiblePaths) {
    if (seen.has(ccp4Path)) continue;
    seen.add(ccp4Path);
    const pythonPath = path.join(ccp4Path, "bin", pythonBinName);
    if (!fs.existsSync(pythonPath)) continue;
    const verdict = assessPython(pythonPath);
    if (verdict.supported) return ccp4Path;
    console.log(`Skipping ${ccp4Path} as the default CCP4: ${verdict.reason}`);
  }

  // Nothing found: say so, rather than naming a location we have just checked
  // and know is not there. Every path above is returned only after confirming it
  // holds bin/ccp4-python, so a returned value always means a real install.
  //
  // The previous fallback returned /Applications/ccp4-9 (macOS) or the platform
  // equivalent regardless, which put a confident, wrong, and on macOS actively
  // unsuitable path in front of every first-time user — it reads as "this is
  // where CCP4 should go" rather than "we could not find CCP4". Empty makes the
  // setup page ask, which is the honest thing on a machine we know nothing about.
  return "";
};

// Config lives under the CCP4i2 home rather than Electron's userData. userData
// also holds Chromium's cache, GPU state and cookies, which have no business in
// a directory we tell users is their CCP4i2 configuration — and unpackaged runs
// put it in ~/.config/Electron, which is not namespaced to this app at all.
export const store = new Store<StoreSchema>({
  cwd: ccp4i2Home(),
  name: "electron",
  defaults: {
    CCP4Dir: getDefaultCCP4Dir(),
    projectRoot: getProjectRoot(), // Computed, not user-configurable
    devMode: false,
    zoomLevel: 0,
    CCP4I2_PROJECTS_DIR: getProjectsDir(),
    theme: "dark",
    autoLaunch: true,
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
    Menu.setApplicationMenu(null);
    setupZoomLevel(store);
    // Set both API_BASE_URL (runtime, for server-side proxy routes) and
    // NEXT_PUBLIC_API_BASE_URL (for any client-side code that needs it)
    process.env.API_BASE_URL = `http://localhost:${djangoServerPort}`;
    process.env.NEXT_PUBLIC_API_BASE_URL = `http://localhost:${djangoServerPort}`;
    nextServer = await startNextServer(isDev, nextServerPort, djangoServerPort);
  })
  .then(async () => {
    // Use /ccp4i2 base path for multi-app integration
    mainWindow = await createWindow(
      `http://localhost:${nextServerPort}/ccp4i2/config`,
      store
    );
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
    // Use /ccp4i2 base path for multi-app integration
    mainWindow = await createWindow(
      `http://localhost:${nextServerPort}/ccp4i2`,
      store
    );
  }
});
