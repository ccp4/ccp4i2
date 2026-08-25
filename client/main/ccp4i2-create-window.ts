import { BrowserWindow } from "electron";
import path from "path";
import { fileURLToPath } from "node:url";
import ElectronStore from "electron-store";
import { StoreSchema } from "../types/store";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

export const createWindow = async (
  url: string,
  store: ElectronStore<StoreSchema>
) => {
  const newWindow = new BrowserWindow({
    width: 1000,
    height: 700,
    webPreferences: {
      preload: path.join(__dirname, "preload.js"),
      contextIsolation: true,
      nodeIntegration: false,
    },
  });
  // Intercept window.open and create a new Electron window instead
  newWindow.webContents.setWindowOpenHandler(({ url }) => {
    createWindow(url, store);
    return { action: "deny" }; // Prevent default behavior
  });
  // The app clears the native application menu (ccp4i2-master.ts) so it doesn't
  // duplicate the in-app menu bar — which also takes the standard View-menu
  // accelerators with it, leaving Reload and Toggle Developer Tools greyed out
  // and their shortcuts dead. Re-register them per window, so they work in
  // report pop-ups (opened through the handler above) as well as the main one.
  newWindow.webContents.on("before-input-event", (event, input) => {
    if (input.type !== "keyDown") return;
    const isMac = process.platform === "darwin";
    const key = input.key.toLowerCase();
    const primary = isMac ? input.meta : input.control;

    if (key === "f12" || (primary && input.shift && key === "i") ||
        (isMac && input.meta && input.alt && key === "i")) {
      newWindow.webContents.toggleDevTools();
      event.preventDefault();
      return;
    }
    if (primary && key === "r") {
      if (input.shift) {
        newWindow.webContents.reloadIgnoringCache();
      } else {
        newWindow.webContents.reload();
      }
      event.preventDefault();
    }
  });

  newWindow.webContents.on("did-finish-load", () => {
    // Restore the persisted zoom level, but guard against corrupted/out-of-range
    // values. A stray accumulated value (e.g. -10.5) renders the whole UI
    // unusably tiny; and because this runs on *every* load — including dev
    // hot-reloads — such a value also defeats menu/keyboard zoom, which gets
    // reset on the next reload. Reset anything outside a sane range back to 0.
    let zoom = store.get("zoomLevel");
    if (typeof zoom !== "number" || !Number.isFinite(zoom) || zoom < -3 || zoom > 3) {
      zoom = 0;
      store.set("zoomLevel", zoom);
    }
    newWindow.webContents.setZoomLevel(zoom);
  });
  // Relay find-in-page results back to renderer
  newWindow.webContents.on("found-in-page", (_event, result) => {
    newWindow.webContents.send("message-from-main", {
      message: "found-in-page",
      activeMatchOrdinal: result.activeMatchOrdinal,
      matches: result.matches,
    });
  });
  setTimeout(() => newWindow?.loadURL(url), 1500);
  return newWindow;
};
