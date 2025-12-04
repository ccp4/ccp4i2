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
  newWindow.webContents.on("did-finish-load", () => {
    newWindow.webContents.setZoomLevel(store.get("zoomLevel"));
  });
  setTimeout(() => newWindow?.loadURL(url), 1500);
  return newWindow;
};
