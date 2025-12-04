import { globalShortcut, BrowserWindow } from "electron";
import ElectronStore from "electron-store";
import { StoreSchema } from "../types/store";

export function setupZoomLevel(store: ElectronStore<StoreSchema>) {
  globalShortcut.register("CommandOrControl+Plus", () => {
    BrowserWindow.getAllWindows().forEach((win) => {
      const current = win.webContents.getZoomLevel();
      win.webContents.setZoomLevel(current + 1);
      store.set("zoomLevel", current + 1);
    });
  });

  globalShortcut.register("CommandOrControl+-", () => {
    BrowserWindow.getAllWindows().forEach((win) => {
      const current = win.webContents.getZoomLevel();
      win.webContents.setZoomLevel(current - 1);
      store.set("zoomLevel", current - 1);
    });
  });

  globalShortcut.register("CommandOrControl+0", () => {
    BrowserWindow.getAllWindows().forEach((win) => {
      const current = win.webContents.getZoomLevel();
      win.webContents.setZoomLevel(0);
      store.set("zoomLevel", 0);
    });
  });
}
