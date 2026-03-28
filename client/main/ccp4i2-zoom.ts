/*
 * Copyright (C) 2025 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
