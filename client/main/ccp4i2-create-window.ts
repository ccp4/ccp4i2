/*
 * Copyright (C) 2025-2026 Newcastle University
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
