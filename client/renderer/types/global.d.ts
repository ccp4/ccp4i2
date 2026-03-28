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
import { ipcRenderer } from "electron";

declare global {
  interface Window {
    electronAPI: {
      onMessage: (
        channel: string,
        callback: (event: Electron.IpcRendererEvent, ...args: any[]) => void
      ) => typeof ipcRenderer.on;
      removeMessageListener: (
        channel: string,
        callback: (event: Electron.IpcRendererEvent, ...args: any[]) => void
      ) => typeof ipcRenderer.off;
      sendMessage: (channel: string, ...args: any[]) => typeof ipcRenderer.send;
      sendSync: (channel: string, ...args: any[]) => any;
    };
  }
}
export {};
