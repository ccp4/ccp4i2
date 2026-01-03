// electron-global.d.ts
import { ipcRenderer } from "electron";

declare global {
  // Adds the Electron API to globalThis
  var electron: {
    onMessage: (
      channel: string,
      callback: (event: Electron.IpcRendererEvent, ...args: any[]) => void
    ) => typeof ipcRenderer.on;
    removeMessageListener: (
      channel: string,
      callback: (event: Electron.IpcRendererEvent, ...args: any[]) => void
    ) => typeof ipcRenderer.on;
    sendMessage: (channel: string, ...args: any[]) => typeof ipcRenderer.send;
  };
}
