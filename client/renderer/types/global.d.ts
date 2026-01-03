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
    };
  }
}
export {};
