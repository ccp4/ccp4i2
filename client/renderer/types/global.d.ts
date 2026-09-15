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
      invoke?: (channel: string, data?: any) => Promise<any>;
      /** Absolute local path of a File the renderer holds, or "" if it has none.
       *  Desktop only (webUtils.getPathForFile); undefined in the web build. */
      getPathForFile?: (file: File) => string;
    };
  }
}
export {};
