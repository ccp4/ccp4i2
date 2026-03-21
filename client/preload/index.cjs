import { contextBridge, ipcRenderer } from "electron";

contextBridge.exposeInMainWorld("electronAPI", {
  sendMessage: (channel, data) => ipcRenderer.send(channel, data),
  sendSync: (channel, data) => ipcRenderer.sendSync(channel, data),
  onMessage: (channel, callback) => ipcRenderer.on(channel, callback),
  removeMessageListener: (channel, callback) => {
    ipcRenderer.removeListener(channel, callback);
  },
});
console.log("[Preload] Script loaded");
