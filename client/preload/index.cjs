import { contextBridge, ipcRenderer, webUtils } from "electron";

contextBridge.exposeInMainWorld("electronAPI", {
  sendMessage: (channel, data) => ipcRenderer.send(channel, data),
  sendSync: (channel, data) => ipcRenderer.sendSync(channel, data),
  onMessage: (channel, callback) => ipcRenderer.on(channel, callback),
  removeMessageListener: (channel, callback) => {
    ipcRenderer.removeListener(channel, callback);
  },
  // Promise-returning IPC (ipcMain.handle). Used e.g. by the Program locations
  // preferences panel to open a native file/folder picker and get the path back.
  invoke: (channel, data) => ipcRenderer.invoke(channel, data),
  // The absolute local path of a File the renderer holds (from a file input or
  // drag-drop). Electron 32 removed File.path; webUtils.getPathForFile is the
  // sanctioned replacement, and must be called here in the preload. Returns "" for
  // a File with no backing path (a synthesized Blob). Lets the desktop import a
  // large file by path (server-side copy) instead of uploading its bytes -- see
  // #512. Absent in the web build (no preload).
  getPathForFile: (file) => {
    try {
      return webUtils.getPathForFile(file) || "";
    } catch {
      return "";
    }
  },
});

// LocalSession surface — exposes the per-launch token and the
// associated user email set by Electron main on process.env.
// Read-only by virtue of being a frozen object literal handed to
// contextBridge; the renderer can read but not change.
// In the web build there is no preload, so window.ccp4i2LocalSession
// is undefined there and the renderer falls back to MSAL.
const localSessionToken = process.env.CCP4I2_LOCAL_SESSION_TOKEN;
const localSessionUserEmail = process.env.CCP4I2_LOCAL_USER_EMAIL;
if (typeof localSessionToken === "string" && localSessionToken.length > 0) {
  contextBridge.exposeInMainWorld("ccp4i2LocalSession", {
    token: localSessionToken,
    userEmail: localSessionUserEmail,
  });
}

console.log("[Preload] Script loaded");
