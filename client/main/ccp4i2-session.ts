import { app, BrowserWindow, dialog } from "electron";
import path from "path";

/**
 * Installs a handler for the "will-download" event on the provided Electron session.
 * This handler intercepts download events, prompts the user to select a save location,
 * and sets the file's save path accordingly. If the user cancels the save dialog,
 * the download is canceled.
 *
 * @param session - The Electron session to attach the "will-download" event handler to.
 */
export function installWillDownloadHandler(session: Electron.Session) {
  // Intercept downloads
  session.on("will-download", (event, item, webContents) => {
    const browserWindow = BrowserWindow.fromWebContents(webContents);
    if (!browserWindow) {
      console.error("Browser window not found");
      return;
    }

    // Prompt user for a save path synchronously
    const filePath = dialog.showSaveDialogSync(browserWindow, {
      title: "Save File",
      defaultPath: path.join(app.getPath("downloads"), item.getFilename()),
      buttonLabel: "Save",
    });

    if (filePath) {
      console.log("Download will be saved to:", filePath);
      item.setSavePath(filePath);
    } else {
      console.log("Download canceled by user.");
      item.cancel();
    }

    // Optional: Monitor download result
    item.once("done", (event, state) => {
      if (state === "completed") {
        console.log("Download completed successfully:", item.getSavePath());
      } else {
        console.error("Download failed or was cancelled:", state);
      }
    });
  });
}
