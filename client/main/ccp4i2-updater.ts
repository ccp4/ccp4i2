/**
 * In-app auto-update for the packaged desktop app (alpha bootstrap).
 *
 * The app auto-updates itself from the ccp4/ccp4i2 GitHub Releases via
 * electron-updater; the *backend* then follows automatically because the app
 * exact-pins `ccp4i2==<its version>` and the readiness probe pip-installs the
 * matching wheel on next launch (see ccp4i2-server-version.ts and the Install
 * flow in ccp4i2-ipc.ts). So one app update pulls its partner backend with it —
 * the lockstep we already have is the backend's update trigger. This is the
 * alpha form of the plan in docs/UPDATE_MECHANISM_PLAN.md; GA decouples them.
 *
 * Platform coverage (deliberate, documented):
 *   - Windows (NSIS) and Linux AppImage: full auto-update.
 *   - macOS: only when the build is SIGNED (Squirrel.Mac validates the
 *     signature; an unsigned app downloads but cannot apply). Gated on
 *     ENABLE_MAC_SIGNING in release.yml.
 *   - Linux .deb: electron-updater has no .deb mechanism — we skip the check
 *     entirely and leave those users to CCP4 UM / apt / a manual download.
 *
 * electron-updater is CommonJS; its ESM named export is unreliable, so we
 * default-import and destructure (the documented interop pattern).
 */
import { app, dialog } from "electron";
import type { BrowserWindow } from "electron";
import electronUpdater from "electron-updater";

const { autoUpdater } = electronUpdater;

/**
 * Wire and kick off the auto-updater. Safe to call unconditionally: it no-ops
 * in dev and on package types electron-updater cannot update, and every failure
 * is swallowed (a broken update check must never block using the app).
 */
export function initAutoUpdater(
  getMainWindow: () => BrowserWindow | null
): void {
  // Dev (unpacked) has no installed artifact to replace.
  if (!app.isPackaged) return;

  // On Linux, electron-updater only supports the AppImage, which sets $APPIMAGE.
  // Its absence means a .deb (or an unpacked run) — nothing to update, so skip
  // the check rather than emit a confusing "no published versions" error.
  if (process.platform === "linux" && !process.env.APPIMAGE) {
    console.log("[updater] Linux non-AppImage build — auto-update skipped.");
    return;
  }

  // Our releases are GitHub pre-releases (…-a52); without this electron-updater
  // ignores them and never offers an alpha→alpha update.
  autoUpdater.allowPrerelease = true;
  // Fetch in the background; prompt only once it is ready to apply.
  autoUpdater.autoDownload = true;
  // If the user defers the restart, apply the update on the next quit anyway.
  autoUpdater.autoInstallOnAppQuit = true;

  autoUpdater.on("update-available", (info) => {
    console.log(`[updater] update available: ${info?.version}`);
  });

  autoUpdater.on("update-not-available", () => {
    console.log("[updater] no update available.");
  });

  autoUpdater.on("update-downloaded", async (info) => {
    console.log(`[updater] update downloaded: ${info?.version}`);
    const win = getMainWindow();
    const opts = {
      type: "info" as const,
      buttons: ["Restart now", "Later"],
      defaultId: 0,
      cancelId: 1,
      title: "Update ready",
      message: `CCP4i2x ${info?.version} has been downloaded.`,
      detail:
        "Restart to apply it. The matching backend will be installed on the " +
        "next launch. You can also keep working — it will apply when you quit.",
    };
    const { response } = win
      ? await dialog.showMessageBox(win, opts)
      : await dialog.showMessageBox(opts);
    if (response === 0) {
      // quitAndInstall on macOS needs a signed app; if it cannot apply, the
      // error handler below swallows it and autoInstallOnAppQuit still tries.
      autoUpdater.quitAndInstall();
    }
  });

  autoUpdater.on("error", (err) => {
    // Never surface to the user: an unsigned mac build, an offline machine, or
    // an unsupported package all land here and none should interrupt work.
    console.log(`[updater] check failed (non-fatal): ${err?.message ?? err}`);
  });

  autoUpdater.checkForUpdates().catch((err) => {
    console.log(`[updater] checkForUpdates threw (non-fatal): ${err?.message ?? err}`);
  });
}
