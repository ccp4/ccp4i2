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
 * Two entry points:
 *   - a silent check at launch (initAutoUpdater), and
 *   - a manual "Check for Updates…" from the Help menu, which sends the
 *     "check-for-updates" IPC and gets a visible result dialog either way.
 *
 * Platform coverage (deliberate, documented):
 *   - Windows (NSIS) and Linux AppImage: full auto-update.
 *   - macOS: only when the build is SIGNED (Squirrel.Mac validates the
 *     signature; an unsigned app downloads but cannot apply). Gated on
 *     ENABLE_MAC_SIGNING in release.yml. A signed app cannot auto-update FROM an
 *     unsigned one, so the first signed build must be installed by hand once.
 *   - Linux .deb: electron-updater has no .deb mechanism — we skip the check
 *     entirely and leave those users to CCP4 UM / apt / a manual download.
 *
 * electron-updater is CommonJS; its ESM named export is unreliable, so we
 * default-import and destructure (the documented interop pattern).
 *
 * All update activity is logged to a file via electron-log — mac
 * ~/Library/Logs/ccp4i2x/main.log, Windows %APPDATA%/ccp4i2x/logs/main.log,
 * Linux ~/.config/ccp4i2x/logs/main.log — so a failed check can be diagnosed
 * by opening a file rather than relaunching from a terminal (which is what the
 * bring-up of this feature actually needed). electron-updater's own internal
 * logs go there too (autoUpdater.logger = log).
 */
import { app, dialog, ipcMain } from "electron";
import type { BrowserWindow } from "electron";
import electronUpdater from "electron-updater";
import log from "electron-log/main";

const { autoUpdater } = electronUpdater;

let getWin: () => BrowserWindow | null = () => null;
let handlersWired = false;

/** Why (if at all) auto-update cannot run for this install. */
function updaterUnavailableReason(): string | null {
  if (!app.isPackaged) return "Updates apply to installed builds only (this is a dev run).";
  // On Linux, electron-updater only supports the AppImage, which sets $APPIMAGE.
  if (process.platform === "linux" && !process.env.APPIMAGE) {
    return "This is a .deb install; update it with your package manager.";
  }
  return null;
}

/** Configure the updater and register its event handlers exactly once. */
function wireOnce(): void {
  if (handlersWired) return;
  handlersWired = true;

  // Route electron-updater's own logs to the file (and console) so a failed
  // check — "No published versions on GitHub", a 404 on the channel file, a
  // download error — lands in main.log instead of vanishing with stdout.
  autoUpdater.logger = log;
  log.transports.file.level = "info";

  // Our releases are pre-releases (…-alpha.56); without this electron-updater
  // ignores them and never offers an alpha→alpha update.
  autoUpdater.allowPrerelease = true;
  // Fetch in the background; prompt only once it is ready to apply.
  autoUpdater.autoDownload = true;
  // If the user defers the restart, apply the update on the next quit anyway.
  autoUpdater.autoInstallOnAppQuit = true;

  autoUpdater.on("update-available", (info) => {
    log.info(`[updater] update available: ${info?.version}`);
  });
  autoUpdater.on("update-not-available", () => {
    log.info("[updater] no update available.");
  });

  autoUpdater.on("update-downloaded", async (info) => {
    log.info(`[updater] update downloaded: ${info?.version}`);
    const win = getWin();
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
    // A silent (launch) check must never interrupt work — an unsigned mac
    // build, an offline machine, or an unsupported package all land here. The
    // manual path reports its own errors via the promise, not this handler.
    log.warn(`[updater] error (non-fatal): ${err?.message ?? err}`);
  });
}

/**
 * Run a check. `interactive` = triggered from the menu, so always end with a
 * visible dialog (found / up to date / unavailable / failed). A non-interactive
 * (launch) check is silent unless an update is found and downloaded.
 */
async function runCheck(interactive: boolean): Promise<void> {
  const win = () => getWin();
  const info = (message: string, detail?: string) => {
    if (!interactive) return;
    const opts = { type: "info" as const, title: "Software Update", message, detail };
    const w = win();
    if (w) dialog.showMessageBox(w, opts);
    else dialog.showMessageBox(opts);
  };

  const reason = updaterUnavailableReason();
  if (reason) {
    info("Automatic updates aren't available here.", reason);
    return;
  }

  wireOnce();

  try {
    const result = await autoUpdater.checkForUpdates();
    const available = result?.isUpdateAvailable;
    if (available) {
      info(
        "An update is available and downloading.",
        "You'll be asked to restart when it's ready.",
      );
    } else {
      info(`You're up to date.`, `CCP4i2x ${app.getVersion()} is the latest version.`);
    }
  } catch (err) {
    log.warn(`[updater] checkForUpdates failed: ${(err as Error)?.message ?? err}`);
    info(
      "Couldn't check for updates.",
      "Please check your connection and try again later.",
    );
  }
}

/**
 * Wire the updater and do a silent check at launch. Also registers the
 * "check-for-updates" IPC so the Help-menu item works on every platform (it
 * reports "not available here" rather than doing nothing on unsupported
 * installs). Safe to call unconditionally; failures never block the app.
 */
export function initAutoUpdater(
  getMainWindow: () => BrowserWindow | null
): void {
  getWin = getMainWindow;

  // Register the manual trigger unconditionally so the menu item always gives
  // feedback, even in dev or on a .deb where the auto-check is skipped.
  ipcMain.removeAllListeners("check-for-updates");
  ipcMain.on("check-for-updates", () => {
    void runCheck(true);
  });

  // Silent check at launch only where updates can actually apply.
  if (updaterUnavailableReason()) return;
  void runCheck(false);
}
