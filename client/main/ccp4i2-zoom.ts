import { BrowserWindow, Input, Menu, MenuItemConstructorOptions } from "electron";
import ElectronStore from "electron-store";
import { StoreSchema } from "../types/store";

/**
 * Zoom, and the native application menu.
 *
 * Zoom was bound with `globalShortcut`, which registers OS-wide hotkeys: on
 * Linux under Wayland (Ubuntu 24.04's default session) they are refused, so
 * Ctrl+/Ctrl- did nothing; "Plus" meant Shift+=, so Ctrl+= never zoomed; and
 * while the app ran it took Ctrl/Cmd +, - and 0 from every other application.
 * The keys are now read per window (`zoomStepForKey`, from the window's
 * before-input-event handler), and the in-app View menu, the macOS menu and
 * the keys all go through `zoomAll`, so the level is kept the same way
 * whichever was used.
 *
 * The app clears the native menu so it does not duplicate the in-app menu
 * bar. On Windows and Linux that removes the bar. macOS always has a menu
 * bar, and clearing it left Electron's default View menu standing with every
 * item disabled; since macOS hands key equivalents to the menu bar before the
 * window, the dead "Toggle Developer Tools" swallowed Cmd+Opt+I (and Cmd+R,
 * Cmd+Q were greyed). macOS gets a real application menu instead.
 */

export type ZoomStep = "in" | "out" | "reset";

/** Change the zoom of every window and remember it. */
export function zoomAll(store: ElectronStore<StoreSchema>, step: ZoomStep) {
  BrowserWindow.getAllWindows().forEach((win) => {
    const current = win.webContents.getZoomLevel();
    const next = step === "in" ? current + 1 : step === "out" ? current - 1 : 0;
    win.webContents.setZoomLevel(next);
    store.set("zoomLevel", next);
  });
}

/** The zoom a key press asks for, or null. Ctrl (Cmd on macOS) with +, = or
 *  numpad + zooms in, with - or numpad - out, with 0 back to actual size;
 *  Shift does not matter, so Ctrl+= and Ctrl+Shift+= both zoom in. */
export function zoomStepForKey(
  input: Pick<Input, "type" | "key" | "control" | "meta" | "alt">,
  isMac: boolean
): ZoomStep | null {
  if (input.type !== "keyDown") return null;
  const primary = isMac ? input.meta : input.control;
  if (!primary || input.alt) return null;
  switch (input.key) {
    case "+":
    case "=":
      return "in";
    case "-":
    case "_":
      return "out";
    case "0":
      return "reset";
    default:
      return null;
  }
}

/** The native application menu: a standard one on macOS, none elsewhere. */
export function installApplicationMenu(store: ElectronStore<StoreSchema>) {
  if (process.platform !== "darwin") {
    Menu.setApplicationMenu(null);
    return;
  }
  const template: MenuItemConstructorOptions[] = [
    { role: "appMenu" },
    { role: "editMenu" },
    {
      label: "View",
      submenu: [
        { role: "reload" },
        { role: "forceReload" },
        { role: "toggleDevTools" },
        { type: "separator" },
        { label: "Actual Size", accelerator: "CmdOrCtrl+0", click: () => zoomAll(store, "reset") },
        { label: "Zoom In", accelerator: "CmdOrCtrl+Plus", click: () => zoomAll(store, "in") },
        { label: "Zoom Out", accelerator: "CmdOrCtrl+-", click: () => zoomAll(store, "out") },
        { type: "separator" },
        { role: "togglefullscreen" },
      ],
    },
    { role: "windowMenu" },
  ];
  Menu.setApplicationMenu(Menu.buildFromTemplate(template));
}
