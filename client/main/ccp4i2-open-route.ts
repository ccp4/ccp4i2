/**
 * "--open-route <path>": open one of the app's own routes in a new window.
 *
 * The one thing outside the renderer that needs to open a window is a
 * process on the same machine: the plugin of a recorded Moorhen session
 * started from `i2run moorhen ...`, which has a job but no window. It
 * spawns the running app's own executable with this argument; the
 * single-instance lock forwards the argv to the running instance, which
 * opens the route. Nothing Moorhen-specific here: any route the app
 * serves can be opened this way. See docs/moorhen-task-design.md.
 *
 * Pure functions only (no electron import), so they can be unit tested.
 */

export const OPEN_ROUTE_FLAG = "--open-route";

/** A route the app serves: an absolute path under /ccp4i2/, no scheme, no
 *  host, no protocol-relative trick, no whitespace. Anything else is
 *  refused rather than loaded. */
export function isOpenableRoute(route: string): boolean {
  if (typeof route !== "string") return false;
  if (!/^\/ccp4i2\/[A-Za-z0-9._~\-\/?=&%+]*$/.test(route)) return false;
  if (route.startsWith("//")) return false;
  return true;
}

/** The route named by `--open-route <path>` or `--open-route=<path>` in an
 *  argv, or null if there is none or it is not openable. */
export function parseOpenRoute(argv: readonly string[]): string | null {
  for (let i = 0; i < argv.length; i++) {
    const arg = argv[i];
    let candidate: string | undefined;
    if (arg === OPEN_ROUTE_FLAG) candidate = argv[i + 1];
    else if (arg.startsWith(`${OPEN_ROUTE_FLAG}=`)) candidate = arg.slice(OPEN_ROUTE_FLAG.length + 1);
    if (candidate === undefined) continue;
    return isOpenableRoute(candidate) ? candidate : null;
  }
  return null;
}

/** The command that starts (or, under the single-instance lock, signals)
 *  this app: the executable plus, in development, the app directory that
 *  `electron` needs as its first argument. Exported to the Django child as
 *  CCP4I2_DESKTOP_LAUNCH (a JSON array) so a job process can append
 *  `--open-route <path>` and spawn it. */
export function desktopLaunchCommand(execPath: string, isPackaged: boolean, appPath: string): string[] {
  return isPackaged ? [execPath] : [execPath, appPath];
}
