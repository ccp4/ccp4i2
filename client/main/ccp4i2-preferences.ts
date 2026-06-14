/**
 * Shared CCP4i2 preferences — the Electron-side reader/writer for
 * `~/.ccp4i2/preferences.json`.
 *
 * This is the TypeScript mirror of `server/ccp4i2/config/preferences.py`: same
 * location, same schema, same atomic-write discipline. It exists so the desktop
 * GUI writes the *shared* bootstrap keys (CCP4 install, projects dir, …) to the
 * one file that the Django control plane and the `i2`/`i2run` CLI also read —
 * instead of burying them in `electron-store`'s `userData` where the CLI never
 * looks.
 *
 * Only *shared* preferences belong here. GUI-only chrome (window geometry, zoom,
 * theme) stays in `electron-store`.
 *
 * Resolution precedence on the Python side is `env var > preferences.json >
 * default`; this module is the file layer's writer (and a reader for display).
 */
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

const PREFERENCES_VERSION = 1;

export interface CCP4i2Preferences {
  version?: number;
  ccp4Dir?: string;
  projectsDir?: string;
  database?: string;
  ccp4i2Root?: string;
  userPreferences?: Record<string, unknown>;
}

/** The CCP4i2 user home (`CCP4I2_HOME` env var, else `~/.ccp4i2`). */
export function ccp4i2Home(): string {
  const override = process.env.CCP4I2_HOME;
  return override ? path.resolve(override) : path.join(os.homedir(), ".ccp4i2");
}

/** Full path to `preferences.json` inside the CCP4i2 user home. */
export function preferencesPath(): string {
  return path.join(ccp4i2Home(), "preferences.json");
}

/** Load preferences; return `{}` if the file is absent or unreadable. */
export function loadPreferences(): CCP4i2Preferences {
  try {
    const text = fs.readFileSync(preferencesPath(), "utf-8");
    const data = JSON.parse(text);
    return data && typeof data === "object" ? (data as CCP4i2Preferences) : {};
  } catch {
    return {};
  }
}

/**
 * Atomically write `preferences.json` (temp file + rename), stamping the schema
 * version. Last writer wins — adequate for this small, rarely-contended file.
 */
export function savePreferences(prefs: CCP4i2Preferences): void {
  const home = ccp4i2Home();
  fs.mkdirSync(home, { recursive: true });
  const payload = { version: PREFERENCES_VERSION, ...prefs };
  const tmp = path.join(home, `.preferences-${process.pid}-${Date.now()}.tmp`);
  fs.writeFileSync(tmp, JSON.stringify(payload, null, 2) + "\n", "utf-8");
  fs.renameSync(tmp, preferencesPath());
}

/**
 * Merge a patch of top-level keys into the file, preserving everything else
 * (including the `userPreferences` bag). Returns the merged preferences.
 */
export function updatePreferences(patch: CCP4i2Preferences): CCP4i2Preferences {
  const merged = { ...loadPreferences(), ...patch };
  savePreferences(merged);
  return merged;
}

/**
 * Build a `settings.py`-parseable SQLite URL for an absolute DB path.
 *
 * Produces `sqlite:///<path>` with forward slashes. The Python side
 * (`config/settings.py`) strips the leading slash before a Windows drive letter
 * (`/C:/...` -> `C:/...`), so this works on both POSIX and Windows:
 *   POSIX   /data/proj/db.sqlite3   -> sqlite:///data/proj/db.sqlite3
 *   Windows C:\proj\db.sqlite3      -> sqlite:///C:/proj/db.sqlite3
 */
export function sqliteUrl(dbPath: string): string {
  const forward = dbPath.replace(/\\/g, "/").replace(/^\/+/, "");
  return "sqlite:///" + forward;
}
