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

export const HOME_DIR_NAME = ".ccp4i2x";
export const LEGACY_HOME_DIR_NAME = ".ccp4i2-django";

/**
 * The CCP4i2 user home: everything per-user lives under one root.
 *
 * 1. `CCP4I2_HOME` wins outright.
 * 2. An existing `~/.ccp4i2x`.
 * 3. An existing `~/.ccp4i2-django` — adopted in place, so installs from
 *    3.1.0a26 and earlier keep working exactly where they are.
 * 4. `~/.ccp4i2x` for a fresh install.
 *
 * `~/.ccp4i2x` and not `~/.ccp4i2`: the latter differs from the legacy Qt-era
 * home `~/.CCP4I2` only by case, so on a case-insensitive filesystem (the macOS
 * default) they collide and the new app would plant its db and preferences in
 * the tree it is migrating from.
 *
 * MUST match the Python side (server/ccp4i2/config/preferences.py). When these
 * two disagree the app and the server use different databases and different
 * project stores — which is exactly what happened before this was shared.
 */
export function ccp4i2Home(): string {
  const override = process.env.CCP4I2_HOME;
  if (override) return path.resolve(override);

  const home = os.homedir();
  const current = path.join(home, HOME_DIR_NAME);
  if (isDirectory(current)) return current;

  const legacy = path.join(home, LEGACY_HOME_DIR_NAME);
  if (isDirectory(legacy)) return legacy;

  return current;
}

function isDirectory(candidate: string): boolean {
  try {
    return fs.statSync(candidate).isDirectory();
  } catch {
    return false;
  }
}

/**
 * Where new projects go unless the user says otherwise: `<home>/projects`.
 * An adopted pre-a27 home keeps its existing `CCP4X_PROJECTS`, because renaming
 * it would strand every absolute path recorded inside those projects.
 */
export function defaultProjectsDir(): string {
  const home = ccp4i2Home();
  const legacy = path.join(home, "CCP4X_PROJECTS");
  if (isDirectory(legacy) && !isDirectory(path.join(home, "projects"))) {
    return legacy;
  }
  return path.join(home, "projects");
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
