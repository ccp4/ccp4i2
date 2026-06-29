/**
 * Single source of truth for which ccp4i2 backend the *packaged* app expects.
 *
 * The packaged Electron app installs the django backend from PyPI and treats
 * anything older than this floor as "not ready" (offering an upgrade). The
 * value is a floor, not an exact pin: install does
 *   pip install --upgrade "ccp4i2>=<floor>"
 * so users automatically pick up newer compatible releases without a rebuild,
 * while still guaranteeing a minimum tested backend.
 *
 * This is deliberately decoupled from the Electron app's own package.json
 * version (currently 0.0.1) — the app and the python package version on
 * independent cadences.
 *
 * Dev mode ignores this entirely: an unpacked dev electron uses the local
 * `-e ./server` checkout regardless of version (mode is keyed strictly on
 * app.isPackaged).
 *
 * To bump the floor for a release, either edit CCP4I2_SERVER_VERSION_FLOOR
 * below, or set CCP4I2_SERVER_VERSION_FLOOR in the build environment so the
 * build can stamp it without a code edit.
 */
// 3.0.1 is the floor, not 3.0.0: the published 3.0.0 wheel predates the bundled
// requirements-runtime.txt, so its two-step install can't find the lock and the
// dep closure (modern django/asgiref) is never applied. 3.0.1 is the first wheel
// that ships the lock. Requiring >=3.0.1 ensures a packaged app installs/upgrades
// to a version that can actually complete the install.
export const CCP4I2_SERVER_VERSION_FLOOR =
  process.env.CCP4I2_SERVER_VERSION_FLOOR || "3.0.1";

/**
 * Compare two dotted numeric version strings (e.g. "3.0.10" vs "3.0.2").
 * Returns a negative number if a < b, 0 if equal, positive if a > b.
 * Non-numeric / missing components are treated as 0. Pre-release suffixes
 * (e.g. "3.0.0rc1") are ignored beyond the numeric prefix of each component,
 * which is sufficient for the floor check.
 */
export function compareVersions(a: string, b: string): number {
  const parse = (v: string) =>
    v
      .split(".")
      .map((part) => parseInt(part, 10))
      .map((n) => (Number.isNaN(n) ? 0 : n));
  const pa = parse(a);
  const pb = parse(b);
  const len = Math.max(pa.length, pb.length);
  for (let i = 0; i < len; i++) {
    const da = pa[i] ?? 0;
    const db = pb[i] ?? 0;
    if (da !== db) return da - db;
  }
  return 0;
}

/**
 * Is the installed ccp4i2 version acceptable for a packaged app?
 * True when it meets or exceeds the floor. An unparseable / missing version
 * is treated as acceptable (the ASGI app imported, so the server will boot —
 * we don't want a parse hiccup to block a working install).
 */
export function meetsServerVersionFloor(installed: string | undefined): boolean {
  if (!installed) return true;
  return compareVersions(installed, CCP4I2_SERVER_VERSION_FLOOR) >= 0;
}
