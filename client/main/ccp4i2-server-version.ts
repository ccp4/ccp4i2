/**
 * Single source of truth for which ccp4i2 backend the *packaged* app expects.
 *
 * During the pre-release (alpha) phase the packaged app pins the backend
 * EXACTLY: it installs and requires one specific ccp4i2 version, e.g.
 *   pip install --no-deps --upgrade "ccp4i2==<version>"
 * so an alpha desktop build and its backend are strictly bound. An alpha app
 * refuses to run against any other version — including the un-suffixed 3.0.x
 * releases that escaped before this discipline. Exact pinning trades
 * "auto-pick-up-newer" for "cannot silently mix versions", which is what we want
 * while the software is not yet world-ready. Relax to a floor / compatible-range
 * (a >= check via compareVersions) once a stable line is declared.
 *
 * Deliberately decoupled from the Electron app's own package.json version
 * (currently 0.0.1) — the app and the python package version on independent
 * cadences.
 *
 * Dev mode ignores this entirely: an unpacked dev electron uses the local
 * `-e ./server` checkout regardless of version (keyed strictly on
 * app.isPackaged).
 *
 * The value is stamped at build time from the release tag via the
 * CCP4I2_SERVER_VERSION_FLOOR build-env var. The env var name is kept for
 * build-workflow compatibility, but it now carries the EXACT required version,
 * not a floor.
 */
export const CCP4I2_REQUIRED_SERVER_VERSION =
  process.env.CCP4I2_SERVER_VERSION_FLOOR || "3.1.0a20";

/**
 * Compare two dotted numeric version strings (e.g. "3.0.10" vs "3.0.2").
 * Returns a negative number if a < b, 0 if equal, positive if a > b.
 * Non-numeric / missing components are treated as 0.
 *
 * NOTE: this ignores PEP 440 pre-release suffixes (it parses only the numeric
 * dotted prefix), so it CANNOT distinguish 3.1.0a1 from 3.1.0a2. It is kept for
 * a future floor/range check; the exact-pin gate below uses string equality
 * instead. Do not use this to compare pre-release versions.
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
 * Normalise a version for exact comparison: lower-case and drop the separators
 * pip treats as equivalent, so "3.1.0-a1" and "3.1.0a1" compare equal. (Our own
 * __version__ always emits the canonical "3.1.0a1" form; this just makes the
 * check tolerant of hand-typed / stamped variants.)
 */
function normalizeVersion(v: string): string {
  return v.trim().toLowerCase().replace(/[-_]/g, "");
}

/**
 * Does the installed ccp4i2 EXACTLY match the version this app is pinned to?
 *
 * Strict during the alpha: the app accepts ONLY its exact partner version.
 * A missing / unreadable installed version counts as a MISMATCH, not a pass —
 * the ASGI app importing proves the backend *boots*, but the pin is about
 * version *identity*, and "I can't tell what's installed" is not "the right
 * thing is installed". Treating unknown as acceptable was the hole that let a
 * wrong/undetectable backend through; here it fails closed so the UI offers to
 * install the exact partner. (Dev/unpacked mode never calls this — see caller.)
 */
export function meetsServerVersionRequirement(
  installed: string | undefined,
): boolean {
  if (!installed) return false;
  return (
    normalizeVersion(installed) ===
    normalizeVersion(CCP4I2_REQUIRED_SERVER_VERSION)
  );
}
