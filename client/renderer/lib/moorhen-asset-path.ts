/**
 * Where the Moorhen runtime assets live, and why the URL carries a version.
 *
 * The API route that serves these files answers with
 * `Cache-Control: public, max-age=31536000, immutable`, and it has to: each
 * pthread worker is started with `new Worker("moorhen.js")`, so a single window
 * loads that one script 33 times. Revalidating would put 33 conditional
 * requests through a route that reads the file from disk each time.
 *
 * But `immutable` is a promise that the bytes at a URL never change, and until
 * this module the URL was the fixed `/api/moorhen/MoorhenAssets`. Every Moorhen
 * bump changes moorhen.wasm, moorhen.js, CootWorker.js and the data archives,
 * so a browser that had opened the viewer once kept the old worker, loader and
 * WASM for a year while the Moorhen React library -- which arrives in
 * content-hashed Next chunks -- was always current. The two then disagree about
 * their message protocol. Nothing fails loudly; only a hard refresh cures it.
 * (1.0.1-dev also raised the baked pthread pool from 8 to 32, and 8 is the
 * value that hung the viewer in July, so a cached browser kept the bad one.)
 *
 * Putting the Moorhen version in the path makes the promise true by
 * construction: a bump changes every URL, the old entries are never asked for
 * again and age out on their own.
 *
 * Electron is unaffected -- it serves the same files from disk at
 * /MoorhenAssets with no HTTP cache in between, so it keeps the bare prefix.
 */

/**
 * The Moorhen package version, injected at build time by next.config.ts.
 *
 * Read as a plain property access so the bundler can substitute the literal.
 * Empty when the build did not set it (a stale .env, an unusual entry point),
 * which is why callers fall back to the unversioned path rather than emitting
 * a URL with `undefined` in it.
 */
export const MOORHEN_VERSION = process.env.NEXT_PUBLIC_MOORHEN_VERSION || "";

/** Path segment that carries the version, e.g. `v/1.0.1-dev.g10d4c0b00/`. */
export function moorhenVersionSegment(version: string = MOORHEN_VERSION): string {
  // encodeURIComponent so an unexpected character in a version string cannot
  // alter the shape of the path the route then has to parse.
  return version ? `v/${encodeURIComponent(version)}/` : "";
}

/**
 * The urlPrefix to hand Moorhen for its runtime assets.
 *
 * @param isElectron serve from disk (no HTTP cache, so no version needed)
 */
export function moorhenUrlPrefix(isElectron: boolean): string {
  if (isElectron) return "/MoorhenAssets";
  return `/api/moorhen/${moorhenVersionSegment()}MoorhenAssets`;
}

/**
 * True when this looks like an Electron window.
 *
 * Kept here so both wrappers agree; `window` is guarded for SSR.
 */
export function isElectronWindow(): boolean {
  return typeof window !== "undefined" && !!(window as any).electronAPI;
}
