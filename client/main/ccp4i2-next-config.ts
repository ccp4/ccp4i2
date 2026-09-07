import fs from "node:fs";
import path from "node:path";

/**
 * Hand the packaged Next server the configuration it was built with.
 *
 * The desktop app ships `renderer/.next` and `renderer/public` but not
 * `renderer/next.config.ts` (see "build.files" in package.json), and a
 * TypeScript config could not be loaded at runtime anyway. So when the app is
 * packaged, `next()` finds no config file and runs on Next's defaults. Most
 * settings are baked into the build manifests and survive that; the ones read
 * at request time do not. The one that bit was
 * `experimental.middlewareClientMaxBodySize`: every request that passes
 * through middleware has its body cloned with a size cap, the default cap is
 * 10 MB, and a body over the cap is silently truncated - so importing a
 * project zip larger than 10 MB failed with a bare "fetch failed" from the
 * proxy, whatever next.config.ts said.
 *
 * `next()`'s `conf` option is not a way out: in Next 15 the custom-server
 * entry point starts the router server, which loads its own config from disk
 * and ignores `conf`. What the router server does honour is the
 * `__NEXT_PRIVATE_STANDALONE_CONFIG` environment variable - the mechanism
 * Next's own generated standalone `server.js` uses to pass the build-time
 * config to a server that has no config file. `next build` records that
 * config, defaults applied and sizes normalised, in
 * `.next/required-server-files.json`; this reads it back and sets the
 * variable before Next starts.
 *
 * Development is untouched: there the config file is on disk and Next loads
 * it as usual.
 */
export const REQUIRED_SERVER_FILES = path.join(".next", "required-server-files.json");
export const STANDALONE_CONFIG_ENV = "__NEXT_PRIVATE_STANDALONE_CONFIG";

export function readBuildTimeNextConfig(rendererDir: string): Record<string, unknown> | null {
  const manifestPath = path.join(rendererDir, REQUIRED_SERVER_FILES);
  let raw: string;
  try {
    raw = fs.readFileSync(manifestPath, "utf8");
  } catch {
    return null;
  }
  try {
    const manifest = JSON.parse(raw);
    return manifest && typeof manifest.config === "object" && manifest.config !== null
      ? (manifest.config as Record<string, unknown>)
      : null;
  } catch {
    return null;
  }
}

/**
 * Sets the standalone-config variable from the build manifest, unless the
 * caller (or an operator) has already set it. Returns true when the variable
 * is set on exit, so the caller can log what the server will run with.
 */
export function applyBuildTimeNextConfig(
  rendererDir: string,
  env: NodeJS.ProcessEnv = process.env
): boolean {
  if (env[STANDALONE_CONFIG_ENV]) return true;
  const config = readBuildTimeNextConfig(rendererDir);
  if (!config) return false;
  env[STANDALONE_CONFIG_ENV] = JSON.stringify(config);
  return true;
}
