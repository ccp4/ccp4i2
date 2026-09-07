/**
 * The packaged desktop app has no next.config.ts on disk, so the Electron
 * main process hands Next the config recorded by `next build` through the
 * same environment variable Next's own standalone server uses. Without it the
 * middleware body-clone cap falls back to 10 MB and larger uploads are
 * silently truncated.
 */
import { describe, expect, it } from "vitest";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import {
  applyBuildTimeNextConfig,
  readBuildTimeNextConfig,
  REQUIRED_SERVER_FILES,
  STANDALONE_CONFIG_ENV,
} from "../../main/ccp4i2-next-config";

function rendererDirWithManifest(manifest: unknown): string {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), "ccp4i2-next-config-"));
  const file = path.join(dir, REQUIRED_SERVER_FILES);
  fs.mkdirSync(path.dirname(file), { recursive: true });
  fs.writeFileSync(file, typeof manifest === "string" ? manifest : JSON.stringify(manifest));
  return dir;
}

const manifest = {
  version: 1,
  config: {
    trailingSlash: true,
    experimental: { middlewareClientMaxBodySize: 2147483648 },
  },
  files: [],
};

describe("readBuildTimeNextConfig", () => {
  it("returns the config block of the build manifest", () => {
    const dir = rendererDirWithManifest(manifest);
    expect(readBuildTimeNextConfig(dir)).toEqual(manifest.config);
  });

  it("returns null when the manifest is absent, unreadable or has no config", () => {
    expect(readBuildTimeNextConfig(path.join(os.tmpdir(), "no-such-renderer"))).toBeNull();
    expect(readBuildTimeNextConfig(rendererDirWithManifest("not json"))).toBeNull();
    expect(readBuildTimeNextConfig(rendererDirWithManifest({ version: 1 }))).toBeNull();
  });
});

describe("applyBuildTimeNextConfig", () => {
  it("sets the standalone-config variable Next reads instead of a config file", () => {
    const env = {} as NodeJS.ProcessEnv;
    expect(applyBuildTimeNextConfig(rendererDirWithManifest(manifest), env)).toBe(true);
    expect(JSON.parse(env[STANDALONE_CONFIG_ENV] as string)).toEqual(manifest.config);
  });

  it("leaves a variable an operator has already set alone", () => {
    const env = { [STANDALONE_CONFIG_ENV]: '{"trailingSlash":false}' } as unknown as NodeJS.ProcessEnv;
    expect(applyBuildTimeNextConfig(rendererDirWithManifest(manifest), env)).toBe(true);
    expect(env[STANDALONE_CONFIG_ENV]).toBe('{"trailingSlash":false}');
  });

  it("reports false and sets nothing when there is no manifest", () => {
    const env = {} as NodeJS.ProcessEnv;
    expect(applyBuildTimeNextConfig(path.join(os.tmpdir(), "no-such-renderer"), env)).toBe(false);
    expect(env[STANDALONE_CONFIG_ENV]).toBeUndefined();
  });
});
