/**
 * The Moorhen asset URL carries the Moorhen version.
 *
 * The route serving these files answers `immutable` for a year, which it must:
 * one window starts 33 pthread workers, each loading moorhen.js. But the URL
 * used to be fixed, so after a Moorhen bump a returning browser kept the old
 * worker, loader and WASM against the new React library and they disagreed
 * about their message protocol -- silently, curable only by a hard refresh.
 */
import { describe, expect, it } from "vitest";
import { moorhenUrlPrefix, moorhenVersionSegment } from "../lib/moorhen-asset-path";

describe("moorhenVersionSegment", () => {
  it("builds a path segment from a version", () => {
    expect(moorhenVersionSegment("1.0.1-dev.g10d4c0b00")).toBe("v/1.0.1-dev.g10d4c0b00/");
  });

  it("is empty when the build injected no version", () => {
    // Falls back to the unversioned path rather than emitting `v/undefined/`.
    expect(moorhenVersionSegment("")).toBe("");
  });

  it("escapes a version that could otherwise reshape the path", () => {
    expect(moorhenVersionSegment("1.0/../../etc")).not.toContain("../");
  });
});

describe("moorhenUrlPrefix", () => {
  it("versions the web path so each Moorhen bump busts the immutable cache", () => {
    const a = moorhenUrlPrefix(false);
    expect(a).toMatch(/^\/api\/moorhen\//);
    expect(a).toMatch(/MoorhenAssets$/);
  });

  it("leaves Electron on the bare on-disk path", () => {
    // Served from disk with no HTTP cache in between, so nothing to bust.
    expect(moorhenUrlPrefix(true)).toBe("/MoorhenAssets");
  });

  it("gives the two builds different prefixes", () => {
    expect(moorhenUrlPrefix(true)).not.toBe(moorhenUrlPrefix(false));
  });
});

/**
 * The route strips `v/<version>/` before its allow-list runs. These cases pin
 * that the stripping cannot be used to smuggle a path past the allow-list --
 * the check and the file read both use the stripped value.
 */
describe("version-segment stripping (as the route does it)", () => {
  const strip = (requested: string) => {
    const versioned = requested.match(/^v\/([^/]+)\/(.+)$/);
    return versioned ? versioned[2] : requested;
  };

  it("removes the version and keeps the asset path", () => {
    expect(strip("v/1.0.1-dev.g10d4c0b00/MoorhenAssets/data.tar.gz")).toBe(
      "MoorhenAssets/data.tar.gz"
    );
    expect(strip("v/1.0.1/moorhen.wasm")).toBe("moorhen.wasm");
  });

  it("leaves an unversioned path alone, so older clients keep working", () => {
    expect(strip("MoorhenAssets/data.tar.gz")).toBe("MoorhenAssets/data.tar.gz");
    expect(strip("moorhen.wasm")).toBe("moorhen.wasm");
  });

  it("consumes only one segment, so a traversal still faces the allow-list", () => {
    // The point is that what survives stripping is what gets checked; it is
    // not silently treated as already-validated.
    expect(strip("v/1.0.1/../../../etc/passwd")).toBe("../../../etc/passwd");
    expect(strip("v/1.0.1/../secrets.env")).toBe("../secrets.env");
  });

  it("does not treat a bare 'v' directory as a version", () => {
    expect(strip("v/onlyonesegment")).toBe("v/onlyonesegment");
  });
});

/**
 * The version has to survive the *build*, not just the source.
 *
 * The first cut of this read `require("moorhen/package.json")` in
 * next.config.ts. Moorhen's package has an `exports` map that does not expose
 * ./package.json, so that threw ERR_PACKAGE_PATH_NOT_EXPORTED, a defensive
 * catch turned it into "", and the built bundle quietly emitted the old
 * unversioned path -- the fix compiled, passed its unit tests, and did
 * nothing. This checks the resolution strategy the config actually uses.
 */
describe("build-time version resolution", () => {
  it("reads Moorhen's package.json off disk, since its exports map hides it", async () => {
    const { existsSync, readFileSync } = await import("fs");
    const nodePath = await import("path");

    let dir = nodePath.resolve(__dirname, "..");
    let found = "";
    for (let up = 0; up < 5; up++) {
      const candidate = nodePath.join(dir, "node_modules", "moorhen", "package.json");
      if (existsSync(candidate)) {
        found = JSON.parse(readFileSync(candidate, "utf8")).version;
        break;
      }
      const parent = nodePath.dirname(dir);
      if (parent === dir) break;
      dir = parent;
    }

    expect(found, "moorhen package.json should be reachable by walking up").toBeTruthy();
    expect(found).toMatch(/\d+\.\d+\.\d+/);
  });

  it("cannot be reached through the package exports map", async () => {
    // Pins *why* the walk exists. If Moorhen ever exports ./package.json this
    // fails, and the simpler require() becomes available.
    const { createRequire } = await import("module");
    const requireFrom = createRequire(__filename);
    expect(() => requireFrom.resolve("moorhen/package.json")).toThrow();
  });
});
