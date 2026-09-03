/**
 * An unconfigured CCP4i2 must not "find" an interpreter.
 *
 * When no CCP4 installation has been located, CCP4Dir is the empty string.
 * `path.join("", "bin", "ccp4-python")` is the RELATIVE path "bin/ccp4-python",
 * which existsSync happily resolves against the current working directory — so
 * launching the app from a directory that happens to contain bin/ccp4-python
 * would silently pick up an unrelated interpreter instead of asking the user
 * where CCP4 is.
 */
import { describe, it, expect, beforeAll, afterAll } from "vitest";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

import { findPython } from "../../main/ccp4i2-django-server";

let workdir: string;
let originalCwd: string;

beforeAll(() => {
  originalCwd = process.cwd();
  workdir = fs.mkdtempSync(path.join(os.tmpdir(), "ccp4i2-findpython-"));
  // A decoy: exactly what an empty CCP4Dir would resolve to, relative to cwd.
  const bin = path.join(workdir, "bin");
  fs.mkdirSync(bin, { recursive: true });
  const name = process.platform === "win32" ? "ccp4-python.bat" : "ccp4-python";
  fs.writeFileSync(path.join(bin, name), "#!/bin/sh\n");
  fs.chmodSync(path.join(bin, name), 0o755);
  process.chdir(workdir);
});

afterAll(() => {
  process.chdir(originalCwd);
  fs.rmSync(workdir, { recursive: true, force: true });
});

describe("findPython", () => {
  it("returns null when no CCP4 directory has been chosen", () => {
    // The decoy under cwd/bin must not be mistaken for a CCP4 installation.
    expect(findPython("", "")).toBeNull();
  });

  it("finds ccp4-python under a real CCP4 directory", () => {
    const found = findPython(workdir, "");
    expect(found).not.toBeNull();
    expect(found).toContain("bin");
  });

  it("returns null for a directory that is not a CCP4 installation", () => {
    const empty = fs.mkdtempSync(path.join(os.tmpdir(), "ccp4i2-nothing-"));
    try {
      expect(findPython(empty, "")).toBeNull();
    } finally {
      fs.rmSync(empty, { recursive: true, force: true });
    }
  });
});
