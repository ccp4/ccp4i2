/**
 * The guard that keeps a CCP4 9 from being offered, probed as "not installed",
 * and then "installed into" (2026-09-04 demo). Pure judgement is tested on
 * facts; the inspector is exercised against real CCP4 trees when this machine
 * has them.
 */
import { describe, expect, it } from "vitest";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { inspectPython, judgePython, listCcp4Dirs } from "../../main/ccp4i2-python-suitability";
import { CCP4I2_MIN_PYTHON } from "../../main/ccp4i2-server-version";

const facts = (version: [number, number, number] | null, writable: boolean | null = true) => ({
  version,
  sitePackages: "/x/site-packages",
  writable,
});

describe("judgePython", () => {
  it("a CCP4 9 (Python 3.9) is unsupported, and not installable", () => {
    const v = judgePython(facts([3, 9, 18]));
    expect(v.supported).toBe(false);
    expect(v.installable).toBe(false);
    expect(v.reason).toMatch(/3\.9/);
    expect(v.reason).toMatch(/CCP4 2026/);
  });

  it("the floor itself is supported", () => {
    const v = judgePython(facts([CCP4I2_MIN_PYTHON[0], CCP4I2_MIN_PYTHON[1], 0]));
    expect(v.supported).toBe(true);
    expect(v.installable).toBe(true);
  });

  it("a read-only site-packages is supported but not installable", () => {
    const v = judgePython(facts([3, 11, 14], false));
    expect(v.supported).toBe(true);
    expect(v.installable).toBe(false);
    expect(v.reason).toMatch(/not writable/);
    expect(v.reason).toMatch(/administrator/);
  });

  it("an interpreter that could not be run is unsupported", () => {
    const v = judgePython({ version: null, sitePackages: null, writable: null, error: "ENOENT" });
    expect(v.supported).toBe(false);
    expect(v.reason).toMatch(/Could not run/);
  });
});

describe("listCcp4Dirs", () => {
  it("orders by number, newest first, so a date beats a major version", () => {
    const root = fs.mkdtempSync(path.join(os.tmpdir(), "ccp4roots-"));
    for (const n of ["ccp4-9", "ccp4-20260702", "ccp4-8", "ccp4-20251105", "notccp4", "ccp4"]) {
      fs.mkdirSync(path.join(root, n));
    }
    expect(listCcp4Dirs(root)).toEqual(["ccp4-20260702", "ccp4-20251105", "ccp4-9", "ccp4-8"]);
    fs.rmSync(root, { recursive: true, force: true });
  });

  it("an unreadable root is simply empty", () => {
    expect(listCcp4Dirs("/no/such/root")).toEqual([]);
  });
});

// Against the real thing, where present on the developer's machine.
const ccp49 = "/Applications/ccp4-9/bin/ccp4-python";
const ccp426 = "/Users/nmemn/Developer/ccp4-20260702/bin/ccp4-python";

describe("inspectPython on real CCP4 trees", () => {
  it.skipIf(!fs.existsSync(ccp49))("CCP4 9 is judged unsupported", () => {
    const v = judgePython(inspectPython(ccp49));
    expect(v.facts.version?.slice(0, 2)).toEqual([3, 9]);
    expect(v.supported).toBe(false);
  });

  it.skipIf(!fs.existsSync(ccp426))("CCP4 2026 is judged supported and installable", () => {
    const v = judgePython(inspectPython(ccp426));
    expect(v.facts.version?.[0]).toBe(3);
    expect(v.supported).toBe(true);
    expect(v.installable).toBe(true);
  });

  it("a missing interpreter reports an error rather than throwing", () => {
    const f = inspectPython("/no/such/python");
    expect(f.version).toBeNull();
    expect(f.error).toBeTruthy();
  });
});
