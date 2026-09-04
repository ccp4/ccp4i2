/**
 * Can this ccp4-python host the CCP4i2 backend at all?
 *
 * The readiness probe asks whether ccp4i2 is installed and importable, and
 * offers Install when it is not. That is the right question for a CCP4 that
 * *could* hold the backend, and the wrong one for a CCP4 9: its Python is 3.9,
 * its site-packages is root-owned, and the "ccp4i2" already in it is the
 * Qt-era package. Pointed there, the app said "not installed", offered
 * Install, and the install could only fail — which is what happened in a demo
 * on 2026-09-04, after the user cleared their settings to show a clean start
 * and the default fell back to /Applications/ccp4-9.
 *
 * Two facts decide it, and they are kept apart because they mean different
 * things:
 *
 *   - the Python version, against CCP4I2_MIN_PYTHON. Below the floor nothing
 *     will work, installed or not: the environment is UNSUPPORTED and the
 *     only remedy is a different CCP4.
 *   - whether site-packages is writable. A read-only one cannot take an
 *     install, but a backend an administrator put there for everyone runs
 *     fine — so this gates Install, not Launch.
 */
import { execFileSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { CCP4I2_MIN_PYTHON } from "./ccp4i2-server-version";

export interface PythonFacts {
  /** sys.version_info[:3], or null when the interpreter could not be run. */
  version: [number, number, number] | null;
  sitePackages: string | null;
  /** os.access(site-packages, W_OK); null when unknown. */
  writable: boolean | null;
  /** Why the interpreter could not be inspected, when it could not. */
  error?: string;
}

export interface Suitability {
  /** Above the Python floor: the backend can run here (if installed). */
  supported: boolean;
  /** Supported AND site-packages writable: Install can succeed here. */
  installable: boolean;
  /** One sentence for the user when supported or installable is false. */
  reason?: string;
  facts: PythonFacts;
}

const PROBE_SOURCE = [
  "import json, os, site, sys",
  "sp = None",
  "try:",
  "    sp = (site.getsitepackages() or [None])[0]",
  "except Exception:",
  "    pass",
  "print(json.dumps({",
  '    "version": list(sys.version_info[:3]),',
  '    "sitePackages": sp,',
  '    "writable": bool(sp and os.access(sp, os.W_OK)),',
  "}))",
  "",
].join("\n");

/** Run the interpreter once and report what it is. Synchronous and short. */
export function inspectPython(pythonPath: string, timeoutMs = 10000): PythonFacts {
  const probePath = path.join(
    os.tmpdir(),
    `ccp4i2-pyfacts-${process.pid}-${Date.now()}.py`
  );
  try {
    fs.writeFileSync(probePath, PROBE_SOURCE);
    // ccp4-python.bat on Windows needs a shell; quote for cmd.exe the way the
    // other spawns in this app do. Elsewhere, no shell.
    const win = process.platform === "win32";
    const q = (s: string) => (win && /\s/.test(s) ? `"${s}"` : s);
    const out = execFileSync(q(pythonPath), [q(probePath)], {
      encoding: "utf8",
      timeout: timeoutMs,
      shell: win,
      stdio: ["ignore", "pipe", "pipe"],
      env: { ...process.env, MPLBACKEND: "Agg" },
    });
    const info = JSON.parse(out.trim().split("\n").pop() || "{}");
    const v = Array.isArray(info.version) ? info.version : null;
    return {
      version: v && v.length >= 2 ? [v[0], v[1], v[2] ?? 0] : null,
      sitePackages: typeof info.sitePackages === "string" ? info.sitePackages : null,
      writable: typeof info.writable === "boolean" ? info.writable : null,
    };
  } catch (error) {
    return {
      version: null,
      sitePackages: null,
      writable: null,
      error: error instanceof Error ? error.message : String(error),
    };
  } finally {
    try {
      fs.unlinkSync(probePath);
    } catch {
      /* best-effort */
    }
  }
}

const fmt = (v: readonly number[]) => v.slice(0, 2).join(".");

/** Pure: judge facts against the floor. Separated so it can be tested. */
export function judgePython(
  facts: PythonFacts,
  minPython: readonly [number, number] = CCP4I2_MIN_PYTHON
): Suitability {
  if (!facts.version) {
    return {
      supported: false,
      installable: false,
      reason:
        `Could not run this CCP4's Python` +
        (facts.error ? ` (${facts.error.split("\n")[0]})` : "") +
        `. Choose a CCP4 installation whose bin/ccp4-python runs.`,
      facts,
    };
  }
  const [maj, min] = facts.version;
  const below = maj < minPython[0] || (maj === minPython[0] && min < minPython[1]);
  if (below) {
    return {
      supported: false,
      installable: false,
      reason:
        `This CCP4's Python is ${fmt(facts.version)}; CCP4i2 needs ${fmt(minPython)} ` +
        `or later. CCP4 9 and earlier ship Python 3.9 and cannot host it — ` +
        `choose a CCP4 2026 (or later) installation.`,
      facts,
    };
  }
  if (facts.writable === false) {
    return {
      supported: true,
      installable: false,
      reason:
        `This CCP4's Python packages directory is not writable by you` +
        (facts.sitePackages ? ` (${facts.sitePackages})` : "") +
        `, so CCP4i2 cannot be installed into it. If an administrator has ` +
        `installed it for everyone it will still run; otherwise choose a CCP4 ` +
        `installation you own.`,
      facts,
    };
  }
  return { supported: true, installable: true, facts };
}

/** Inspect and judge in one call. */
export function assessPython(pythonPath: string): Suitability {
  return judgePython(inspectPython(pythonPath));
}

/** ccp4-* directories directly under `root`, newest first by their number. */
export const listCcp4Dirs = (root: string): string[] => {
  let names: string[] = [];
  try {
    names = fs.readdirSync(root).filter((n) => /^ccp4-\d+/.test(n));
  } catch {
    return [];
  }
  const num = (n: string) => parseInt(n.replace(/^ccp4-/, ""), 10) || 0;
  return names.sort((a, b) => num(b) - num(a));
};

