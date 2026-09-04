/**
 * Quitting must take the whole server tree with it, not one process of it.
 * Exercised on a real tree: a shell that starts two sleepers, spawned the way
 * the app now spawns uvicorn (its own process group).
 */
import { describe, expect, it } from "vitest";
import { spawn, execFileSync } from "node:child_process";
import {
  parentPidEnvironment,
  processGroupOptions,
  terminateProcessTree,
} from "../../main/ccp4i2-process-tree";

const posix = process.platform !== "win32";
const marker = `ccp4i2-tree-test-${process.pid}`;

// pgrep called directly, not through a shell: a `sh -c "pgrep -f <marker>"`
// wrapper carries the marker in its own command line and pgrep reports it,
// which made this read "alive" on Linux after the tree was dead.
const treeAlive = () => {
  try {
    return execFileSync("pgrep", ["-f", marker], { encoding: "utf8" }).trim().length > 0;
  } catch {
    return false; // pgrep exits 1 when nothing matches
  }
};

describe("terminateProcessTree", () => {
  it.skipIf(!posix)("kills the leader and every descendant", async () => {
    // A leader that starts two children and waits, all tagged so pgrep can
    // find whatever survives.
    const child = spawn(
      "sh",
      ["-c", `sh -c "sleep 60" ${marker} & sh -c "sleep 60" ${marker} & wait`],
      { stdio: "ignore", ...processGroupOptions() }
    );
    await new Promise((r) => setTimeout(r, 300));
    expect(treeAlive()).toBe(true);

    terminateProcessTree(child, 2000);
    await new Promise((r) => setTimeout(r, 300));

    expect(treeAlive()).toBe(false);
  });

  it("tolerates a child that has already gone, and no child at all", () => {
    expect(() => terminateProcessTree(null)).not.toThrow();
    const done = spawn("sh", ["-c", "true"], { stdio: "ignore", ...processGroupOptions() });
    return new Promise<void>((resolve) =>
      done.on("exit", () => {
        expect(() => terminateProcessTree(done)).not.toThrow();
        resolve();
      })
    );
  });
});

describe("what the server is told", () => {
  it("carries our pid so the server can watch us", () => {
    expect(parentPidEnvironment()).toEqual({ CCP4I2_PARENT_PID: String(process.pid) });
  });

  it.skipIf(!posix)("puts the server in its own process group on POSIX", () => {
    expect(processGroupOptions()).toEqual({ detached: true });
  });
});
