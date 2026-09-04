/**
 * Stop a server we started, and everything it started, on every way out.
 *
 * uvicorn with --workers 2 is a tree: a supervisor and two workers (and on
 * Windows, the .bat's shell above them). `child.kill()` signals one process.
 * The app's quit handler did that and nothing else, and it ran only on a
 * clean quit — never on a crash, a debugger stop, or SIGKILL. A Linux tester
 * accumulated a stack of orphaned backends this way, each bound to the next
 * free port, eating memory (tracker: "Zombie uvicorn processes").
 *
 * Here: POSIX children are spawned into their own process group so the
 * whole tree can be signalled at once; termination is SIGTERM, then SIGKILL
 * after a grace period if anything is still there; Windows uses taskkill /T.
 * `registerExitHandlers` hooks the paths a quit handler misses. The
 * crash-and-SIGKILL case, where no handler can run, is covered from the
 * other side by the Python parent watchdog (ccp4i2.config.parent_watchdog),
 * which the server starts when CCP4I2_PARENT_PID is in its environment.
 */
import { ChildProcess, execFileSync } from "node:child_process";

const isWindows = process.platform === "win32";

/** Spawn options that put a POSIX child in its own process group. */
export const processGroupOptions = (): { detached?: boolean } =>
  isWindows ? {} : { detached: true };

/** Environment the server needs to watch us. */
export const parentPidEnvironment = (): Record<string, string> => ({
  CCP4I2_PARENT_PID: String(process.pid),
});

const signalGroup = (pid: number, signal: NodeJS.Signals): boolean => {
  try {
    // A negative pid addresses the process group the child leads.
    process.kill(-pid, signal);
    return true;
  } catch (err) {
    if ((err as NodeJS.ErrnoException).code === "ESRCH") return false; // gone
    // Not a group leader (spawned without detached): fall back to the pid.
    try {
      process.kill(pid, signal);
      return true;
    } catch {
      return false;
    }
  }
};

const groupIsAlive = (pid: number): boolean => signalGroup(pid, 0 as unknown as NodeJS.Signals);

/**
 * Terminate `child` and its descendants. Synchronous by design: the callers
 * are exit paths, where nothing asynchronous is guaranteed to run.
 */
export function terminateProcessTree(
  child: ChildProcess | null | undefined,
  graceMs = 1500
): void {
  if (!child || child.pid == null) return;
  const pid = child.pid;

  if (isWindows) {
    try {
      execFileSync("taskkill", ["/PID", String(pid), "/T", "/F"], {
        stdio: "ignore",
        timeout: 10000,
      });
    } catch {
      /* already gone, or taskkill unavailable: nothing more we can do */
    }
    return;
  }

  if (!signalGroup(pid, "SIGTERM")) return;
  // Give uvicorn its orderly shutdown, then insist. Busy-wait: this runs
  // inside exit handlers where timers will not fire.
  const deadline = Date.now() + graceMs;
  while (Date.now() < deadline) {
    if (!groupIsAlive(pid)) return;
    Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0, 50);
  }
  if (groupIsAlive(pid)) signalGroup(pid, "SIGKILL");
}

/**
 * Run `cleanup` on every way out of the process that still runs JavaScript:
 * a normal exit, the signals a terminal or a debugger sends, and an uncaught
 * exception. Each fires `cleanup` at most once.
 */
export function registerExitHandlers(cleanup: () => void): void {
  let done = false;
  const once = () => {
    if (done) return;
    done = true;
    try {
      cleanup();
    } catch {
      /* exiting anyway */
    }
  };
  process.on("exit", once);
  for (const sig of ["SIGINT", "SIGTERM", "SIGHUP"] as NodeJS.Signals[]) {
    process.on(sig, () => {
      once();
      process.exit(0);
    });
  }
  process.on("uncaughtException", (err) => {
    console.error("Uncaught exception; shutting the backend down:", err);
    once();
    process.exit(1);
  });
}
