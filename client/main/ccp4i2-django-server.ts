import { ccp4_setup_sh } from "./ccp4i2-setup-sh";
import { ccp4_setup_windows } from "./ccp4i2-setup-windows";
import path from "node:path";
import { spawn } from "node:child_process";
import { fileURLToPath } from "node:url";
import { execSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import { ccp4i2Home } from "./ccp4i2-preferences";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

/**
 * Gets the next available run number for log files
 */
function getNextRunNumber(logDir: string): number {
  if (!fs.existsSync(logDir)) {
    fs.mkdirSync(logDir, { recursive: true });
    return 1;
  }

  const files = fs.readdirSync(logDir);
  const logFiles = files.filter((file) => file.match(/^uvicorn-(\d+)\.log$/));

  if (logFiles.length === 0) {
    return 1;
  }

  const runNumbers = logFiles.map((file) => {
    const match = file.match(/^uvicorn-(\d+)\.log$/);
    return match ? parseInt(match[1], 10) : 0;
  });

  return Math.max(...runNumbers) + 1;
}

/**
 * Finds the Python executable, preferring ccp4-python from the CCP4 installation.
 * Falls back to project's virtual environment if ccp4-python is not available.
 *
 * @param CCP4Dir - The path to the CCP4 installation directory.
 * @param projectRoot - The path to the project root (fallback for .venv in dev mode).
 * @returns The path to the Python executable, or null if not found.
 */
export function findPython(CCP4Dir: string, projectRoot: string): string | null {
  const isWindows = process.platform === "win32";

  // Prefer ccp4-python from CCP4 installation. Skipped when no CCP4 has been
  // chosen: path.join("", "bin", ...) yields the RELATIVE "bin/ccp4-python",
  // which existsSync resolves against the current working directory — so an
  // unconfigured app could "find" an unrelated interpreter under wherever it
  // happened to be launched from.
  if (CCP4Dir && CCP4Dir.length > 0) {
    const ccp4PythonBin = isWindows ? "ccp4-python.bat" : "ccp4-python";
    const ccp4PythonPath = path.join(CCP4Dir, "bin", ccp4PythonBin);
    if (fs.existsSync(ccp4PythonPath)) {
      return ccp4PythonPath;
    }
  }

  // Fallback to project's virtual environment (dev mode only)
  if (projectRoot) {
    const pythonBin = isWindows ? "python.exe" : "python";
    const binDir = isWindows ? "Scripts" : "bin";

    const venvPaths = [
      path.join(projectRoot, ".venv", binDir, pythonBin),
      path.join(projectRoot, ".venv.py311", binDir, pythonBin),
      path.join(projectRoot, ".venv.py39", binDir, pythonBin),
    ];

    for (const venvPath of venvPaths) {
      if (fs.existsSync(venvPath)) {
        return venvPath;
      }
    }
  }

  return null;
}

/**
 * Starts the Django server using Uvicorn with the specified configuration.
 *
 * In packaged mode, ccp4i2 is pip-installed into ccp4-python — no bundled
 * Python code in the Electron resources. Uvicorn loads the ASGI app via
 * the module path ccp4i2.config.asgi:application.
 *
 * In dev mode, the server runs from the local server/ directory with
 * ccp4-python (which has ccp4i2 installed via pip install -e .).
 *
 * @param CCP4Dir - The path to the CCP4 directory (provides ccp4-python and env vars like $CLIBD).
 * @param projectRoot - The path to the ccp4i2 project root (dev mode .venv fallback).
 * @param UVICORN_PORT - The port number for the Uvicorn server.
 * @param NEXT_PORT - The port number for the Next.js server.
 * @param isDev - A boolean flag indicating whether the server is running in development mode.
 * @param CCP4I2_PROJECTS_DIR - The directory where CCP4i2 projects are stored.
 * @returns The spawned Python process running the Uvicorn server.
 */
export async function startDjangoServer(
  CCP4Dir: string,
  projectRoot: string,
  UVICORN_PORT: number,
  NEXT_PORT: number,
  isDev: boolean = false,
  CCP4I2_PROJECTS_DIR: string = ""
) {
  // Set up CCP4 environment variables (still needed for $CLIBD, $CBIN, etc.).
  // Guarded: with no CCP4 located, CCP4Dir is "" and the setup helpers would
  // derive $CCP4="" and $CCP4I_TCLTK="/bin" from it, i.e. point the environment
  // at the filesystem root. Better to leave the environment alone and fail on
  // the interpreter below, which can say something useful.
  if (CCP4Dir && CCP4Dir.length > 0) {
    if (process.platform === "win32") {
      ccp4_setup_windows(CCP4Dir);
    } else {
      ccp4_setup_sh(CCP4Dir);
    }
  }

  // Find Python interpreter (prefers ccp4-python, falls back to .venv)
  const PYTHON_PATH = findPython(CCP4Dir, projectRoot);
  if (!PYTHON_PATH) {
    throw new Error(
      CCP4Dir && CCP4Dir.length > 0
        ? `Could not find ccp4-python in ${CCP4Dir}/bin. Check that this is a ` +
          `CCP4 installation root, and choose a different one on the setup page ` +
          `if not.`
        : `No CCP4 installation has been chosen yet. Use "Locate" on the setup ` +
          `page to point CCP4i2 at your CCP4 directory (the one containing bin/).`
    );
  }

  if (CCP4I2_PROJECTS_DIR.length > 0) {
    process.env.CCP4I2_PROJECTS_DIR = CCP4I2_PROJECTS_DIR;
  }
  // Deliberately NOT setting CCP4I2_DB_FILE from the projects directory.
  //
  // Doing so tied the database to wherever projects happened to live, so
  // choosing a new projects folder silently produced a NEW, EMPTY database and
  // the previous one simply stopped being consulted — projects appearing to
  // vanish, with nothing lost but nothing visible either. Qt-era CCP4i2 kept
  // one database in the user's CCP4i2 directory regardless of where projects
  // sat, and for alpha and beta we match that: the server default is
  // <ccp4i2 home>/db.sqlite3.
  //
  // An installation that already has an explicit "database" in preferences.json
  // keeps it — that file is honoured ahead of the default, so existing testers
  // carry on using the database they have rather than being switched to an
  // empty one.
  //
  // The longer-term intent is different again: self-contained SETS of projects,
  // each with its own database alongside them, chosen explicitly rather than as
  // a side effect of picking a folder.

  console.log(`🚀 Next.js running on http://localhost:${NEXT_PORT}`);
  console.log(`🐍 Using Python: ${PYTHON_PATH}`);

  // In dev mode, use server/ as the working directory
  const serverCwd = isDev
    ? path.join(process.cwd(), "..", "server")
    : undefined;

  // Typed as ProcessEnv, not Record<string, string | undefined>: Next's global
  // augmentation makes NODE_ENV a required member of ProcessEnv, so the looser
  // type does not satisfy the env option of spawn/execSync. That mismatch made
  // four call overloads fail and their results widen to `never`, producing
  // seven type errors in this file — invisible until something typechecked it.
  const pythonEnv: NodeJS.ProcessEnv = {
    ...process.env,
    DJANGO_SETTINGS_MODULE: "ccp4i2.config.settings",
    UVICORN_PORT: `${UVICORN_PORT}`,
    NEXT_ADDRESS: `http://localhost:${NEXT_PORT}`,
    // Force local execution mode for Electron app
    EXECUTION_MODE: "local",
    MPLBACKEND: "Agg", // Force matplotlib to use non-GUI backend
    // Windows-specific DLL path fixes
    ...(process.platform === "win32" && {
      PATH: `${path.join(CCP4Dir, "bin")};${process.env.PATH}`,
      CCP4: CCP4Dir,
      // Force matplotlib to avoid Qt completely
      MPLCONFIGDIR: path.join(os.tmpdir(), "matplotlib-config"),
    }),
  };

  // Run migrations
  try {
    execSync(`"${PYTHON_PATH}" -m django migrate`, {
      env: pythonEnv,
      stdio: "inherit",
      ...(serverCwd && { cwd: serverCwd }),
    });
    console.log(`🐍 Migration completed successfully`);
  } catch (error) {
    console.error(`🐍 Migration failed:`, error);
    // Try alternative approach for Windows
    if (process.platform === "win32") {
      try {
        execSync(
          `"${PYTHON_PATH}" -c "import os; os.environ['MPLBACKEND']='Agg'; import django; django.setup(); from django.core.management import execute_from_command_line; execute_from_command_line(['manage.py', 'migrate'])"`,
          {
            env: pythonEnv,
            ...(serverCwd && { cwd: serverCwd }),
          }
        );
        console.log(`🐍 Alternative migration completed`);
      } catch (altError) {
        console.error(`🐍 Alternative migration also failed:`, altError);
        throw altError;
      }
    } else {
      throw error;
    }
  }

  // Write the on-disk recovery snapshot for any project that lacks one.
  //
  // Snapshots are normally kept current by database signals, but a project
  // that predates the mechanism never triggers one -- and a finished, dormant
  // project is both the kind most worth protecting and the kind nothing will
  // ever touch again. --missing-only makes this a no-op once every project is
  // covered, so it is cheap to run on every launch.
  //
  // Never fatal: failing to write a safety net must not stop the app starting.
  try {
    execSync(`"${PYTHON_PATH}" -m django snapshot_projects`, {
      env: pythonEnv,
      stdio: "inherit",
      ...(serverCwd && { cwd: serverCwd }),
    });
  } catch (error) {
    console.error(`🐍 Could not write project recovery snapshots:`, error);
  }

  // Setup logging for production
  let logStream: fs.WriteStream | null = null;
  if (!isDev) {
    // Under the one user home, not a ~/.ccp4i2 of its own: server logs are
    // per-user state like everything else, and a tester asked to "send me your
    // log" should not have to be told which of several directories to look in.
    const logDir = path.join(ccp4i2Home(), "logs");

    const runNumber = getNextRunNumber(logDir);
    const logFile = path.join(logDir, `uvicorn-${runNumber}.log`);
    fs.mkdirSync(logDir, { recursive: true });
    logStream = fs.createWriteStream(logFile, { flags: "a" });
    console.log(`🐍 Uvicorn logs will be written to: ${logFile}`);
  }

  // Start Python process — ccp4i2 is pip-installed so uvicorn uses module path
  // Use 2 workers for concurrent requests (no --reload, requires manual restart)
  //
  // --ws none: the ccp4i2 backend is a plain HTTP Django ASGI app
  // (get_asgi_application; no channels / websocket routing). uvicorn's default
  // "--ws auto" eagerly imports its websockets protocol impl at startup, which
  // on newer uvicorn (>=0.30) is the sansio impl that needs websockets>=13.
  // If the environment ships an older websockets (e.g. the CCP4 bundle's 10.4),
  // that import raises "cannot import name 'ServerProtocol'" and every worker
  // crashes on boot. Since we never serve websockets, disable the protocol
  // entirely rather than depend on a specific websockets version.
  const uvicornArgs = [
    "-m", "uvicorn", "ccp4i2.config.asgi:application", "--workers", "2",
    "--ws", "none",
    "--log-level", "warning",  // Suppress per-request INFO access logs
  ];

  // shell:true is REQUIRED on Windows to launch ccp4-python.bat (a .bat cannot be
  // spawned without a shell on modern Node — see spawnPython in ccp4i2-ipc.ts).
  // With a shell the executable is NOT auto-escaped, so quote it when the CCP4
  // install path contains whitespace — matching the quoted `migrate` execSync
  // above, which would otherwise succeed while this launch broke. uvicornArgs are
  // all space-free literals and need no quoting.
  const launcher = /\s/.test(PYTHON_PATH) ? `"${PYTHON_PATH}"` : PYTHON_PATH;
  const pythonProcess = spawn(launcher, uvicornArgs, {
    env: pythonEnv,
    shell: true,
    ...(serverCwd && { cwd: serverCwd }),
  });

  console.log(`🚀 Uvicorn running on http://localhost:${UVICORN_PORT}`);

  pythonProcess.stdout.on("data", (data) => {
    if (isDev) {
      console.log(`🐍 Python Output: ${data}`);
    } else {
      logStream?.write(`[STDOUT] ${new Date().toISOString()}: ${data}`);
    }
  });

  pythonProcess.stderr.on("data", (data) => {
    if (isDev) {
      // Uvicorn writes all access logs and Django warnings to stderr.
      // Only label genuine errors (tracebacks, ERROR level) as errors;
      // routine output uses console.log to avoid alarming users.
      const text = String(data);
      const isError = /\bERROR\b|Traceback|Exception/i.test(text);
      if (isError) {
        console.error(`🐍 Python Error: ${data}`);
      } else {
        console.log(`🐍 Python: ${data}`);
      }
    } else {
      logStream?.write(`[STDERR] ${new Date().toISOString()}: ${data}`);
    }
  });

  pythonProcess.on("close", (code) => {
    const message = `🐍 Python process exited with code ${code}`;
    if (isDev) {
      console.log(message);
    } else {
      logStream?.write(`[EXIT] ${new Date().toISOString()}: ${message}\n`);
      logStream?.end();
    }
  });

  return pythonProcess;
}
