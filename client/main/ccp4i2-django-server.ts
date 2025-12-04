import { ccp4_setup_sh } from "./ccp4i2-setup-sh";
import { ccp4_setup_windows } from "./ccp4i2-setup-windows";
import path from "node:path";
import { spawn } from "node:child_process";
import { fileURLToPath } from "node:url";
import { execSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";

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
 * @param projectRoot - The path to the project root (fallback for .venv).
 * @returns The path to the Python executable, or null if not found.
 */
function findPython(CCP4Dir: string, projectRoot: string): string | null {
  const isWindows = process.platform === "win32";

  // Prefer ccp4-python from CCP4 installation
  const ccp4PythonBin = isWindows ? "ccp4-python.bat" : "ccp4-python";
  const ccp4PythonPath = path.join(CCP4Dir, "bin", ccp4PythonBin);
  if (fs.existsSync(ccp4PythonPath)) {
    return ccp4PythonPath;
  }

  // Fallback to project's virtual environment
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

  return null;
}

/**
 * Starts the Django server using Uvicorn with the specified configuration.
 * Prefers ccp4-python from the CCP4 installation, falls back to .venv if not available.
 *
 * @param CCP4Dir - The path to the CCP4 directory (provides ccp4-python and env vars like $CLIBD).
 * @param projectRoot - The path to the ccp4i2 project root (fallback .venv location).
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
  // Set up CCP4 environment variables (still needed for $CLIBD, $CBIN, etc.)
  if (process.platform === "win32") {
    ccp4_setup_windows(CCP4Dir);
  } else {
    ccp4_setup_sh(CCP4Dir);
  }

  // Find Python interpreter (prefers ccp4-python, falls back to .venv)
  const PYTHON_PATH = findPython(CCP4Dir, projectRoot);
  if (!PYTHON_PATH) {
    throw new Error(
      `Could not find Python interpreter. ` +
        `Please ensure CCP4 is installed with ccp4-python in ${CCP4Dir}/bin, ` +
        `or create a .venv in ${projectRoot}.`
    );
  }

  if (CCP4I2_PROJECTS_DIR.length > 0) {
    process.env.CCP4I2_PROJECTS_DIR = CCP4I2_PROJECTS_DIR;
    process.env.CCP4I2_DB_FILE = path.join(CCP4I2_PROJECTS_DIR, "db.sqlite3");
  }

  // Set CCP4I2_ROOT to the project root (used by backend for plugin discovery)
  process.env.CCP4I2_ROOT = projectRoot;

  console.log(`🚀 Next.js running on http://localhost:${NEXT_PORT}`);
  console.log(`🐍 Using Python: ${PYTHON_PATH}`);
  console.log(`📁 Project root (CCP4I2_ROOT): ${projectRoot}`);

  const oldCWD = process.cwd();
  const oldPythonPath = process.env.PYTHONPATH || "";
  let serverSrcPath: string;

  if (isDev) {
    serverSrcPath = path.join(process.cwd(), "..", "server");
    process.chdir(path.join(process.cwd(), "..", "server"));
  } else {
    // process.resourcesPath is an Electron-specific property
    const resourcesPath = (process as any).resourcesPath as string;
    serverSrcPath = path.join(resourcesPath, "server");
    process.chdir(path.join(resourcesPath, "server"));
  }

  // Set environment variables BEFORE any Python operations
  // PYTHONPATH needs both server/ (for Django) and projectRoot (for core/ module)
  const pythonPathSeparator = process.platform === "win32" ? ";" : ":";
  const pythonPath = [serverSrcPath, projectRoot].join(pythonPathSeparator);

  const pythonEnv: Record<string, string | undefined> = {
    ...process.env,
    UVICORN_PORT: `${UVICORN_PORT}`,
    NEXT_ADDRESS: `http://localhost:${NEXT_PORT}`,
    // Force local execution mode for Electron app
    EXECUTION_MODE: "local",
    MPLBACKEND: "Agg", // Force matplotlib to use non-GUI backend
    PYTHONPATH: pythonPath,
    CCP4I2_ROOT: projectRoot,
    // Windows-specific DLL path fixes
    ...(process.platform === "win32" && {
      PATH: `${path.join(CCP4Dir, "bin")};${process.env.PATH}`,
      CCP4: CCP4Dir,
      // Force matplotlib to avoid Qt completely
      MPLCONFIGDIR: path.join(os.tmpdir(), "matplotlib-config"),
    }),
  };

  console.log(`🐍`, pythonEnv);
  // Run migrations with the updated environment
  try {
    execSync(`"${PYTHON_PATH}" manage.py migrate`, {
      env: pythonEnv,
      stdio: "inherit", // Show migration output
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
            cwd: serverSrcPath,
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

  console.log(`🐍 serverSrcPath: ${serverSrcPath} ${typeof serverSrcPath}`);

  // Setup logging for production
  let logStream: fs.WriteStream | null = null;
  if (!isDev) {
    const homeDir = os.homedir();
    let logDir = path.join(homeDir, ".ccp4x");
    if (!fs.existsSync(logDir)) {
      logDir = path.join(homeDir, "ccp4x");
    }

    const runNumber = getNextRunNumber(logDir);
    const logFile = path.join(logDir, `uvicorn-${runNumber}.log`);
    fs.mkdirSync(logDir, { recursive: true });
    logStream = fs.createWriteStream(logFile, { flags: "a" });
    console.log(`🐍 Uvicorn logs will be written to: ${logFile}`);
  }

  // Start Python process with the corrected environment
  let pythonProcess: any;
  const uvicornArgs = isDev
    ? ["-m", "uvicorn", "asgi:application", "--reload"]
    : ["-m", "uvicorn", "asgi:application"];

  pythonProcess = spawn(PYTHON_PATH, uvicornArgs, {
    env: pythonEnv,
    shell: true,
    cwd: serverSrcPath, // Ensure we're in the right directory
  });

  console.log(`🚀 Uvicorn running on http://localhost:${UVICORN_PORT}`);

  // Restore original environment
  process.env.PYTHONPATH = oldPythonPath;
  process.chdir(oldCWD);

  pythonProcess.stdout.on("data", (data) => {
    if (isDev) {
      console.log(`🐍 Python Output: ${data}`);
    } else {
      logStream?.write(`[STDOUT] ${new Date().toISOString()}: ${data}`);
    }
  });

  pythonProcess.stderr.on("data", (data) => {
    if (isDev) {
      console.error(`🐍 Python Error: ${data}`);
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
