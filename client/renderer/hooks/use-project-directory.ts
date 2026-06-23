import { useEffect, useMemo, useState } from "react";
import path from "path";

/**
 * Manages directory selection state and Electron IPC for the project directory
 * field on both the new-project ("create") and edit-project ("move") forms.
 *
 * Create mode: shows a Default/Custom toggle; uses ""
 *   so the parent-directory choice is persisted to the Electron store.
 * Move mode: shows the current directory read-only; uses "select-directory"
 *   for a one-off pick that does NOT alter the global default.
 */
export function useProjectDirectory(name: string, mode: "create" | "move") {
  const [electronAPIAvailable, setElectronAPIAvailable] = useState(false);

  // The parent directory used to compute the target path.
  // Create mode: initialised with a placeholder, then overwritten by get-config.
  // Move mode: empty string = no move requested.
  const [parentDirectory, setParentDirectory] = useState<string>(
    mode === "create" ? "/home/user/CCP4X_PROJECTS" : ""
  );

  // The configured projects root (CCP4I2_PROJECTS_DIR), populated via get-config
  // for both modes so the component can detect whether the project is already
  // at its default location.
  const [projectsDir, setProjectsDir] = useState("");

  // Create mode only: whether the user has switched to "Custom Directory".
  const [customMode, setCustomMode] = useState(false);

  // Move mode only: three-way choice for how to handle the directory.
  // "keep"    – no move; "default" – move to server's default path;
  // "custom"  – move to a user-picked path.
  const [moveMode, setMoveModeState] = useState<"keep" | "default" | "custom">("keep");

  // Map of path → exists, keyed by the checked path so responses from
  // concurrent check-file-exists calls (default dir vs. custom dir) don't
  // overwrite each other.
  // In create mode the initial value for the computed path is true so the
  // submit button is blocked until the first IPC round-trip completes.
  const [existsByPath, setExistsByPath] = useState<Record<string, boolean>>({});

  useEffect(() => {
    if (!window.electronAPI) {
      // Web mode: no file-exists checks possible, leave existsByPath empty.
      return;
    }
    setElectronAPIAvailable(true);

    // Always fetch get-config so we can compute the default directory path
    // in both create and move modes.
    window.electronAPI.sendMessage("get-config");

    const handler = (_event: Electron.IpcRendererEvent, data: any) => {
      if (data.message === "get-config") {
        const dir = data.config?.CCP4I2_PROJECTS_DIR ?? "";
        setProjectsDir(dir);
        // In create mode the parent directory IS the projects root.
        if (mode === "create") setParentDirectory(dir);
      }
      if (data.message === "select-directory" && mode === "move") {
        setParentDirectory(data.path);
      }
      if (data.message === "check-file-exists") {
        setExistsByPath((prev) => ({ ...prev, [data.path]: data.exists }));
      }
    };

    window.electronAPI.onMessage("message-from-main", handler);
    return () => {
      window.electronAPI.removeMessageListener("message-from-main", handler);
    };
  }, [mode]);

  /**
   * What the project directory WOULD be at the default location.
   * Also triggers a check-file-exists so we know if that path is taken.
   */
  const defaultDirectory = useMemo<string>(() => {
    const result = projectsDir ? path.join(projectsDir, name.toLowerCase()) : "";
    if (result && typeof window !== "undefined" && window.electronAPI) {
      window.electronAPI.sendMessage("check-file-exists", { path: result });
    }
    return result;
  }, [projectsDir, name]);

  /**
   * The full computed target directory path.
   * Sends a check-file-exists IPC message each time it changes so
   * directoryExists stays current.
   */
  const computedDirectory = useMemo<string>(() => {
    const result = path.join(parentDirectory || "", name.toLocaleLowerCase());
    if (typeof window !== "undefined" && window.electronAPI) {
      window.electronAPI.sendMessage("check-file-exists", { path: result });
    }
    return result;
  }, [parentDirectory, name]);

  /** Validation error to display on the directory field. */
  const directoryError = useMemo<string>(() => {
    if (mode === "create") {
      if (customMode && computedDirectory.length === 0)
        return "Directory is required";
      // In create mode, block if the target exists (true until IPC confirms otherwise).
      const exists = existsByPath[computedDirectory] ?? true;
      if (computedDirectory.length > 0 && exists)
        return "Directory already exists";
    } else if (moveMode === "default") {
      if (defaultDirectory && existsByPath[defaultDirectory])
        return "Directory already exists";
    } else if (moveMode === "custom") {
      if (!parentDirectory)
        return "Directory is required";
      if (parentDirectory && computedDirectory && existsByPath[computedDirectory])
        return "Directory already exists";
    }
    return "";
  }, [mode, customMode, moveMode, computedDirectory, defaultDirectory, parentDirectory, existsByPath]);

  /**
   * The directory value to pass to the API:
   * - create + default mode     → null           (server picks its own default)
   * - create + custom mode      → computedDirectory
   * - move  + "keep"            → null           (no move)
   * - move  + "custom" + parent → computedDirectory
   * - move  + "custom" + empty  → null           (picker not yet used)
   */
  const effectiveDirectory = useMemo<string | null>(() => {
    if (mode === "create") return customMode ? computedDirectory : null;
    if (moveMode === "custom" && parentDirectory) return computedDirectory;
    return null;
  }, [mode, customMode, moveMode, computedDirectory, parentDirectory]);

  /** Clear parent and reset to "keep" when switching away from "custom". */
  function setMoveMode(v: "keep" | "default" | "custom") {
    setMoveModeState(v);
    if (v !== "custom") setParentDirectory("");
  }

  function handleBrowse() {
    if (!window.electronAPI) return;
    window.electronAPI.sendMessage(
      mode === "create" ? "locate-ccp4i2-project-directory" : "select-directory"
    );
  }

  return {
    // UI props (pass to ProjectDirectoryField)
    electronAPIAvailable,
    parentDirectory,
    computedDirectory,
    defaultDirectory,
    customMode,
    setCustomMode: (v: boolean) => setCustomMode(v),
    moveMode,
    setMoveMode,
    directoryError,
    handleBrowse,
    clearParent: () => setParentDirectory(""),
    // Form submission helpers
    effectiveDirectory,
    directoryHasError: directoryError.length > 0,
  };
}

export type UseProjectDirectoryReturn = ReturnType<typeof useProjectDirectory>;
