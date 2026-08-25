/**
 * Native path picker.
 *
 * Some task parameters name a path on the machine that runs the job rather
 * than a file to upload — a directory of diffraction images, a BORGES library,
 * an existing run directory. In the desktop build the server is localhost, so
 * Electron's native dialog returns a path the server can open. In the web
 * build there is no such dialog and callers should fall back to letting the
 * user type the path.
 */

export type BrowseMode = "file" | "directory";

export interface BrowsePathOptions {
  mode?: BrowseMode;
  title?: string;
  message?: string;
  buttonLabel?: string;
}

/** True in the Electron desktop build, where the native picker is available. */
export const canBrowsePath = (): boolean =>
  typeof window !== "undefined" && Boolean(window.electronAPI?.invoke);

/**
 * Open the native picker. Resolves to the chosen absolute path, or null if the
 * user cancelled or the build has no native dialog.
 */
export const browsePath = async (
  options: BrowsePathOptions = {}
): Promise<string | null> => {
  const api = typeof window !== "undefined" ? window.electronAPI : undefined;
  if (!api?.invoke) return null;
  try {
    return (await api.invoke("browse-path", options)) ?? null;
  } catch {
    return null;
  }
};
