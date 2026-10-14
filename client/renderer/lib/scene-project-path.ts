/**
 * The `projectPath` scene ref: a file named by its path in the project
 * directory, for files the database does not record — above all what a
 * running job has written so far. It loads through the project's file
 * endpoint, and that loader URL becomes the molecule's or map's uniqueId, so
 * the same URL shape is how a loaded file is recognised again (re-applying a
 * scene, lifting the live view).
 */

export const projectPathUrl = (projectPk: number, path: string): string =>
  `/api/proxy/ccp4i2/projects/${projectPk}/project_file?path=${encodeURIComponent(path)}`;

/** The project-relative path a loader URL fetched, or null for any other URL. */
export function extractProjectPathFromUniqueId(uniqueId: string): string | null {
  const match = uniqueId.match(/\/projects\/\d+\/project_file\/?\?path=([^&#]+)/);
  if (!match) return null;
  try {
    return decodeURIComponent(match[1]);
  } catch {
    return null;
  }
}
