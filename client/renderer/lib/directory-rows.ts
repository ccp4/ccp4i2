/**
 * Rows for the lazy directory browser.
 *
 * The browser holds one listing per open folder, keyed by the folder's path
 * relative to the root it browses ("" is the root), and draws the open part
 * of the tree as a flat list of rows so that a virtualiser can draw only the
 * rows on screen. A folder that is open but not yet listed draws one
 * "loading" row in place of its children.
 */

export interface DirectoryEntry {
  path: string;
  name: string;
  type: "directory" | "file";
  size?: number;
  mtime?: number;
  link?: boolean;
}

export type DirectoryRow =
  | { kind: "entry"; entry: DirectoryEntry; relative: string; depth: number; open: boolean }
  | { kind: "loading"; relative: string; depth: number };

const trailingInt = (name: string): number | null => {
  const match = name.match(/_(\d+)$/);
  return match ? parseInt(match[1], 10) : null;
};

/** Folders first; job_2 before job_10; otherwise by name. */
export function compareEntries(a: DirectoryEntry, b: DirectoryEntry): number {
  if (a.type !== b.type) return a.type === "directory" ? -1 : 1;
  const aInt = trailingInt(a.name);
  const bInt = trailingInt(b.name);
  if (
    aInt !== null &&
    bInt !== null &&
    a.name.replace(/_\d+$/, "") === b.name.replace(/_\d+$/, "")
  ) {
    return aInt - bInt;
  }
  return a.name.localeCompare(b.name);
}

export const childPath = (parent: string, name: string) =>
  parent ? `${parent}/${name}` : name;

/**
 * The visible rows, depth first. ``fileFilter`` hides files only: a folder
 * cannot be judged empty of matches until it is opened.
 */
export function directoryRows(
  listings: Record<string, DirectoryEntry[] | undefined>,
  open: Set<string>,
  fileFilter?: (entry: DirectoryEntry) => boolean
): DirectoryRow[] {
  const rows: DirectoryRow[] = [];
  const walk = (folder: string, depth: number) => {
    const entries = listings[folder];
    if (!entries) {
      rows.push({ kind: "loading", relative: folder, depth });
      return;
    }
    const shown = entries
      .filter((e) => e.type === "directory" || !fileFilter || fileFilter(e))
      .sort(compareEntries);
    for (const entry of shown) {
      const relative = childPath(folder, entry.name);
      const isOpen = entry.type === "directory" && open.has(relative);
      rows.push({ kind: "entry", entry, relative, depth, open: isOpen });
      if (isOpen) walk(relative, depth + 1);
    }
  };
  walk("", 0);
  return rows;
}
