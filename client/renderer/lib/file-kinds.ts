/**
 * What can be done with a file found by browsing, judged from its name.
 *
 * Files in a job directory carry no database type, so the extension is all
 * there is. Some extensions are ambiguous (``.cif`` is coordinates or a
 * monomer dictionary); those are settled by whoever opens the file, not here.
 */

/** How the preview dialog shows a file (its ContentSpecification language). */
export type PreviewLanguage =
  | "text"
  | "json"
  | "xml"
  | "cif"
  | "csv"
  | "mtz"
  | "clustalw"
  | "image"
  | "molblock-raw";

/** What Moorhen makes of a file, or of a ``.cif`` that could be either. */
export type MoorhenKind = "coordinates" | "map" | "cif";

export interface FileKind {
  /** In-app preview, or null when the file has no useful one. */
  preview: PreviewLanguage | null;
  /** Label for the preview item: "Preview", "MTZ header", "Table", ... */
  previewLabel: string;
  moorhen: MoorhenKind | null;
  /** Opens in a new tab, served from the project with a read grant. */
  newTab: "html" | "pdf" | null;
  coot: boolean;
  ccp4mg: boolean;
  viewhkl: boolean;
}

export const extensionOf = (name: string): string => {
  const lower = name.toLowerCase();
  // PanDDA and others name maps "<x>.native.ccp4"; the last part decides.
  const dot = lower.lastIndexOf(".");
  return dot > 0 ? lower.slice(dot + 1) : "";
};

const TEXT = new Set([
  "log", "txt", "out", "err", "inp", "com", "sh", "py", "phil", "cards",
  "def", "eff", "params", "yaml", "yml", "dict", "lib", "seq", "pir",
  "fasta", "fa", "ent", "pdb", "mmcif", "md", "tsv", "html", "htm", "dat",
]);
const IMAGES = new Set(["png", "jpg", "jpeg", "gif", "svg", "webp", "bmp"]);
const COORDINATES = new Set(["pdb", "ent", "mmcif"]);
const MAPS = new Set(["ccp4", "map", "mrc"]);
/** Binary formats with no in-app view: offer no text preview of them. */
const BINARY = new Set([
  "ccp4", "map", "mrc", "mtz", "gz", "zip", "tar", "tgz", "bz2", "xz",
  "h5", "nxs", "cbf", "img", "pkl", "npy", "npz", "sqlite", "db", "so",
  "dylib", "ps", "pdf", ...IMAGES,
]);

export function fileKind(name: string): FileKind {
  const ext = extensionOf(name);
  const kind: FileKind = {
    preview: null,
    previewLabel: "Preview",
    moorhen: null,
    newTab: null,
    coot: false,
    ccp4mg: false,
    viewhkl: false,
  };
  if (ext === "json") kind.preview = "json";
  else if (ext === "xml") kind.preview = "xml";
  else if (ext === "cif") kind.preview = "cif";
  else if (ext === "csv") {
    kind.preview = "csv";
    kind.previewLabel = "Table";
  } else if (ext === "aln") kind.preview = "clustalw";
  else if (ext === "mtz") {
    kind.preview = "mtz";
    kind.previewLabel = "MTZ header";
  } else if (IMAGES.has(ext)) {
    kind.preview = "image";
    kind.previewLabel = "View image";
  } else if (ext === "mol" || ext === "sdf") {
    kind.preview = "molblock-raw";
    kind.previewLabel = "View structure";
  } else if (TEXT.has(ext) || !BINARY.has(ext)) {
    // Unknown extensions get a text preview too: the dialog sniffs the bytes
    // and says so when the file turns out to be binary.
    kind.preview = "text";
  }
  if (ext === "html" || ext === "htm") kind.newTab = "html";
  if (ext === "pdf") kind.newTab = "pdf";

  if (COORDINATES.has(ext)) kind.moorhen = "coordinates";
  else if (MAPS.has(ext)) kind.moorhen = "map";
  else if (ext === "cif") kind.moorhen = "cif";

  const structural = COORDINATES.has(ext) || ext === "cif" || ext === "mtz" || MAPS.has(ext);
  kind.coot = structural;
  kind.ccp4mg = structural;
  kind.viewhkl = ext === "mtz";
  return kind;
}

/** 1.2 KB, 3.4 MB, 1.1 GB. */
export function formatSize(size?: number): string {
  if (size === undefined || size === null) return "";
  if (size < 1024) return `${size} B`;
  const units = ["KB", "MB", "GB", "TB"];
  let value = size / 1024;
  let unit = 0;
  while (value >= 1024 && unit < units.length - 1) {
    value /= 1024;
    unit += 1;
  }
  return `${value.toFixed(1)} ${units[unit]}`;
}

/**
 * When a file was last written, short enough for a row: "just now",
 * "4 min ago", "3 h ago", then the date. ``mtime`` is in seconds, as the
 * server sends it.
 */
export function formatModified(mtime?: number, now: number = Date.now()): string {
  if (mtime === undefined || mtime === null) return "";
  const seconds = Math.max(0, now / 1000 - mtime);
  if (seconds < 60) return "just now";
  if (seconds < 3600) return `${Math.floor(seconds / 60)} min ago`;
  if (seconds < 86400) return `${Math.floor(seconds / 3600)} h ago`;
  const date = new Date(mtime * 1000);
  const sameYear = date.getFullYear() === new Date(now).getFullYear();
  return date.toLocaleDateString(undefined, {
    day: "numeric",
    month: "short",
    ...(sameYear ? {} : { year: "numeric" }),
  });
}

/** The full local date and time, for a tooltip or the menu header. */
export const formatTimestamp = (mtime?: number): string =>
  mtime === undefined || mtime === null ? "" : new Date(mtime * 1000).toLocaleString();
