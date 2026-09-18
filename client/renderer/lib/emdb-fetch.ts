/**
 * Fetching maps from EMDB into a parameter. The server does the work
 * (`repositories/emdb/<entry>/` lists what an entry has, and
 * `jobs/{id}/fetch_repository_file/` fetches one file straight into the
 * project); this module holds the client-side choices, kept pure so they
 * can be tested. See docs/emdb-map-fetch-plan.md.
 */

export interface EmdbFile {
  kind: "map" | "half_map" | "mask";
  file: string;
  sub_type: number;
  label: string;
  index: number;
  pixel_spacing: number | null;
  dimensions: number[] | null;
  size_kbytes: number | null;
  contour_level: number | null;
}

export interface EmdbEntrySummary {
  repository: "emdb";
  entry: string;
  title: string | null;
  resolution: number | null;
  files: EmdbFile[];
  pdb_ids: string[];
}

/** A parameter's requiredSubType qualifier as a list of numbers; the
 *  def.xml value arrives as a number, a string, a comma list or an array. */
export function requiredSubTypes(qualifier: unknown): number[] {
  if (qualifier == null || qualifier === "") return [];
  const raw = Array.isArray(qualifier) ? qualifier : String(qualifier).split(",");
  return raw.map((v) => Number(v)).filter((n) => Number.isFinite(n) && n > 0);
}

/** The file to preselect: the first whose subtype the parameter requires,
 *  else the main map, else the first listed. */
export function pickDefaultEmdbFile(files: EmdbFile[], qualifier: unknown): EmdbFile | null {
  if (files.length === 0) return null;
  const wanted = requiredSubTypes(qualifier);
  const bySubtype = files.find((f) => wanted.includes(f.sub_type));
  return bySubtype ?? files.find((f) => f.kind === "map") ?? files[0];
}

/** For a chosen half map, the other half of the pair, if the entry has it. */
export function otherHalfMap(files: EmdbFile[], chosen: EmdbFile | null): EmdbFile | null {
  if (!chosen || chosen.kind !== "half_map") return null;
  return files.find((f) => f.kind === "half_map" && f.file !== chosen.file) ?? null;
}

/** A one-line description for the chooser. */
export function describeEmdbFile(f: EmdbFile): string {
  const bits: string[] = [];
  if (f.size_kbytes != null) bits.push(`${(f.size_kbytes / 1024).toFixed(0)} MB`);
  if (f.pixel_spacing != null) bits.push(`${f.pixel_spacing.toFixed(2)} Å/px`);
  if (f.dimensions) bits.push(new Set(f.dimensions).size === 1 ? `${f.dimensions[0]}³` : f.dimensions.join("×"));
  return bits.join(", ");
}

/** Among a task's inputs, a sibling map parameter that wants a half map,
 *  to receive the other half of a pair. `lookup` is the container's
 *  objectPath → item map; `objectPath` is the parameter being fetched. */
export function halfMapSibling(
  lookup: Record<string, any> | undefined,
  objectPath: string,
): { objectPath: string; name: string } | null {
  if (!lookup) return null;
  const parent = objectPath.split(".").slice(0, -1).join(".");
  for (const [path, item] of Object.entries(lookup)) {
    if (path === objectPath || !path.startsWith(parent + ".") || path.slice(parent.length + 1).includes(".")) continue;
    if (item?._class !== "CMapDataFile") continue;
    if (!requiredSubTypes(item?._qualifiers?.requiredSubType).includes(5)) continue;
    return { objectPath: path, name: path.split(".").at(-1) as string };
  }
  return null;
}
