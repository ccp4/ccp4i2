/**
 * The reflection file a new project's "Define AU contents" job should give
 * Matthews (#609). Dropped sequences and dropped data are imported as
 * separate jobs, and the AU contents job was made without the data, so its
 * report had no Matthews analysis to check the copy numbers against. Any
 * MTZ the data import wrote carries the cell and space group Matthews needs;
 * observed data are preferred.
 */
export interface ImportedFile {
  uuid: string;
  type: string;
  job_param_name?: string;
}

export function matthewsReflectionFile(files: ImportedFile[]): string | null {
  const mtz = files.filter((f) => /^application\/CCP4-mtz/.test(f.type ?? ""));
  const chosen =
    mtz.find((f) => f.type === "application/CCP4-mtz-observed") ?? mtz[0];
  return chosen ? chosen.uuid.replace(/-/g, "") : null;
}
