/**
 * Opening files in Moorhen by their path in the project, for files the
 * database does not know: whatever a job directory holds, including what a
 * running job has written so far (a PanDDA dataset's model and maps, say).
 *
 * A path names a file or a folder; a folder stands for the files directly in
 * it that Moorhen can open. Models load first (the first one centres the
 * view), then maps. A ``.cif`` is read before it is placed: coordinates
 * become a model, a monomer dictionary is attached to the models opened with
 * it, or drawn as monomers when it came alone.
 */
import { fileKind } from "./file-kinds";

export interface PathLoadPlan {
  coordinates: string[];
  maps: string[];
  /** Read each to learn whether it is coordinates or a dictionary. */
  cifs: string[];
}

export function planPathLoad(paths: string[]): PathLoadPlan {
  const plan: PathLoadPlan = { coordinates: [], maps: [], cifs: [] };
  for (const path of [...paths].sort()) {
    const kind = fileKind(path).moorhen;
    if (kind === "coordinates") plan.coordinates.push(path);
    else if (kind === "map") plan.maps.push(path);
    else if (kind === "cif") plan.cifs.push(path);
  }
  return plan;
}

/** Whether a path names a file Moorhen opens; anything else is taken as a folder. */
export const isLoadableFile = (path: string) => fileKind(path).moorhen !== null;

export function classifyCif(text: string): "coordinates" | "dictionary" | null {
  if (/^_atom_site\./m.test(text) || /^\s*loop_\s*\n\s*_atom_site\./m.test(text)) {
    return "coordinates";
  }
  if (/_chem_comp_atom\./.test(text) || /^data_comp_/m.test(text)) return "dictionary";
  return null;
}

export const baseName = (path: string) => path.split("/").pop() || path;

/** Masks are drawn as masks; PanDDA and others say so in the name. */
export const looksLikeMask = (path: string) => /mask/i.test(baseName(path));
