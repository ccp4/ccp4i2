/**
 * Serialise a Moorhen molecule for upload to CCP4i2, shared by "push to
 * CCP4i2" (a new coordinate_selector job) and "save to this job" (the
 * recorded Moorhen session), so the two cannot drift.
 */
import { moorhen } from "moorhen/types/moorhen";

export type CoordinateFormat = "pdb" | "mmcif" | "unknown";

export function detectCoordinateFormat(text: string): CoordinateFormat {
  const trimmedText = text.replace(/^\s+/, "");
  if (/^(HEADER|TITLE|ATOM  |HETATM)/m.test(trimmedText)) return "pdb";
  if (/^data_/m.test(trimmedText) && /_atom_site\./.test(trimmedText)) return "mmcif";
  return "unknown";
}

export const slugifyFilename = (name: string) =>
  name
    .replace(/[/\\?%*:|"<>]/g, "")
    .replace(/\s+/g, "_")
    .replace(/[^a-zA-Z0-9._-]/g, "")
    .replace(/^_+|_+$/g, "") || "model";

export interface SerialisedMolecule {
  blob: Blob;
  filename: string;
  format: CoordinateFormat;
  text: string;
}

/** The molecule's current coordinates as a typed Blob with a safe filename,
 *  or null if Moorhen returned nothing. */
export async function serialiseMolecule(
  mol: moorhen.Molecule,
): Promise<SerialisedMolecule | null> {
  const text = await mol.getAtoms();
  if (!text) return null;
  const format = detectCoordinateFormat(text);
  const filename = slugifyFilename(mol.name) + (format === "mmcif" ? ".cif" : ".pdb");
  const blob = new Blob([text], {
    type: format === "mmcif" ? "chemical/x-cif" : "chemical/x-pdb",
  });
  return { blob, filename, format, text };
}
