/**
 * Which ligand a campaign dataset is looking for, and which molecule on
 * screen it should be added to. Both answers come from CCP4i2's data model,
 * not from Moorhen (docs/campaign-place-ligand-design.md).
 *
 * The code comes from the job's restraint dictionary, never from the
 * coordinates: at the moment the user asks for the ligand it is not in the
 * coordinates, which is the whole reason they are asking. A refmac LIBOUT
 * may declare standard monomers alongside the fragment, so the dictionary's
 * comp_ids are filtered with the same rule the server's site-scene detector
 * uses (`campaign_scene._is_fragment_code`). The server has gemmi to
 * recognise amino acids, nucleotides and water; this port lists them
 * instead, and copies `COMMON_NON_LIGANDS` verbatim. Keep the two in step.
 * The list is deliberately a constant rather than something fetched: a
 * round trip to learn that GOL is glycerol would be absurd.
 */
import { extractDictCompIds } from "./moorhen-dictionaries";

/** Codes a dictionary may define that are never the fragment being sought. */
export const NON_LIGAND_CODES: ReadonlySet<string> = new Set([
  // Standard amino acids and the variants gemmi tabulates as amino acids.
  "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
  "MSE", "SEC", "PYL", "ASX", "GLX", "UNK",
  // RNA and DNA nucleotides.
  "A", "C", "G", "U", "I", "N", "DA", "DC", "DG", "DT", "DU", "DI", "DN",
  // Water, in every spelling gemmi accepts.
  "HOH", "WAT", "H2O", "DOD", "D2O",
  // campaign_scene.COMMON_NON_LIGANDS: crystallisation additives and ions.
  "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "PG0", "MPD", "BME",
  "SO4", "PO4", "ACT", "ACY", "FMT", "CIT", "FLC", "TLA", "MES", "EPE", "TRS",
  "IPA", "DMS", "NH4", "NO3", "CAC",
  "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CL", "BR", "IOD",
  "CD", "NI", "CO", "CU", "CU1", "HG",
]);

/**
 * The ligand codes a job's dictionaries define, sorted and deduplicated.
 * Several dictionaries or one that defines several fragments both give more
 * than one; the caller must then ask rather than take the first.
 */
export function candidateLigandCodes(dictTexts: string[]): string[] {
  const codes = new Set<string>();
  for (const text of dictTexts) {
    for (const code of extractDictCompIds(text)) {
      const clean = code.trim().toUpperCase();
      if (clean && !NON_LIGAND_CODES.has(clean)) codes.add(clean);
    }
  }
  return [...codes].sort();
}

/** The part of a Moorhen molecule this module relies on. */
export interface LigandHost {
  uniqueId?: string;
  molNo: number | null;
  addLigandOfType: (resType: string, fromMolNo?: number) => Promise<unknown>;
  redraw: () => Promise<void>;
}

/**
 * The molecule loaded from a CCP4i2 file, found by the loader URL that
 * `fetchMolecule` records as `uniqueId`. Never the active or the last
 * molecule: the campaign page can hold extra jobs brought in from the
 * project browser, and merging a fragment into one of those is silent and
 * not obviously wrong on screen.
 */
export function moleculeForFile<T extends { uniqueId?: string }>(
  molecules: T[],
  fileId: number | null | undefined,
): T | undefined {
  if (fileId == null) return undefined;
  const suffix = `/files/${fileId}/download/`;
  return molecules.find((m) => typeof m.uniqueId === "string" && m.uniqueId.endsWith(suffix));
}

/**
 * Add `code` to the molecule loaded from `fileId`, at the view centre.
 * The dictionary is looked up on the same molecule it is merged into,
 * because a dictionary is attached to one molecule and never to Coot's
 * global store (see lib/moorhen-dictionaries). Placed, not fitted.
 */
export async function placeLigand<T extends LigandHost>(
  molecules: T[],
  fileId: number | null | undefined,
  code: string,
): Promise<T> {
  const target = moleculeForFile(molecules, fileId);
  if (!target || target.molNo == null || target.molNo === -1) {
    throw new Error("The dataset's coordinates are no longer loaded");
  }
  await target.addLigandOfType(code, target.molNo);
  await target.redraw();
  return target;
}
