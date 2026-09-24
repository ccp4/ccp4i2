/**
 * Ligand dictionaries in the Moorhen viewer: which ones belong with a
 * molecule, and how they are attached to it.
 *
 * The rule (decided 2026-09-18): a molecule's dictionaries are the ones its
 * JOB took as input or wrote as output, and nothing else. A project
 * routinely holds several ligands with the same residue name (LIG, DRG), so
 * any wider association is a guess that is wrong as soon as there are two.
 * The server answers the question from the database
 * (`jobs/{id}/dictionaries/`, `files/{id}/companion_dictionaries/`), so the
 * association no longer depends on which code path loaded the molecule.
 *
 * Attaching: a dictionary is read against ONE molecule (Moorhen's
 * `addDict`), never into Coot's global store, because a global entry is
 * inherited by every molecule that has none of its own, including another
 * job's ligand of the same name. Moorhen fetches "missing monomers" from
 * the network while it loads coordinates; that fetch is deferred until the
 * job's own dictionaries are attached, so it only ever asks for what they
 * do not cover.
 */
import { apiGet, apiText } from "../api-fetch";

export const DICTIONARY_TYPE = "application/refmac-dictionary";
export const COORDINATE_TYPES = new Set(["chemical/x-pdb", "chemical/x-cif", "chemical/x-mmcif"]);

export interface DictionaryFile {
  id: number;
  name: string;
  annotation?: string;
  role?: "own" | "input";
}

/** A dictionary ready to attach: its text, and where it came from (so a
 *  captured scene can refer to it by file). */
export interface DictionaryToAttach {
  text: string;
  fileId?: number;
  projectId?: string;
}

/** comp_ids defined by a refmac/coot dictionary CIF (its `data_comp_<X>`
 *  blocks, excluding the `data_comp_list` header). */
export function extractDictCompIds(cifText: string): string[] {
  const out: string[] = [];
  const re = /^data_comp_(\S+)/gm;
  let m: RegExpExecArray | null;
  while ((m = re.exec(cifText)) !== null) {
    if (m[1] !== "list") out.push(m[1]);
  }
  return out;
}

function unwrapList(payload: any): DictionaryFile[] {
  const data = payload && typeof payload === "object" && "success" in payload ? payload.data : payload;
  return Array.isArray(data) ? data : [];
}

/** The dictionaries that belong with a job (its own files and its inputs). */
export async function fetchJobDictionaryFiles(jobId: number): Promise<DictionaryFile[]> {
  try {
    return unwrapList(await apiGet(`jobs/${jobId}/dictionaries/`));
  } catch (err) {
    console.warn(`[dictionaries] could not list dictionaries of job ${jobId}:`, err);
    return [];
  }
}

/** The dictionaries that belong with a file: those of the job it belongs to. */
export async function fetchCompanionDictionaryFiles(fileId: number): Promise<DictionaryFile[]> {
  try {
    return unwrapList(await apiGet(`files/${fileId}/companion_dictionaries/`));
  } catch (err) {
    console.warn(`[dictionaries] could not list companions of file ${fileId}:`, err);
    return [];
  }
}

/** Download each dictionary's text. A file that cannot be read is skipped. */
export async function fetchDictionaryTexts(
  files: DictionaryFile[],
  projectId?: string,
): Promise<DictionaryToAttach[]> {
  const out: DictionaryToAttach[] = [];
  for (const f of files) {
    try {
      const text = await apiText(`/api/proxy/ccp4i2/files/${f.id}/download/`);
      if (text) out.push({ text, fileId: f.id, projectId });
    } catch (err) {
      console.warn(`[dictionaries] could not read dictionary ${f.name}:`, err);
    }
  }
  return out;
}

/** The part of a Moorhen molecule this module relies on. */
export interface DictionaryHost {
  molNo: number | null;
  addDict: (text: string) => Promise<void>;
  loadMissingMonomers?: () => Promise<void>;
}

/**
 * Load coordinates into `mol` with its dictionaries attached before Moorhen
 * goes looking for missing monomers.
 *
 * `load` must perform the coordinate load (Moorhen's loadToCootFromString,
 * which calls `loadMissingMonomers` internally). While it runs that call is
 * a no-op; afterwards the dictionaries are attached to this molecule, the
 * real method is restored, and it runs once for whatever is still missing.
 * With no dictionaries, nothing is deferred and `load` runs untouched.
 */
export async function loadWithDictionaries<T extends DictionaryHost>(
  mol: T,
  dictionaries: DictionaryToAttach[],
  load: () => Promise<unknown>,
): Promise<void> {
  if (dictionaries.length === 0 || typeof mol.loadMissingMonomers !== "function") {
    await load();
    for (const d of dictionaries) await mol.addDict(d.text);
    return;
  }
  const hadOwn = Object.prototype.hasOwnProperty.call(mol, "loadMissingMonomers");
  const original = mol.loadMissingMonomers;
  mol.loadMissingMonomers = async () => {};
  try {
    await load();
    if (mol.molNo != null && mol.molNo !== -1) {
      for (const d of dictionaries) {
        try {
          await mol.addDict(d.text);
        } catch (err) {
          console.warn("[dictionaries] addDict failed:", err);
        }
      }
    }
  } finally {
    if (hadOwn) mol.loadMissingMonomers = original;
    else delete (mol as DictionaryHost).loadMissingMonomers;
  }
  if (mol.molNo != null && mol.molNo !== -1) {
    try {
      await mol.loadMissingMonomers?.();
    } catch (err) {
      console.warn("[dictionaries] deferred missing-monomer fetch failed:", err);
    }
  }
}

/** Provenance to remember for a molecule: comp_id -> the file that defines it. */
export function provenanceOf(
  dictionaries: DictionaryToAttach[],
): Map<string, { fileId: number; projectId?: string }> {
  const out = new Map<string, { fileId: number; projectId?: string }>();
  for (const d of dictionaries) {
    if (d.fileId == null) continue;
    for (const compId of extractDictCompIds(d.text)) {
      out.set(compId, { fileId: d.fileId, projectId: d.projectId });
    }
  }
  return out;
}
