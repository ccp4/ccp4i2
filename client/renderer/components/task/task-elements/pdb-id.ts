/**
 * PDB identifier normalization.
 *
 * The PDB is migrating from the legacy 4-character codes (e.g. `1jst`) to the
 * 12-character extended format `pdb_0000XXXX` (e.g. `pdb_00001jst`), because the
 * 4-character namespace is nearly exhausted. The hosts we call behave as follows
 * (verified against the live endpoints):
 *
 *   - PDBe `api/pdb/entry/files` and `api/pdb/entry/molecules` (metadata JSON):
 *     accept BOTH forms, but always KEY the response by the legacy 4-char id
 *     even when queried with the extended id (`pdb_00001jst` -> `{ "1jst": … }`).
 *     The download URLs they echo back mirror the queried form.
 *   - PDBe download host (`entry-files/download`): serves ONLY the legacy 4-char
 *     filename. `pdb_00001jst.cif` 404s; `1jst.cif` works.
 *   - PDB-REDO (`pdb-redo.eu/db`): serves ONLY the legacy 4-char form (extended 500s).
 *   - RCSB FASTA (`rcsb.org/fasta/entry`): accepts ONLY the legacy 4-char form
 *     (extended 404s).
 *   - RCSB downloads (`files.rcsb.org/download`): accept both forms.
 *
 * So the legacy 4-char form is the universally-working one for entries that have
 * one; the extended form is needed only for genuinely-new entries beyond the
 * 4-char namespace (which have no legacy form). Policy: PREFER the legacy form
 * ({@link toFetchPdbId}), falling back to extended only when no legacy form
 * exists. Tolerant key lookup ({@link pdbEntryPayload}) absorbs PDBe keying by
 * the legacy id. Either form may be typed by the user.
 */

/** A legacy 4-char PDB code: a digit followed by three alphanumerics. */
const SHORT_RE = /^[0-9][a-z0-9]{3}$/;

/** An extended id whose 8-char body encodes a legacy 4-char code (`pdb_0000XXXX`). */
const EXTENDED_LEGACY_RE = /^pdb_0000([0-9][a-z0-9]{3})$/;

/**
 * Canonical extended PDB identifier, lowercased.
 * `1jst` -> `pdb_00001jst`; already-extended ids pass through (lowercased);
 * anything unrecognised is returned trimmed/lowercased unchanged.
 */
export function toExtendedPdbId(raw: string): string {
  const id = (raw ?? "").trim().toLowerCase();
  if (id.startsWith("pdb_")) return id;
  if (SHORT_RE.test(id)) return `pdb_0000${id}`;
  return id;
}

/**
 * Legacy 4-char id when one is derivable, else null.
 * `pdb_00001jst` -> `1jst`; `1jst` -> `1jst`; a genuinely-extended id beyond the
 * 4-char namespace (body not starting `0000`) has no short form -> null.
 */
export function toShortPdbId(raw: string): string | null {
  const id = (raw ?? "").trim().toLowerCase();
  if (SHORT_RE.test(id)) return id;
  const m = id.match(EXTENDED_LEGACY_RE);
  return m ? m[1] : null;
}

/**
 * Id to use when issuing a fetch against any PDB host (API query or file
 * download). Prefers the legacy 4-char form — the only form the download/file
 * hosts serve for legacy entries — and falls back to the extended form for
 * genuinely-new entries that have no legacy form. For PDBe metadata queries the
 * queried form also dictates the form of the download URLs echoed back, so this
 * keeps those URLs resolvable.
 */
export function toFetchPdbId(raw: string): string {
  return toShortPdbId(raw) ?? toExtendedPdbId(raw);
}

/**
 * Id best suited for display/filenames: the compact 4-char form when derivable,
 * otherwise the extended form.
 */
export function toDisplayPdbId(raw: string): string {
  return toShortPdbId(raw) ?? toExtendedPdbId(raw);
}

/**
 * Resolve the single entry payload from a PDBe-style response keyed by id.
 * Tries every known id variant, then — since these endpoints return exactly one
 * entry — falls back to the sole value. Returns undefined if nothing matches.
 */
export function pdbEntryPayload<T = any>(
  data: Record<string, T> | null | undefined,
  raw: string,
): T | undefined {
  if (!data || typeof data !== "object") return undefined;
  const candidates = [
    toExtendedPdbId(raw),
    toShortPdbId(raw),
    (raw ?? "").trim().toLowerCase(),
  ].filter(Boolean) as string[];
  for (const key of candidates) {
    if (Object.prototype.hasOwnProperty.call(data, key)) return data[key];
  }
  const values = Object.values(data);
  return values.length === 1 ? values[0] : undefined;
}
