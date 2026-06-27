/**
 * PDB identifier normalization.
 *
 * The PDB is migrating from the legacy 4-character codes (e.g. `1jst`) to the
 * 12-character extended format `pdb_0000XXXX` (e.g. `pdb_00001jst`), because the
 * 4-character namespace is nearly exhausted. The APIs we call behave as follows
 * (verified against the live endpoints):
 *
 *   - PDBe `api/pdb/entry/files` and `api/pdb/entry/molecules`: accept BOTH
 *     forms, but always KEY the JSON response by the legacy 4-char id even when
 *     queried with the extended id. So a request for `pdb_00001jst` returns
 *     `{ "1jst": { ... } }`. Tolerant key lookup is therefore mandatory.
 *   - RCSB downloads (`files.rcsb.org/download`): accept both forms.
 *   - PDB-REDO (`pdb-redo.eu/db`): accept both forms.
 *   - RCSB FASTA (`rcsb.org/fasta/entry`): accepts ONLY the 4-char form;
 *     the extended form 404s. Use {@link toShortPdbId} for that endpoint.
 *
 * Policy: when an endpoint accepts both, we send the extended (future-proof)
 * form; when an endpoint only accepts the legacy form, we send the short form.
 * Either form may be typed by the user.
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
