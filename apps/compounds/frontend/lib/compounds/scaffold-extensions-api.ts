/**
 * Fetcher + types for the ScaffoldExtension management surface.
 *
 * The management page uses the slice-24 substructures listing endpoint
 * (which already returns the combined seed + extension catalog), enriched
 * by slice 25 with `canonical_smarts` and a `created_by` username so the
 * chemist can see exactly what their stored pattern matches with after
 * RDKit aromatic perception.
 */
import { apiDelete, apiGet } from './api';

export interface ScaffoldCatalogEntry {
  /** UUID for ScaffoldExtension rows; null for the curated seed catalog. */
  id: string | null;
  name: string;
  aliases: string[];
  /** Raw SMARTS as the chemist typed it. May be Kekulé SMILES. */
  smarts: string;
  /**
   * Canonical aromatic SMARTS produced server-side via the
   * SMILES-first round trip. Null when neither SMILES nor SMARTS
   * parsing produces something sensible — a likely sign the row is
   * unusable as-stored.
   */
  canonical_smarts: string | null;
  /** "seed" | "extension:project" | "extension:shared" */
  source: string;
  /** Free-text gloss (extension only). */
  notes?: string;
  /** Target id when the row is project-scoped, else null. */
  target_id?: string | null;
  target_name?: string | null;
  /** Username of the creator (extension only). */
  created_by?: string | null;
  created_at?: string | null;
}

interface SubstructuresResponse {
  scaffolds: ScaffoldCatalogEntry[];
}

export const SCAFFOLD_CATALOG_KEY = '/api/proxy/compounds/nlp/substructures';

export async function listScaffoldCatalog(opts?: {
  targetId?: string | null;
}): Promise<ScaffoldCatalogEntry[]> {
  const qs = opts?.targetId ? `?target_id=${encodeURIComponent(opts.targetId)}` : '';
  const body = await apiGet<SubstructuresResponse>(`nlp/substructures${qs}`);
  return body.scaffolds;
}

export const deleteScaffoldExtension = (id: string): Promise<void> =>
  apiDelete(`nlp/scaffold/extensions/${id}`);
