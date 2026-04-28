/**
 * Substructures catalog fetcher — used by the scatter view's chemotype
 * chip strip.
 */
import { apiGet } from './api';
import type { ScaffoldDefinition } from './scaffold-matching';

interface SubstructuresResponse {
  scaffolds: ScaffoldDefinition[];
}

export const SUBSTRUCTURES_KEY = '/api/proxy/compounds/nlp/substructures';

export async function listSubstructures(opts?: {
  targetId?: string | null;
}): Promise<ScaffoldDefinition[]> {
  const qs = opts?.targetId ? `?target_id=${encodeURIComponent(opts.targetId)}` : '';
  const body = await apiGet<SubstructuresResponse>(`nlp/substructures${qs}`);
  return body.scaffolds;
}
