/**
 * Client-side fetcher for the Selection endpoints.
 *
 * A Selection is a server-persisted snapshot of a compound list — created
 * implicitly on every successful NLP query, with a 7-day default expiry.
 * Marking ``is_saved=true`` clears the expiry; the chemist's session list
 * surfaces saved + non-expired ephemeral rows.
 *
 * Mirrors the wire format of ``apps/compounds/registry/selection_views.py``
 * — the ``_serialize`` helper there is the source of truth.
 */
import { apiDelete, apiGet, apiPatch } from './api';

export interface SelectionSummary {
  id: string;
  name: string;
  n_compounds: number;
  created_at: string;
  expires_at: string | null;
  is_saved: boolean;
  target_id: string | null;
  target_name: string | null;
  source_prompt: string;
}

export interface SelectionDetail extends SelectionSummary {
  compound_ids: string[];
}

interface SelectionListResponse {
  selections: SelectionSummary[];
}

// SWR cache key — matches the URL the proxy forwards to Django, so a
// global `mutate(SELECTIONS_LIST_KEY)` from anywhere invalidates the
// session list without prop-drilling.
export const SELECTIONS_LIST_KEY = '/api/proxy/compounds/selections';

export async function listSelections(opts?: {
  isSaved?: boolean;
}): Promise<SelectionSummary[]> {
  const qs =
    opts?.isSaved === undefined ? '' : `?is_saved=${opts.isSaved ? 'true' : 'false'}`;
  const body = await apiGet<SelectionListResponse>(`selections${qs}`);
  return body.selections;
}

export const getSelection = (id: string): Promise<SelectionDetail> =>
  apiGet<SelectionDetail>(`selections/${id}`);

export interface SelectionPatch {
  name?: string;
  is_saved?: boolean;
  target_id?: string | null;
}

export const patchSelection = (
  id: string,
  patch: SelectionPatch,
): Promise<SelectionSummary> =>
  apiPatch<SelectionSummary>(`selections/${id}`, patch);

export const deleteSelection = (id: string): Promise<void> =>
  apiDelete(`selections/${id}`);
