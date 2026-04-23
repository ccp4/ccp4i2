/**
 * Client-side types + fetcher for the NLP compound-selector endpoint.
 *
 * **Pivoted 2026-04-23**: the endpoint now returns a *compound selection*
 * with a pre-built redirect URL to `/assays/aggregate`. The client
 * navigates the user there on confirmation. No table renderer here.
 *
 * Mirrors the dataclasses in apps/compounds/nlp/spec.py — updates must
 * be paired on both sides.
 *
 * Unlike standard compounds endpoints, the NLP endpoint returns
 * meaningful structured bodies on 4xx/5xx (404 disabled, 400 scope/spec
 * error, 429 cap, 502 parse error). The fetcher here returns the JSON
 * body regardless of HTTP status so the UI can render the discriminated
 * union naturally. Throws only on network / JSON-parse failures.
 */
import { authFetch } from './api';

// ---------------------------------------------------------------------------
// Types — mirror apps/compounds/nlp/spec.py exactly
// ---------------------------------------------------------------------------

export interface Threshold {
  op: string;
  value: number;
  unit: string | null;
}

export interface MeasurementFilter {
  protocol_hint: string | null;
  metric: string | null;
  threshold: Threshold | null;
  protocol_id?: string | null;  // set by view on continuation
}

export interface CompoundSelector {
  registration_target_as_typed: string | null;
  assay_target_as_typed: string | null;
  measurement_filters: MeasurementFilter[];
  registration_target_id?: string | null;
  assay_target_id?: string | null;
}

export interface TargetCandidate {
  id: string;
  name: string;
  gene_symbols: string[];
}

export interface ProtocolCandidate {
  id: string;
  name: string;
  n_runs: number;
  last_run: string | null;
  n_compounds: number;
}

export type NLPResponse =
  | {
      status: 'selection';
      redirect_url: string;
      compound_formatted_ids: string[];
      target_names: string[];
      protocol_names: string[];
      n_matched: number;
      scope_sentence: string;
    }
  | {
      status: 'clarify';
      field: string;
      query: string;
      candidates: TargetCandidate[] | ProtocolCandidate[];
      filter_index?: number;           // present for protocol clarify
      partial_selector: CompoundSelector;
    }
  | {
      status: 'miss';
      field: string;
      query: string;
      suggestions: TargetCandidate[] | ProtocolCandidate[];
      filter_index?: number;
    }
  | { status: 'not_a_query'; reason: string }
  | { status: 'error'; kind: string; field?: string; message: string };


// ---------------------------------------------------------------------------
// Field constants — mirror spec.py FIELD_* values
// ---------------------------------------------------------------------------

export const FIELD_REGISTRATION_TARGET = 'registration_target_as_typed';
export const FIELD_ASSAY_TARGET = 'assay_target_as_typed';
export const FIELD_PROTOCOL_HINT = 'protocol_hint';


// ---------------------------------------------------------------------------
// Fetcher — tolerates 4xx/5xx and returns the body
// ---------------------------------------------------------------------------

/**
 * POST /api/proxy/compounds/nlp/query (no trailing slash — Next.js proxy
 * appends when forwarding to Django; see i2remote.py:470 for the same
 * convention).
 */
export async function postNlpQuery(
  body: { prompt: string } | { selector: CompoundSelector },
): Promise<NLPResponse> {
  const res = await authFetch('/api/proxy/compounds/nlp/query', {
    method: 'POST',
    body: JSON.stringify(body),
  });
  let parsed: unknown;
  try {
    parsed = await res.json();
  } catch {
    throw new Error(
      `NLP endpoint returned non-JSON response (HTTP ${res.status})`,
    );
  }
  return parsed as NLPResponse;
}


// ---------------------------------------------------------------------------
// Clarify continuation — apply a picker choice to the partial selector
// ---------------------------------------------------------------------------

/**
 * Pin the user's clarify-picker choice on the appropriate place in the
 * selector and return the next selector to POST as a continuation.
 *
 * - Target clarify: set registration_target_id or assay_target_id on the
 *   selector (depending on `field`).
 * - Protocol clarify: set protocol_id on the specific measurement_filter
 *   identified by `filterIndex`.
 */
export function applyClarifyPick(
  partial: CompoundSelector,
  field: string,
  pickedId: string,
  filterIndex?: number,
): CompoundSelector {
  const next: CompoundSelector = {
    ...partial,
    measurement_filters: partial.measurement_filters.map((f) => ({ ...f })),
  };
  if (field === FIELD_REGISTRATION_TARGET) {
    next.registration_target_id = pickedId;
  } else if (field === FIELD_ASSAY_TARGET) {
    next.assay_target_id = pickedId;
  } else if (field === FIELD_PROTOCOL_HINT) {
    const idx = filterIndex ?? 0;
    if (idx >= 0 && idx < next.measurement_filters.length) {
      next.measurement_filters[idx] = {
        ...next.measurement_filters[idx],
        protocol_id: pickedId,
      };
    }
  }
  return next;
}
