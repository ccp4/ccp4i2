/**
 * Client-side types + fetcher for the NLP query endpoint.
 *
 * Mirrors the dataclasses in apps/compounds/nlp/spec.py. Updates must be
 * paired — if the backend adds a response field, this type gets it too.
 *
 * Unlike the other compounds endpoints, the NLP endpoint returns
 * meaningful structured bodies on 4xx/5xx (e.g. 404 for feature-disabled,
 * 400 for spec errors, 502 for LLM parse failures). The fetcher here
 * returns the JSON body regardless of status so the UI can render the
 * discriminated union naturally. Use `errorOnTransport` for genuine
 * network/JSON failures only.
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

export interface QuerySpec {
  registration_target_as_typed: string | null;
  assay_target_as_typed: string | null;
  protocol_hint: string | null;
  metric: string | null;
  threshold: Threshold | null;
  columns: string | null;
  columns_explicit: string[] | null;
  registration_target_id: string | null;
  assay_target_id: string | null;
  protocol_id: string | null;
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

export interface TableRow {
  compound_id: string;
  formatted_id: string;
  smiles: string | null;
  value: number;
  value_unit: string | null;
  value_in_query_unit: number;
  properties: Record<string, number | null>;
}

export type NLPResponse =
  | {
      status: 'table';
      rows: TableRow[];
      property_columns: string[];
      footer_excluded: Record<string, number>;
      filtered_silent: Record<string, number>;
      scope_sentence: string;
    }
  | {
      status: 'clarify';
      field: string;
      query: string;
      candidates: TargetCandidate[] | ProtocolCandidate[];
      partial_spec: QuerySpec;
    }
  | {
      status: 'miss';
      field: string;
      query: string;
      suggestions: TargetCandidate[] | ProtocolCandidate[];
    }
  | { status: 'not_a_query'; reason: string }
  | { status: 'error'; kind: string; field?: string; message: string };


// ---------------------------------------------------------------------------
// Field constants — mirror the FIELD_* in spec.py. Kept so the UI can tell
// which picker to show (registration target vs assay target vs protocol)
// without string-matching on raw values.
// ---------------------------------------------------------------------------

export const FIELD_REGISTRATION_TARGET = 'registration_target_as_typed';
export const FIELD_ASSAY_TARGET = 'assay_target_as_typed';
export const FIELD_PROTOCOL_HINT = 'protocol_hint';


// ---------------------------------------------------------------------------
// Fetcher — tolerates 4xx/5xx and returns the body.
// ---------------------------------------------------------------------------

/**
 * POST /api/proxy/compounds/nlp/query. Returns the parsed JSON body on
 * any HTTP status (including 404, 400, 429, 502) so the caller can
 * render the discriminated union. Throws only on network / parse failures.
 *
 * Note on URL shape: the Next.js proxy appends `/` when forwarding to
 * Django. Clients MUST NOT include a trailing slash here (see
 * Docker/cli/i2remote.py:470 for the established convention).
 */
export async function postNlpQuery(
  body: { prompt: string } | { spec: QuerySpec },
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
  // We trust the backend to produce a conforming shape — the discriminator
  // parsing happens at the rendering layer.
  return parsed as NLPResponse;
}


// ---------------------------------------------------------------------------
// Convenience helpers
// ---------------------------------------------------------------------------

/**
 * Apply a clarify picker choice — set the pinning ID on the partial_spec
 * for the field the server is asking about, clear any competing pinning
 * (so "re-pick a different target" works).
 */
export function applyClarifyPick(
  partial_spec: QuerySpec,
  field: string,
  pickedId: string,
): QuerySpec {
  const next = { ...partial_spec };
  if (field === FIELD_REGISTRATION_TARGET) {
    next.registration_target_id = pickedId;
  } else if (field === FIELD_ASSAY_TARGET) {
    next.assay_target_id = pickedId;
  } else if (field === FIELD_PROTOCOL_HINT) {
    next.protocol_id = pickedId;
  }
  return next;
}
