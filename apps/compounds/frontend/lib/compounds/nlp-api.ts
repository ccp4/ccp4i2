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

/**
 * Half-open ISO-date filter range [after, before). Either end nullable.
 * Dates are "YYYY-MM-DD" strings.
 */
export interface DateRange {
  after: string | null;
  before: string | null;
}

export interface MeasurementFilter {
  protocol_hint: string | null;
  metric: string | null;
  threshold: Threshold | null;
  assay_date_range?: DateRange | null;
  assayed_by_as_typed?: string | null;
  protocol_id?: string | null;  // set by view on continuation
  assayed_by_id?: string | null;
}

export interface CompoundSelector {
  registration_target_as_typed: string | null;
  assay_target_as_typed: string | null;
  measurement_filters: MeasurementFilter[];
  registered_date_range?: DateRange | null;
  registered_by_as_typed?: string | null;
  scaffold_hints?: string[];
  /**
   * Pinned compound IDs as the user typed them ("NCL-00026007",
   * "26007", "compound 26007"). UNION semantics with the rest — these
   * appear in the final selection regardless of whether they pass the
   * other predicates. Resolved server-side via
   * compounds.formatting.extract_reg_number.
   */
  compound_refs_as_typed?: string[];
  registration_target_id?: string | null;
  assay_target_id?: string | null;
  registered_by_id?: string | null;
  /**
   * Slice 16: when the registered-by clarify chip pinned a Supplier
   * (rather than a User), the supplier id round-trips here. Mutually
   * exclusive with registered_by_id at any one time.
   */
  registered_by_supplier_id?: string | null;
  scaffold_ids?: string[];
}

/**
 * Assay-selection query (sibling of CompoundSelector). The LLM emits
 * exactly one of {CompoundSelector, AssaySelector}; the view dispatches
 * to the right executor based on type.
 */
export interface AssaySelector {
  target_as_typed: string | null;
  protocol_hint?: string | null;
  date_range?: DateRange | null;
  created_by_as_typed?: string | null;
  scaffold_hints?: string[];
  // Pinnings from clarify continuation
  target_id?: string | null;
  protocol_id?: string | null;
  created_by_id?: string | null;
  scaffold_ids?: string[];
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

export interface UserCandidate {
  id: string;
  display: string;
  email: string | null;
  n_compounds: number;
  n_assays: number;
  /**
   * Discriminator (slice 16). Picker chips under
   * ``registered_by_as_typed`` may also carry SupplierCandidate
   * entries; the frontend dispatches on this field to pick the right
   * pin field on continuation. Older payloads without the field
   * default to "user".
   */
  kind?: 'user';
}

export interface SupplierCandidate {
  id: string;
  name: string;
  initials: string | null;
  n_compounds: number;
  is_user_linked: boolean;
  kind: 'supplier';
}

export interface ScaffoldCandidate {
  id: string;        // canonical name — also used as the pinning ID
  name: string;
  aliases: string[];
  smarts: string;
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
      status: 'assay_selection';
      redirect_url: string;
      assay_ids: string[];
      target_names: string[];
      protocol_names: string[];
      n_matched: number;
      scope_sentence: string;
    }
  | {
      status: 'clarify';
      field: string;
      query: string;
      candidates: TargetCandidate[] | ProtocolCandidate[] | UserCandidate[] | ScaffoldCandidate[] | (UserCandidate | SupplierCandidate)[];
      filter_index?: number;           // present for protocol / assayed-by clarify
      scaffold_index?: number;         // present for scaffold clarify
      partial_selector: CompoundSelector;
    }
  | {
      status: 'miss';
      field: string;
      query: string;
      suggestions: TargetCandidate[] | ProtocolCandidate[] | UserCandidate[] | ScaffoldCandidate[] | (UserCandidate | SupplierCandidate)[];
      filter_index?: number;
      scaffold_index?: number;
      /**
       * Present on a MetricMiss (field === FIELD_METRIC). The KPIs that
       * ARE recorded in the filter's scope, so the response can answer
       * itself ("did you mean pIC50?"). Absent on other miss kinds.
       */
      available_metrics?: string[];
    }
  | { status: 'not_a_query'; reason: string }
  | { status: 'error'; kind: string; field?: string; message: string };


// ---------------------------------------------------------------------------
// Field constants — mirror spec.py FIELD_* values
// ---------------------------------------------------------------------------

export const FIELD_REGISTRATION_TARGET = 'registration_target_as_typed';
export const FIELD_ASSAY_TARGET = 'assay_target_as_typed';
export const FIELD_PROTOCOL_HINT = 'protocol_hint';
export const FIELD_REGISTERED_BY = 'registered_by_as_typed';
export const FIELD_ASSAYED_BY = 'assayed_by_as_typed';
export const FIELD_SCAFFOLD_HINT = 'scaffold_hint';
export const FIELD_COMPOUND_REF = 'compound_ref';
export const FIELD_METRIC = 'metric';


// ---------------------------------------------------------------------------
// Fetcher — tolerates 4xx/5xx and returns the body
// ---------------------------------------------------------------------------

/**
 * POST /api/proxy/compounds/nlp/query (no trailing slash — Next.js proxy
 * appends when forwarding to Django; see i2remote.py:470 for the same
 * convention).
 */
export async function postNlpQuery(
  body: { prompt: string } | { selector: CompoundSelector | AssaySelector },
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
 * - Target clarify: set registration_target_id / assay_target_id AND
 *   replace the ambiguous *_as_typed string with the picked target's
 *   canonical name, so the selector now self-documents which programme
 *   was chosen (otherwise the user's original "MYC" would still appear
 *   in the selector even though Myc-Aur was picked).
 * - Protocol clarify: same treatment on the specific measurement_filter
 *   identified by `filterIndex` — set protocol_id and replace the
 *   ambiguous protocol_hint ("HTRF") with the chosen protocol's name
 *   ("CDK4-HTRF cellular rebind" etc.).
 *
 * Raises an error (rather than silently crashing) if `partial` is missing
 * or lacks `measurement_filters`. That indicates a client/server version
 * mismatch — the v1 server emitted `partial_spec` (no measurement_filters);
 * the v2 client expects `partial_selector` (with measurement_filters).
 * Surface the mismatch loudly so it's diagnosable rather than a deep
 * React crash.
 */
function isAssaySelector(
  s: CompoundSelector | AssaySelector,
): s is AssaySelector {
  // Compound selectors always carry a (possibly empty) measurement_filters
  // array; assay selectors don't. One-field discriminator.
  return !Array.isArray((s as CompoundSelector).measurement_filters);
}

function applyScaffoldPin(
  hints: string[] | undefined,
  ids: string[] | undefined,
  index: number,
  canonical: string,
): { hints: string[]; ids: string[] } {
  // Scaffold_hints / _ids are parallel arrays indexed together; pad to
  // length, write into [index], keep the rest as-is.
  const nextHints = [...(hints ?? [])];
  const nextIds = [...(ids ?? [])];
  while (nextHints.length <= index) nextHints.push('');
  while (nextIds.length <= index) nextIds.push('');
  nextHints[index] = canonical;
  nextIds[index] = canonical;   // pinning id for scaffolds is the canonical name
  return { hints: nextHints, ids: nextIds };
}

export function applyClarifyPick(
  partial: CompoundSelector | AssaySelector | undefined,
  field: string,
  candidate: { id: string; name?: string | null; display?: string | null; kind?: string },
  filterIndex?: number,
  scaffoldIndex?: number,
): CompoundSelector | AssaySelector {
  if (
    !partial ||
    (!Array.isArray((partial as CompoundSelector).measurement_filters) &&
      !('target_as_typed' in partial))
  ) {
    throw new Error(
      "Clarify response is missing partial_selector shape — this usually " +
      "means the server is running a build older than the assay-selector " +
      "pivot. Redeploy the server image and retry.",
    );
  }

  // Targets/protocols carry `name`; users carry `display`; scaffolds carry
  // `name` and id == canonical name. `candidate.id` is the pinning value.
  const canonical = candidate.name ?? candidate.display ?? candidate.id;

  if (isAssaySelector(partial)) {
    const next: AssaySelector = { ...partial };
    if (field === FIELD_ASSAY_TARGET) {
      next.target_id = candidate.id;
      if (canonical) next.target_as_typed = canonical;
    } else if (field === FIELD_PROTOCOL_HINT) {
      next.protocol_id = candidate.id;
      if (canonical) next.protocol_hint = canonical;
    } else if (field === FIELD_ASSAYED_BY) {
      next.created_by_id = candidate.id;
      if (canonical) next.created_by_as_typed = canonical;
    } else if (field === FIELD_SCAFFOLD_HINT) {
      const { hints, ids } = applyScaffoldPin(
        next.scaffold_hints, next.scaffold_ids, scaffoldIndex ?? 0, canonical,
      );
      next.scaffold_hints = hints;
      next.scaffold_ids = ids;
    }
    return next;
  }

  // Compound-selector path.
  const next: CompoundSelector = {
    ...partial,
    measurement_filters: partial.measurement_filters.map((f) => ({ ...f })),
  };
  if (field === FIELD_REGISTRATION_TARGET) {
    next.registration_target_id = candidate.id;
    if (canonical) next.registration_target_as_typed = canonical;
  } else if (field === FIELD_ASSAY_TARGET) {
    next.assay_target_id = candidate.id;
    if (canonical) next.assay_target_as_typed = canonical;
  } else if (field === FIELD_REGISTERED_BY) {
    // Slice 16: pin to whichever pool the candidate came from.
    // SupplierCandidate carries `kind: "supplier"`; UserCandidate
    // either carries `kind: "user"` or no kind field at all (back-
    // compat with pre-slice-16 server payloads).
    if (candidate.kind === 'supplier') {
      next.registered_by_supplier_id = candidate.id;
      next.registered_by_id = null;
    } else {
      next.registered_by_id = candidate.id;
      next.registered_by_supplier_id = null;
    }
    if (canonical) next.registered_by_as_typed = canonical;
  } else if (field === FIELD_PROTOCOL_HINT) {
    const idx = filterIndex ?? 0;
    if (idx >= 0 && idx < next.measurement_filters.length) {
      next.measurement_filters[idx] = {
        ...next.measurement_filters[idx],
        protocol_id: candidate.id,
        ...(canonical ? { protocol_hint: canonical } : {}),
      };
    }
  } else if (field === FIELD_ASSAYED_BY) {
    const idx = filterIndex ?? 0;
    if (idx >= 0 && idx < next.measurement_filters.length) {
      next.measurement_filters[idx] = {
        ...next.measurement_filters[idx],
        assayed_by_id: candidate.id,
        ...(canonical ? { assayed_by_as_typed: canonical } : {}),
      };
    }
  } else if (field === FIELD_SCAFFOLD_HINT) {
    const { hints, ids } = applyScaffoldPin(
      next.scaffold_hints, next.scaffold_ids, scaffoldIndex ?? 0, canonical,
    );
    next.scaffold_hints = hints;
    next.scaffold_ids = ids;
  }
  return next;
}
