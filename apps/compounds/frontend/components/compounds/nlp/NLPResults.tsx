'use client';

import {
  Alert,
  Box,
  Button,
  Chip,
  Paper,
  Stack,
  Typography,
} from '@mui/material';
import { ArrowForward } from '@mui/icons-material';
import {
  AssaySelector,
  CompoundSelector,
  FIELD_ASSAYED_BY,
  FIELD_COMPOUND_REF,
  FIELD_METRIC,
  FIELD_PROTOCOL_HINT,
  FIELD_REGISTERED_BY,
  FIELD_SCAFFOLD_HINT,
  NLPResponse,
  ProtocolCandidate,
  ScaffoldCandidate,
  TargetCandidate,
  UserCandidate,
} from '@/lib/compounds/nlp-api';

type CandidateKind = 'target' | 'protocol' | 'user' | 'scaffold' | 'compound' | 'metric';

const KIND_FOR_FIELD: Record<string, CandidateKind> = {
  [FIELD_PROTOCOL_HINT]: 'protocol',
  [FIELD_REGISTERED_BY]: 'user',
  [FIELD_ASSAYED_BY]: 'user',
  [FIELD_SCAFFOLD_HINT]: 'scaffold',
  [FIELD_COMPOUND_REF]: 'compound',
  [FIELD_METRIC]: 'metric',
};

function kindForField(field: string): CandidateKind {
  return KIND_FOR_FIELD[field] ?? 'target';
}

const LABEL_FOR_KIND: Record<CandidateKind, string> = {
  protocol: 'protocol',
  user: 'person',
  scaffold: 'substructure',
  compound: 'compound',
  metric: 'metric',
  target: 'target',
};

const PREVIEW_COMPOUNDS = 5;

interface Props {
  response: NLPResponse;
  onClarifyPick: (
    partial: CompoundSelector | AssaySelector,
    field: string,
    candidate: TargetCandidate | ProtocolCandidate | UserCandidate | ScaffoldCandidate,
    filterIndex?: number,
    scaffoldIndex?: number,
  ) => void;
}

export function NLPResults({ response, onClarifyPick }: Props) {
  switch (response.status) {
    case 'selection':
      return <SelectionView response={response} />;
    case 'assay_selection':
      return <AssaySelectionView response={response} />;
    case 'clarify':
      return <ClarifyView response={response} onPick={onClarifyPick} />;
    case 'miss':
      return <MissView response={response} />;
    case 'not_a_query':
      return <NotAQueryView response={response} />;
    case 'error':
      return <ErrorView response={response} />;
  }
}

// ---------------------------------------------------------------------------
// Selection — the happy path. Shows count + preview + "View as cards →" button
// that redirects to /assays/aggregate.
// ---------------------------------------------------------------------------

function SelectionView({
  response,
}: {
  response: Extract<NLPResponse, { status: 'selection' }>;
}) {
  const { n_matched, compound_formatted_ids, scope_sentence, redirect_url } = response;
  const preview = compound_formatted_ids.slice(0, PREVIEW_COMPOUNDS);
  const overflow = Math.max(0, compound_formatted_ids.length - PREVIEW_COMPOUNDS);

  const handleClick = () => {
    // Use window.location so the redirect respects any NEXT_PUBLIC_ROUTE_PREFIX.
    // next/router would re-render the current nested Layout; for a full
    // navigation to the aggregation page a plain location.assign is cleanest.
    window.location.assign(redirect_url);
  };

  return (
    <Paper elevation={1} sx={{ p: 3 }}>
      <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
        {scope_sentence}
      </Typography>

      <Typography variant="h4" sx={{ mb: 2 }}>
        {n_matched === 0 ? (
          'No matching compounds found'
        ) : (
          <>
            Found <strong>{n_matched.toLocaleString()}</strong>{' '}
            matching {n_matched === 1 ? 'compound' : 'compounds'}
          </>
        )}
      </Typography>

      {preview.length > 0 && (
        <Box sx={{ mb: 3 }}>
          <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mb: 0.5 }}>
            First {preview.length}:
          </Typography>
          <Typography variant="body2" sx={{ fontFamily: 'monospace' }}>
            {preview.join(', ')}
            {overflow > 0 && (
              <span style={{ color: 'rgba(0,0,0,0.6)' }}>
                {', and '}
                {overflow.toLocaleString()} more
              </span>
            )}
          </Typography>
        </Box>
      )}

      {n_matched > 0 && (
        <Box sx={{ display: 'flex', justifyContent: 'flex-end' }}>
          <Button
            variant="contained"
            size="large"
            endIcon={<ArrowForward />}
            onClick={handleClick}
          >
            View as cards
          </Button>
        </Box>
      )}
    </Paper>
  );
}

// ---------------------------------------------------------------------------
// Assay selection — like SelectionView but redirects to the /assays list
// rather than /assays/aggregate.
// ---------------------------------------------------------------------------

function AssaySelectionView({
  response,
}: {
  response: Extract<NLPResponse, { status: 'assay_selection' }>;
}) {
  const { n_matched, assay_ids, scope_sentence, redirect_url } = response;
  const previewCount = Math.min(PREVIEW_COMPOUNDS, assay_ids.length);
  const overflow = Math.max(0, assay_ids.length - previewCount);

  const handleClick = () => {
    window.location.assign(redirect_url);
  };

  return (
    <Paper elevation={1} sx={{ p: 3 }}>
      <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
        {scope_sentence}
      </Typography>

      <Typography variant="h4" sx={{ mb: 2 }}>
        {n_matched === 0 ? (
          'No matching assays found'
        ) : (
          <>
            Found <strong>{n_matched.toLocaleString()}</strong>{' '}
            matching {n_matched === 1 ? 'assay' : 'assays'}
          </>
        )}
      </Typography>

      {n_matched > 0 && (
        <Box sx={{ mb: 3 }}>
          <Typography variant="caption" color="text.secondary">
            {overflow > 0
              ? `(${n_matched.toLocaleString()} total — view the list for all of them)`
              : `(${n_matched.toLocaleString()} ${n_matched === 1 ? 'assay' : 'assays'})`}
          </Typography>
        </Box>
      )}

      {n_matched > 0 && (
        <Box sx={{ display: 'flex', justifyContent: 'flex-end' }}>
          <Button
            variant="contained"
            size="large"
            endIcon={<ArrowForward />}
            onClick={handleClick}
          >
            View assay list
          </Button>
        </Box>
      )}
    </Paper>
  );
}

// ---------------------------------------------------------------------------
// Clarify — chip picker
// ---------------------------------------------------------------------------

function ClarifyView({
  response,
  onPick,
}: {
  response: Extract<NLPResponse, { status: 'clarify' }>;
  onPick: (
    partial: CompoundSelector | AssaySelector,
    field: string,
    candidate: TargetCandidate | ProtocolCandidate | UserCandidate | ScaffoldCandidate,
    filterIndex?: number,
    scaffoldIndex?: number,
  ) => void;
}) {
  const { field, query, candidates, partial_selector, filter_index, scaffold_index } = response;
  const kind = kindForField(field);
  const heading =
    kind === 'protocol'
      ? `Which ${query !== '' ? `"${query}" ` : ''}protocol did you mean?`
      : kind === 'user'
        ? `Which person did you mean by "${query}"?`
        : kind === 'scaffold'
          ? `Which substructure did you mean by "${query}"?`
          : `Which target did you mean by "${query}"?`;

  // If the server didn't include partial_selector (e.g. a pre-pivot server
  // still emitting the old partial_spec shape), picking a chip would crash
  // on `.measurement_filters`. Surface the mismatch legibly instead.
  // AssaySelector is valid too — it lacks measurement_filters but carries
  // target_as_typed, so accept either discriminator.
  const partialSelectorMissing =
    !partial_selector ||
    (!Array.isArray((partial_selector as CompoundSelector).measurement_filters) &&
      !('target_as_typed' in partial_selector));

  return (
    <Paper elevation={1} sx={{ p: 3 }}>
      <Typography variant="h6" sx={{ mb: 2 }}>
        {heading}
      </Typography>
      {partialSelectorMissing && (
        <Alert severity="warning" sx={{ mb: 2 }}>
          The server returned a clarify response without the selector context
          needed to continue. This usually means the server is running a
          pre-pivot build. Redeploy the server image and retry the prompt.
        </Alert>
      )}
      <Stack direction="row" spacing={1.5} sx={{ flexWrap: 'wrap', gap: 1.5 }}>
        {candidates.map((candidate) => (
          <Chip
            key={candidate.id}
            label={<CandidateLabel candidate={candidate} kind={kind} />}
            onClick={
              partialSelectorMissing
                ? undefined
                : () => onPick(
                    partial_selector, field, candidate,
                    filter_index, scaffold_index,
                  )
            }
            clickable={!partialSelectorMissing}
            disabled={partialSelectorMissing}
            sx={{ height: 'auto', '& .MuiChip-label': { display: 'block', py: 1, px: 1.5 } }}
          />
        ))}
      </Stack>
    </Paper>
  );
}

function CandidateLabel({
  candidate,
  kind,
}: {
  candidate: TargetCandidate | ProtocolCandidate | UserCandidate | ScaffoldCandidate;
  kind: CandidateKind;
}) {
  if (kind === 'protocol') {
    const p = candidate as ProtocolCandidate;
    const lastRun = p.last_run ? new Date(p.last_run).toLocaleDateString() : 'never';
    return (
      <Box>
        <Typography variant="body2" sx={{ fontWeight: 500 }}>{p.name}</Typography>
        <Typography variant="caption" color="text.secondary">
          {p.n_runs} runs · {p.n_compounds} compounds · last run {lastRun}
        </Typography>
      </Box>
    );
  }
  if (kind === 'user') {
    const u = candidate as UserCandidate;
    const counts: string[] = [];
    if (u.n_compounds > 0) counts.push(`${u.n_compounds} compounds`);
    if (u.n_assays > 0) counts.push(`${u.n_assays} assays`);
    return (
      <Box>
        <Typography variant="body2" sx={{ fontWeight: 500 }}>{u.display}</Typography>
        <Typography variant="caption" color="text.secondary">
          {u.email ?? ''}
          {u.email && counts.length ? ' · ' : ''}
          {counts.join(' · ')}
        </Typography>
      </Box>
    );
  }
  if (kind === 'scaffold') {
    const s = candidate as ScaffoldCandidate;
    return (
      <Box>
        <Typography variant="body2" sx={{ fontWeight: 500 }}>{s.name}</Typography>
        <Typography variant="caption" color="text.secondary" sx={{ fontFamily: 'monospace' }}>
          {s.smarts}
        </Typography>
      </Box>
    );
  }
  const t = candidate as TargetCandidate;
  return (
    <Box>
      <Typography variant="body2" sx={{ fontWeight: 500 }}>{t.name}</Typography>
      {t.gene_symbols.length > 0 && (
        <Typography variant="caption" color="text.secondary">
          {t.gene_symbols.join(', ')}
        </Typography>
      )}
    </Box>
  );
}

// ---------------------------------------------------------------------------
// Miss — suggestions list
// ---------------------------------------------------------------------------

function MissView({
  response,
}: {
  response: Extract<NLPResponse, { status: 'miss' }>;
}) {
  const { query, suggestions, field, available_metrics } = response;
  const kind = kindForField(field);
  const label = LABEL_FOR_KIND[kind];

  const chipLabel = (
    s: TargetCandidate | ProtocolCandidate | UserCandidate | ScaffoldCandidate,
  ): string => {
    if (kind === 'protocol') return (s as ProtocolCandidate).name;
    if (kind === 'user') return (s as UserCandidate).display;
    if (kind === 'scaffold') return (s as ScaffoldCandidate).name;
    return (s as TargetCandidate).name;
  };

  // MetricMiss is special: no fuzzy suggestions, but the response
  // carries `available_metrics` — the KPIs that ARE recorded in scope.
  // Render those as chips so the user can see what to type instead.
  const showAvailableMetrics =
    kind === 'metric' && (available_metrics?.length ?? 0) > 0;

  // Compound-ID miss is special: compound IDs are deterministic, so
  // there are no fuzzy suggestions — instead steer the user to check
  // the format and registration.
  const helperText =
    kind === 'compound'
      ? 'Check the format (e.g. NCL-00026007 or 26007) and that the compound has been registered.'
      : kind === 'metric'
        ? showAvailableMetrics
          ? 'Available KPIs in scope:'
          : 'No measurements in scope under any KPI — try a different protocol or remove the metric.'
        : suggestions.length > 0
          ? 'Did you mean:'
          : '';

  return (
    <Alert severity="info" sx={{ alignItems: 'flex-start' }}>
      <Typography variant="body2" sx={{ mb: 1 }}>
        No {label} matched <strong>&ldquo;{query}&rdquo;</strong>
        {helperText ? `. ${helperText}` : '.'}
      </Typography>
      {showAvailableMetrics && (
        <Stack direction="row" spacing={1} sx={{ flexWrap: 'wrap', gap: 0.5, mt: 0.5 }}>
          {available_metrics!.map((m) => (
            <Chip key={m} size="small" label={m} />
          ))}
        </Stack>
      )}
      {!showAvailableMetrics && suggestions.length > 0 && (
        <Stack direction="row" spacing={1} sx={{ flexWrap: 'wrap', gap: 0.5, mt: 0.5 }}>
          {suggestions.map((s) => (
            <Chip key={s.id} size="small" label={chipLabel(s)} />
          ))}
        </Stack>
      )}
    </Alert>
  );
}

// ---------------------------------------------------------------------------
// Not-a-query
// ---------------------------------------------------------------------------

function NotAQueryView({
  response,
}: {
  response: Extract<NLPResponse, { status: 'not_a_query' }>;
}) {
  return (
    <Alert severity="info">
      <Typography variant="body2" sx={{ mb: 1 }}>
        {response.reason}
      </Typography>
      <Typography variant="caption" color="text.secondary">
        The command bar selects compounds and sends you to the aggregation
        page to view them. It doesn&apos;t generate narratives, register
        compounds, or answer meta-questions. Try phrasing as a selection —
        for example, <em>&ldquo;CDK4 compounds with HTRF IC50 &lt; 100 nM&rdquo;</em>.
      </Typography>
    </Alert>
  );
}

// ---------------------------------------------------------------------------
// Error (server-structured, not a network failure)
// ---------------------------------------------------------------------------

function ErrorView({
  response,
}: {
  response: Extract<NLPResponse, { status: 'error' }>;
}) {
  const severity = response.kind === 'disabled' ? 'warning' : 'error';
  return (
    <Alert severity={severity}>
      <Typography variant="body2" sx={{ fontWeight: 500 }}>
        {response.message}
      </Typography>
      {response.kind !== 'disabled' && (
        <Typography variant="caption" color="text.secondary">
          ({response.kind}
          {response.field ? ` · ${response.field}` : ''})
        </Typography>
      )}
    </Alert>
  );
}
