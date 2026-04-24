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
  CompoundSelector,
  NLPResponse,
  ProtocolCandidate,
  TargetCandidate,
} from '@/lib/compounds/nlp-api';

const PREVIEW_COMPOUNDS = 5;

interface Props {
  response: NLPResponse;
  onClarifyPick: (
    partial: CompoundSelector,
    field: string,
    candidate: TargetCandidate | ProtocolCandidate,
    filterIndex?: number,
  ) => void;
}

export function NLPResults({ response, onClarifyPick }: Props) {
  switch (response.status) {
    case 'selection':
      return <SelectionView response={response} />;
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
// Clarify — chip picker
// ---------------------------------------------------------------------------

function ClarifyView({
  response,
  onPick,
}: {
  response: Extract<NLPResponse, { status: 'clarify' }>;
  onPick: (
    partial: CompoundSelector,
    field: string,
    candidate: TargetCandidate | ProtocolCandidate,
    filterIndex?: number,
  ) => void;
}) {
  const { field, query, candidates, partial_selector, filter_index } = response;
  const isProtocol = field === 'protocol_hint';
  const heading = isProtocol
    ? `Which ${query !== '' ? `"${query}" ` : ''}protocol did you mean?`
    : `Which target did you mean by "${query}"?`;

  // If the server didn't include partial_selector (e.g. a pre-pivot server
  // still emitting the old partial_spec shape), picking a chip would crash
  // on `.measurement_filters`. Surface the mismatch legibly instead.
  const partialSelectorMissing =
    !partial_selector || !Array.isArray(partial_selector.measurement_filters);

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
            label={<CandidateLabel candidate={candidate} isProtocol={isProtocol} />}
            onClick={
              partialSelectorMissing
                ? undefined
                : () => onPick(partial_selector, field, candidate, filter_index)
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
  isProtocol,
}: {
  candidate: TargetCandidate | ProtocolCandidate;
  isProtocol: boolean;
}) {
  if (isProtocol) {
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
  const { query, suggestions, field } = response;
  const isProtocol = field === 'protocol_hint';
  const kind = isProtocol ? 'protocol' : 'target';

  return (
    <Alert severity="info" sx={{ alignItems: 'flex-start' }}>
      <Typography variant="body2" sx={{ mb: 1 }}>
        No {kind} matched <strong>&ldquo;{query}&rdquo;</strong>.
        {suggestions.length > 0 ? ' Did you mean:' : ''}
      </Typography>
      {suggestions.length > 0 && (
        <Stack direction="row" spacing={1} sx={{ flexWrap: 'wrap', gap: 0.5, mt: 0.5 }}>
          {suggestions.map((s) => (
            <Chip
              key={s.id}
              size="small"
              label={isProtocol ? (s as ProtocolCandidate).name : (s as TargetCandidate).name}
            />
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
