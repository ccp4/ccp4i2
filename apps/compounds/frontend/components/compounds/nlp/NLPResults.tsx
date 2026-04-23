'use client';

import {
  Alert,
  Box,
  Chip,
  Link as MuiLink,
  Paper,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow as MuiTableRow,
  Tooltip,
  Typography,
} from '@mui/material';
import Link from 'next/link';
import { routes } from '@/lib/compounds/routes';
import { PROPERTY_LABELS } from '@/components/compounds/aggregation/shared';
import {
  NLPResponse,
  ProtocolCandidate,
  QuerySpec,
  TableRow as NLPTableRow,
  TargetCandidate,
} from '@/lib/compounds/nlp-api';

// Human-readable label for each FILTER_*/EXCLUDE_* reason — matches the
// constants in apps/compounds/nlp/spec.py.
const EXCLUDE_REASON_LABELS: Record<string, string> = {
  unit_unknown: 'no unit recorded',
  unit_type_mismatch: 'unit-type mismatch',
  unit_incompatible: 'incompatible unit family',
  query_missing_unit: 'threshold needs a unit',
};

interface Props {
  response: NLPResponse;
  onClarifyPick: (
    partialSpec: QuerySpec,
    field: string,
    pickedId: string,
  ) => void;
}

export function NLPResults({ response, onClarifyPick }: Props) {
  switch (response.status) {
    case 'table':
      return <TableView response={response} />;
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
// Table — the happy path
// ---------------------------------------------------------------------------

function TableView({
  response,
}: {
  response: Extract<NLPResponse, { status: 'table' }>;
}) {
  const { rows, property_columns, footer_excluded, scope_sentence } = response;

  const formatValue = (v: number, unit: string | null) => {
    if (!Number.isFinite(v)) return '–';
    const abs = Math.abs(v);
    const formatted =
      abs >= 100 ? v.toFixed(1) : abs >= 1 ? v.toFixed(2) : v.toPrecision(3);
    return unit ? `${formatted} ${unit}` : formatted;
  };

  const formatProperty = (v: number | null | undefined) => {
    if (v === null || v === undefined || !Number.isFinite(v)) return '–';
    return Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(2);
  };

  return (
    <Stack spacing={2}>
      <Typography variant="body2" color="text.secondary">
        {scope_sentence} — {rows.length} {rows.length === 1 ? 'row' : 'rows'}
      </Typography>

      <TableContainer component={Paper} sx={{ maxHeight: '70vh' }}>
        <Table stickyHeader size="small">
          <TableHead>
            <MuiTableRow>
              <TableCell sx={{ fontWeight: 600 }}>Compound</TableCell>
              <TableCell sx={{ fontWeight: 600 }}>Value</TableCell>
              {property_columns.map((prop) => (
                <TableCell
                  key={prop}
                  align="right"
                  sx={{ fontWeight: 600 }}
                >
                  <Tooltip title={prop}>
                    <span>{PROPERTY_LABELS[prop as keyof typeof PROPERTY_LABELS] ?? prop}</span>
                  </Tooltip>
                </TableCell>
              ))}
            </MuiTableRow>
          </TableHead>
          <TableBody>
            {rows.length === 0 ? (
              <MuiTableRow>
                <TableCell colSpan={2 + property_columns.length} align="center">
                  <Typography variant="body2" color="text.secondary" sx={{ py: 3 }}>
                    No rows matched.
                  </Typography>
                </TableCell>
              </MuiTableRow>
            ) : (
              rows.map((row) => <NLPRowRender key={row.compound_id + ':' + row.value} row={row} formatValue={formatValue} formatProperty={formatProperty} propertyColumns={property_columns} />)
            )}
          </TableBody>
        </Table>
      </TableContainer>

      {Object.keys(footer_excluded).length > 0 && (
        <Paper variant="outlined" sx={{ p: 2 }}>
          <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mb: 1 }}>
            Excluded rows (not shown in the table above):
          </Typography>
          <Stack direction="row" spacing={1} sx={{ flexWrap: 'wrap', gap: 0.5 }}>
            {Object.entries(footer_excluded).map(([reason, count]) => (
              <Chip
                key={reason}
                size="small"
                label={`${count} × ${EXCLUDE_REASON_LABELS[reason] ?? reason}`}
              />
            ))}
          </Stack>
        </Paper>
      )}
    </Stack>
  );
}

function NLPRowRender({
  row,
  formatValue,
  formatProperty,
  propertyColumns,
}: {
  row: NLPTableRow;
  formatValue: (v: number, unit: string | null) => string;
  formatProperty: (v: number | null | undefined) => string;
  propertyColumns: string[];
}) {
  return (
    <MuiTableRow hover>
      <TableCell>
        <MuiLink
          component={Link}
          href={routes.registry.compound(row.compound_id)}
          underline="hover"
        >
          {row.formatted_id}
        </MuiLink>
      </TableCell>
      <TableCell>{formatValue(row.value, row.value_unit)}</TableCell>
      {propertyColumns.map((prop) => (
        <TableCell key={prop} align="right">
          {formatProperty(row.properties[prop])}
        </TableCell>
      ))}
    </MuiTableRow>
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
  onPick: (partialSpec: QuerySpec, field: string, pickedId: string) => void;
}) {
  const { field, query, candidates, partial_spec } = response;
  const isProtocol = field === 'protocol_hint';
  const heading = isProtocol
    ? `Which ${query !== '' ? `"${query}" ` : ''}protocol did you mean?`
    : `Which target did you mean by "${query}"?`;

  return (
    <Paper elevation={1} sx={{ p: 3 }}>
      <Typography variant="h6" sx={{ mb: 2 }}>
        {heading}
      </Typography>
      <Stack direction="row" spacing={1.5} sx={{ flexWrap: 'wrap', gap: 1.5 }}>
        {candidates.map((candidate) => (
          <Chip
            key={candidate.id}
            label={<CandidateLabel candidate={candidate} isProtocol={isProtocol} />}
            onClick={() => onPick(partial_spec, field, candidate.id)}
            clickable
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
        The NLP query bar returns tables of compounds and assay results. It
        doesn&apos;t generate narratives, register compounds, or answer
        meta-questions. Try rephrasing as a tabular query, for example:
        &ldquo;Show HTRF IC50 values for CDK4 compounds&rdquo;.
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
