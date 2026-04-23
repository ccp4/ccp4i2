'use client';

import { useCallback, useMemo, useState } from 'react';
import {
  Alert,
  Autocomplete,
  Box,
  Button,
  FormControl,
  IconButton,
  InputLabel,
  MenuItem,
  Paper,
  Select,
  Stack,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
  Tooltip,
  Typography,
} from '@mui/material';
import { Add, Delete, DataObject, ViewModule } from '@mui/icons-material';
import type {
  Protocol,
  ScorecardAxis,
  ScorecardAxisKind,
  ScorecardConfig,
  ScorecardLipinskiAxis,
  ScorecardProtocolAxis,
  ScorecardRatioAxis,
  ScorecardWorstOfAxis,
} from '@/types/compounds/models';

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

interface Props {
  value: ScorecardConfig | null | undefined;
  protocols: Protocol[];
  onChange: (next: ScorecardConfig) => void;
  /** Surface structural errors from the server or local validation. */
  error?: string | null;
  disabled?: boolean;
}

/** Editor primitive — controlled. Renders axes as cards with kind-specific
 *  fields, plus an optional raw-JSON escape hatch for expert edits. */
export function ScorecardEditor({ value, protocols, onChange, error, disabled }: Props) {
  const axes = value?.axes ?? [];
  const [mode, setMode] = useState<'form' | 'json'>('form');
  const [jsonDraft, setJsonDraft] = useState<string>(() =>
    JSON.stringify(value ?? { axes: [] }, null, 2),
  );
  const [jsonError, setJsonError] = useState<string | null>(null);

  const setAxes = useCallback(
    (next: ScorecardAxis[]) => {
      onChange({ axes: next });
    },
    [onChange],
  );

  const updateAxis = useCallback(
    (i: number, patch: Partial<ScorecardAxis>) => {
      const next = axes.slice();
      next[i] = { ...next[i], ...patch } as ScorecardAxis;
      setAxes(next);
    },
    [axes, setAxes],
  );

  const removeAxis = useCallback(
    (i: number) => setAxes(axes.filter((_, j) => j !== i)),
    [axes, setAxes],
  );

  const addAxis = useCallback(
    (kind: ScorecardAxisKind) => setAxes([...axes, emptyAxis(kind)]),
    [axes, setAxes],
  );

  const handleModeChange = (_: unknown, next: 'form' | 'json' | null) => {
    if (next === null) return;
    if (next === 'json') {
      setJsonDraft(JSON.stringify(value ?? { axes: [] }, null, 2));
      setJsonError(null);
    }
    setMode(next);
  };

  const applyJsonDraft = () => {
    try {
      const parsed = JSON.parse(jsonDraft);
      if (!parsed || typeof parsed !== 'object' || !Array.isArray(parsed.axes)) {
        throw new Error('Must be an object with an "axes" array.');
      }
      onChange(parsed as ScorecardConfig);
      setJsonError(null);
      setMode('form');
    } catch (e) {
      setJsonError(e instanceof Error ? e.message : 'Invalid JSON');
    }
  };

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
        <Typography variant="h6" sx={{ flex: 1 }}>
          Axes
        </Typography>
        <ToggleButtonGroup
          exclusive
          size="small"
          value={mode}
          onChange={handleModeChange}
          disabled={disabled}
        >
          <ToggleButton value="form">
            <ViewModule fontSize="small" sx={{ mr: 0.5 }} /> Form
          </ToggleButton>
          <ToggleButton value="json">
            <DataObject fontSize="small" sx={{ mr: 0.5 }} /> JSON
          </ToggleButton>
        </ToggleButtonGroup>
      </Box>

      {error && <Alert severity="error">{error}</Alert>}

      {mode === 'form' ? (
        <>
          <Stack spacing={2}>
            {axes.map((axis, i) => (
              <AxisCard
                key={i}
                axis={axis}
                protocols={protocols}
                onChange={(patch) => updateAxis(i, patch)}
                onRemove={() => removeAxis(i)}
                disabled={disabled}
              />
            ))}
            {axes.length === 0 && (
              <Alert severity="info">
                No axes yet. Add one to start building the scorecard.
              </Alert>
            )}
          </Stack>
          <AddAxisMenu onAdd={addAxis} disabled={disabled} />
        </>
      ) : (
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
          <Typography variant="caption" color="text.secondary">
            Raw JSON editor. Useful for copy-paste between targets or bulk edits.
          </Typography>
          {jsonError && <Alert severity="error">{jsonError}</Alert>}
          <TextField
            multiline
            minRows={12}
            maxRows={30}
            value={jsonDraft}
            onChange={(e) => setJsonDraft(e.target.value)}
            slotProps={{ input: { sx: { fontFamily: 'monospace', fontSize: '0.85rem' } } }}
            disabled={disabled}
          />
          <Box sx={{ display: 'flex', justifyContent: 'flex-end' }}>
            <Button variant="outlined" onClick={applyJsonDraft} disabled={disabled}>
              Apply JSON
            </Button>
          </Box>
        </Box>
      )}
    </Box>
  );
}

// ---------------------------------------------------------------------------
// Axis card
// ---------------------------------------------------------------------------

const KIND_LABELS: Record<ScorecardAxisKind, string> = {
  protocol: 'Protocol value',
  ratio: 'Ratio of two protocols',
  worst_of: 'Worst of several protocols',
  lipinski: 'Lipinski compliance (0–4)',
};

const KIND_HELP: Record<ScorecardAxisKind, string> = {
  protocol: 'Axis value is the geomean of the selected protocol.',
  ratio: 'Axis value = numerator.geomean / denominator.geomean. Use for selectivity ratios.',
  worst_of: 'Axis value is the worst (least-favourable) geomean across the listed protocols.',
  lipinski: 'Axis value is the count (0–4) of Lipinski rule-of-5 criteria this compound passes.',
};

function AxisCard({
  axis,
  protocols,
  onChange,
  onRemove,
  disabled,
}: {
  axis: ScorecardAxis;
  protocols: Protocol[];
  onChange: (patch: Partial<ScorecardAxis>) => void;
  onRemove: () => void;
  disabled?: boolean;
}) {
  const handleKindChange = (kind: ScorecardAxisKind) => {
    // Kind switch resets kind-specific fields but keeps label + thresholds.
    onChange(
      emptyAxis(kind, {
        label: axis.label,
        target_value: axis.target_value,
        poor_value: axis.poor_value,
        threshold_scale: axis.threshold_scale,
      }),
    );
  };

  return (
    <Paper variant="outlined" sx={{ p: 2, display: 'flex', flexDirection: 'column', gap: 1.5 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
        <Typography variant="subtitle1" fontWeight={600} sx={{ flex: 1 }}>
          {axis.label || <em style={{ color: '#999' }}>(unnamed axis)</em>}
        </Typography>
        <Tooltip title="Remove axis" arrow>
          <span>
            <IconButton onClick={onRemove} disabled={disabled} size="small">
              <Delete fontSize="small" />
            </IconButton>
          </span>
        </Tooltip>
      </Box>

      <Box sx={{ display: 'flex', gap: 1.5, flexWrap: 'wrap' }}>
        <FormControl size="small" sx={{ minWidth: 220 }} disabled={disabled}>
          <InputLabel>Axis kind</InputLabel>
          <Select
            value={axis.kind}
            label="Axis kind"
            onChange={(e) => handleKindChange(e.target.value as ScorecardAxisKind)}
          >
            {(Object.keys(KIND_LABELS) as ScorecardAxisKind[]).map((k) => (
              <MenuItem key={k} value={k}>
                {KIND_LABELS[k]}
              </MenuItem>
            ))}
          </Select>
        </FormControl>
        <TextField
          size="small"
          label="Label"
          value={axis.label}
          onChange={(e) => onChange({ label: e.target.value })}
          sx={{ flex: 1, minWidth: 220 }}
          disabled={disabled}
        />
      </Box>

      <Typography variant="caption" color="text.secondary">
        {KIND_HELP[axis.kind]}
      </Typography>

      <KindSpecificFields axis={axis} protocols={protocols} onChange={onChange} disabled={disabled} />

      <Box sx={{ display: 'flex', gap: 1.5, flexWrap: 'wrap' }}>
        <TextField
          size="small"
          label="Excellent value"
          type="number"
          value={axis.target_value ?? ''}
          onChange={(e) => onChange({ target_value: parseNumber(e.target.value) })}
          inputProps={{ step: 'any' }}
          sx={{ flex: 1, minWidth: 160 }}
          disabled={disabled}
          helperText="Value considered excellent"
        />
        <TextField
          size="small"
          label="Poor value"
          type="number"
          value={axis.poor_value ?? ''}
          onChange={(e) => onChange({ poor_value: parseNumber(e.target.value) })}
          inputProps={{ step: 'any' }}
          sx={{ flex: 1, minWidth: 160 }}
          disabled={disabled}
          helperText="Value considered poor"
        />
        <FormControl size="small" sx={{ minWidth: 120 }} disabled={disabled}>
          <InputLabel>Scale</InputLabel>
          <Select
            value={axis.threshold_scale ?? (axis.kind === 'lipinski' ? 'linear' : 'log')}
            label="Scale"
            onChange={(e) => onChange({ threshold_scale: e.target.value as 'log' | 'linear' })}
          >
            <MenuItem value="log">Log</MenuItem>
            <MenuItem value="linear">Linear</MenuItem>
          </Select>
        </FormControl>
      </Box>
    </Paper>
  );
}

// ---------------------------------------------------------------------------
// Kind-specific fields
// ---------------------------------------------------------------------------

function KindSpecificFields({
  axis,
  protocols,
  onChange,
  disabled,
}: {
  axis: ScorecardAxis;
  protocols: Protocol[];
  onChange: (patch: Partial<ScorecardAxis>) => void;
  disabled?: boolean;
}) {
  switch (axis.kind) {
    case 'protocol':
      return (
        <ProtocolPicker
          label="Protocol"
          value={axis.protocol_id}
          protocols={protocols}
          onChange={(id) => onChange({ protocol_id: id } as Partial<ScorecardProtocolAxis>)}
          disabled={disabled}
        />
      );
    case 'ratio':
      return (
        <Box sx={{ display: 'flex', gap: 1.5, flexWrap: 'wrap' }}>
          <ProtocolPicker
            label="Numerator"
            value={axis.numerator_id}
            protocols={protocols}
            onChange={(id) =>
              onChange({ numerator_id: id } as Partial<ScorecardRatioAxis>)
            }
            disabled={disabled}
            sx={{ flex: 1, minWidth: 240 }}
          />
          <ProtocolPicker
            label="Denominator"
            value={axis.denominator_id}
            protocols={protocols}
            onChange={(id) =>
              onChange({ denominator_id: id } as Partial<ScorecardRatioAxis>)
            }
            disabled={disabled}
            sx={{ flex: 1, minWidth: 240 }}
          />
        </Box>
      );
    case 'worst_of':
      return (
        <ProtocolMultiPicker
          label="Protocols"
          value={axis.protocol_ids}
          protocols={protocols}
          onChange={(ids) =>
            onChange({ protocol_ids: ids } as Partial<ScorecardWorstOfAxis>)
          }
          disabled={disabled}
        />
      );
    case 'lipinski':
      return null; // No kind-specific fields
  }
}

// ---------------------------------------------------------------------------
// Pickers
// ---------------------------------------------------------------------------

function ProtocolPicker({
  label,
  value,
  protocols,
  onChange,
  disabled,
  sx,
}: {
  label: string;
  value: string | null | undefined;
  protocols: Protocol[];
  onChange: (id: string) => void;
  disabled?: boolean;
  sx?: object;
}) {
  const selected = useMemo(
    () => protocols.find((p) => p.id === value) ?? null,
    [protocols, value],
  );
  return (
    <Autocomplete
      options={protocols}
      value={selected}
      getOptionLabel={(p) => p.name}
      onChange={(_, p) => p && onChange(p.id)}
      disabled={disabled}
      size="small"
      sx={sx ?? { minWidth: 280 }}
      renderInput={(params) => <TextField {...params} label={label} />}
    />
  );
}

function ProtocolMultiPicker({
  label,
  value,
  protocols,
  onChange,
  disabled,
}: {
  label: string;
  value: string[];
  protocols: Protocol[];
  onChange: (ids: string[]) => void;
  disabled?: boolean;
}) {
  const selected = useMemo(
    () => protocols.filter((p) => value.includes(p.id)),
    [protocols, value],
  );
  return (
    <Autocomplete
      multiple
      options={protocols}
      value={selected}
      getOptionLabel={(p) => p.name}
      onChange={(_, picked) => onChange(picked.map((p) => p.id))}
      disabled={disabled}
      size="small"
      renderInput={(params) => <TextField {...params} label={label} />}
    />
  );
}

function AddAxisMenu({
  onAdd,
  disabled,
}: {
  onAdd: (kind: ScorecardAxisKind) => void;
  disabled?: boolean;
}) {
  return (
    <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
      {(Object.keys(KIND_LABELS) as ScorecardAxisKind[]).map((k) => (
        <Button
          key={k}
          variant="outlined"
          startIcon={<Add />}
          onClick={() => onAdd(k)}
          disabled={disabled}
          size="small"
        >
          {KIND_LABELS[k]}
        </Button>
      ))}
    </Box>
  );
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function emptyAxis(
  kind: ScorecardAxisKind,
  base?: {
    label?: string;
    target_value?: number | null;
    poor_value?: number | null;
    threshold_scale?: 'log' | 'linear';
  },
): ScorecardAxis {
  const shared = {
    label: base?.label ?? '',
    target_value: base?.target_value ?? null,
    poor_value: base?.poor_value ?? null,
    threshold_scale: base?.threshold_scale,
  };
  switch (kind) {
    case 'protocol':
      return { ...shared, kind, protocol_id: '' } as ScorecardProtocolAxis;
    case 'ratio':
      return { ...shared, kind, numerator_id: '', denominator_id: '' } as ScorecardRatioAxis;
    case 'worst_of':
      return { ...shared, kind, protocol_ids: [] } as ScorecardWorstOfAxis;
    case 'lipinski':
      return {
        ...shared,
        kind,
        target_value: base?.target_value ?? 4,
        poor_value: base?.poor_value ?? 1,
        threshold_scale: base?.threshold_scale ?? 'linear',
      } as ScorecardLipinskiAxis;
  }
}

function parseNumber(s: string): number | null {
  if (s.trim() === '') return null;
  const n = Number(s);
  return Number.isFinite(n) ? n : null;
}

