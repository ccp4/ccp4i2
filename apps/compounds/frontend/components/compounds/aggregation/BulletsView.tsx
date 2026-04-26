'use client';

import { Fragment, useMemo, useState } from 'react';
import {
  Box,
  Paper,
  Typography,
  Tooltip,
  Alert,
  FormControl,
  Select,
  MenuItem,
  IconButton,
  ToggleButton,
  ToggleButtonGroup,
} from '@mui/material';
import { ArrowDownward, ArrowUpward } from '@mui/icons-material';
import { useRouter } from 'next/navigation';
import type {
  CompactAggregationResponse,
  CompactRow,
  ConcentrationDisplayMode,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import type { ScorecardAxis, ScorecardConfig } from '@/types/compounds/models';
import {
  evaluateScorecard,
  groupAxesBySector,
  sectorColour,
  type AxisEvaluation,
} from '@/lib/compounds/scorecard';
import { formatConcentrationValue } from '@/lib/compounds/aggregation-api';
import { MoleculeChip } from '../MoleculeView';
import { routes } from '@/lib/compounds/routes';

interface Props {
  data: CompactAggregationResponse;
  scorecardConfig: ScorecardConfig | null | undefined;
  fillHeight?: boolean;
  searchTerm?: string;
  /** Used to honour units in the bullet cell text (per-protocol axes
   *  inherit their protocol's kpi_unit). */
  concentrationDisplay?: ConcentrationDisplayMode;
}

type SortKey = 'compound' | `axis_${number}`;
type Order = 'asc' | 'desc';
type ThumbSize = 'sm' | 'md' | 'lg';

const THUMB_SIZES: Record<ThumbSize, number> = {
  sm: 64,
  md: 110,
  lg: 180,
};

/**
 * Bullet-chart small-multiples view: one row per compound, one cell per
 * scorecard axis. Each cell is a horizontal bar showing the compound's
 * normalised score against that axis, on a poor→excellent gradient.
 *
 * Designed for SAR-glance reading across many compounds at once — the
 * Cards-view spider gets tangled past ~6 compounds, but a bullets table
 * scales linearly. Same scoring logic, different layout.
 */
export function BulletsView({
  data,
  scorecardConfig,
  fillHeight,
  searchTerm = '',
  concentrationDisplay = 'natural',
}: Props) {
  const router = useRouter();
  const [sortKey, setSortKey] = useState<SortKey>('compound');
  // Best-first by default for axis sorts: a high `t` means "near target"
  // regardless of whether the underlying axis is lower-better or higher-
  // better. Compound-ID sorts default to ascending alphabetical.
  const [order, setOrder] = useState<Order>('asc');
  // Chemists vary on how big they want the structure thumbnail. Default
  // to medium (110px) — generous enough to read substituents at a glance
  // but doesn't push compound count off the screen. Local state for now.
  const [thumbSize, setThumbSize] = useState<ThumbSize>('md');
  const thumbPx = THUMB_SIZES[thumbSize];

  const filteredRows = useMemo(() => {
    const all = data.data as CompactRow[];
    if (!searchTerm) return all;
    const term = searchTerm.toLowerCase();
    return all.filter((r) => r.formatted_id?.toLowerCase().includes(term));
  }, [data.data, searchTerm]);

  // Sort rows. For axis sorts we use the normalised score t so 'best
  // first' makes sense regardless of axis direction. Nulls always sink
  // to the bottom regardless of order direction. NB: this useMemo MUST
  // stay above the early-return guard below — declaring it after makes
  // the component's hook count vary with scorecard state (React #310).
  const rows = useMemo(() => {
    const sortValue = (row: CompactRow): number | string | null => {
      if (sortKey === 'compound') return row.formatted_id ?? '';
      const idx = Number(sortKey.slice(5));
      if (!scorecardConfig) return null;
      const evals = evaluateScorecard(scorecardConfig, row);
      return evals[idx]?.t ?? null;
    };
    const sign = order === 'asc' ? 1 : -1;
    return [...filteredRows].sort((a, b) => {
      const va = sortValue(a);
      const vb = sortValue(b);
      if (va == null && vb == null) return 0;
      if (va == null) return 1; // nulls always last
      if (vb == null) return -1;
      if (typeof va === 'number' && typeof vb === 'number') return sign * (va - vb);
      return sign * String(va).localeCompare(String(vb));
    });
  }, [filteredRows, sortKey, order, scorecardConfig]);

  if (!scorecardConfig?.axes?.length) {
    return (
      <Paper sx={{ p: 3 }}>
        <Alert severity="info">
          The bullets view needs a scorecard. Configure one on the target&apos;s
          page (Scorecard button) and pick this format again.
        </Alert>
      </Paper>
    );
  }

  // Reorder axes so same-sector axes are adjacent, and capture the original
  // index so we can map evaluations back to their original axis position.
  const sectorGroups = groupAxesBySector(scorecardConfig.axes);
  const orderedEntries = sectorGroups.flatMap((g) => g.items);
  const axes = orderedEntries.map((entry) => entry.item);
  const orderedIndices = orderedEntries.map((entry) => entry.index);
  const protocols = data.protocols;

  const sortOptions: { value: SortKey; label: string }[] = [
    { value: 'compound', label: 'Compound ID' },
    ...scorecardConfig.axes.map((axis, i) => ({
      value: `axis_${i}` as SortKey,
      label: `★ ${axis.label || `axis ${i + 1}`}`,
    })),
  ];

  return (
    <Paper
      sx={{
        p: 2,
        overflow: 'auto',
        ...(fillHeight && { height: '100%' }),
      }}
    >
      {/* Sort controls */}
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, mb: 1.5 }}>
        <Typography variant="body2" color="text.secondary">
          Sort by:
        </Typography>
        <FormControl size="small" sx={{ minWidth: 200 }}>
          <Select
            value={sortKey}
            onChange={(e) => setSortKey(e.target.value as SortKey)}
            sx={{ fontSize: '0.875rem' }}
          >
            {sortOptions.map((opt) => (
              <MenuItem key={opt.value} value={opt.value}>
                {opt.label}
              </MenuItem>
            ))}
          </Select>
        </FormControl>
        <Tooltip title={order === 'asc' ? 'Ascending' : 'Descending'} arrow>
          <IconButton
            size="small"
            onClick={() => setOrder(order === 'asc' ? 'desc' : 'asc')}
          >
            {order === 'asc' ? <ArrowUpward fontSize="small" /> : <ArrowDownward fontSize="small" />}
          </IconButton>
        </Tooltip>
        <Typography variant="caption" color="text.secondary" sx={{ ml: 1 }}>
          {sortKey.startsWith('axis_')
            ? order === 'desc'
              ? 'best first'
              : 'worst first'
            : ''}
        </Typography>
        <Box sx={{ flex: 1 }} />
        <Typography variant="body2" color="text.secondary">
          Structure:
        </Typography>
        <ToggleButtonGroup
          exclusive
          size="small"
          value={thumbSize}
          onChange={(_, v) => v && setThumbSize(v as ThumbSize)}
        >
          <ToggleButton value="sm" sx={{ px: 1, py: 0.25, fontSize: '0.75rem' }}>S</ToggleButton>
          <ToggleButton value="md" sx={{ px: 1, py: 0.25, fontSize: '0.75rem' }}>M</ToggleButton>
          <ToggleButton value="lg" sx={{ px: 1, py: 0.25, fontSize: '0.75rem' }}>L</ToggleButton>
        </ToggleButtonGroup>
      </Box>

      <Box
        sx={{
          display: 'grid',
          gridTemplateColumns: `${thumbPx}px minmax(140px, max-content) repeat(${axes.length}, minmax(120px, 1fr))`,
          rowGap: 0.5,
          columnGap: 1.5,
          alignItems: 'center',
        }}
      >
        {/* Sector banner row — one tinted band per sector spanning its
            grouped axis columns. Hidden if no axis has a sector. */}
        {sectorGroups.some((g) => g.sector !== null) && (
          <>
            <Box />
            <Box />
            {sectorGroups.map((group, gi) => (
              <Box
                key={`sector-${gi}`}
                sx={{
                  gridColumn: `span ${group.items.length}`,
                  bgcolor: group.sector ? sectorColour(group.sector) : 'transparent',
                  color: group.sector ? '#fff' : 'text.secondary',
                  textAlign: 'center',
                  fontSize: '0.7rem',
                  fontWeight: 700,
                  textTransform: 'uppercase',
                  py: 0.4,
                  borderRadius: 0.5,
                  // Faint tint on the unsectored band so it reads as a band too.
                  border: group.sector ? 'none' : '1px dashed',
                  borderColor: 'divider',
                  whiteSpace: 'nowrap',
                  overflow: 'hidden',
                  textOverflow: 'ellipsis',
                }}
              >
                {group.sector ?? '—'}
              </Box>
            ))}
          </>
        )}

        {/* Axis-label header row */}
        <Box />
        <Box />
        {axes.map((axis, i) => (
          <Tooltip key={i} title={axis.label} arrow>
            <Typography
              variant="caption"
              sx={{
                fontWeight: 600,
                color: 'text.secondary',
                whiteSpace: 'nowrap',
                overflow: 'hidden',
                textOverflow: 'ellipsis',
                pb: 0.5,
                borderBottom: 1,
                borderColor: 'divider',
              }}
            >
              {axis.label}
            </Typography>
          </Tooltip>
        ))}

        {/* Body rows */}
        {rows.map((row) => {
          const evals = evaluateScorecard(scorecardConfig, row);
          return (
            <Fragment key={row.compound_id}>
              <Box
                sx={{ cursor: 'pointer' }}
                onClick={() => router.push(routes.registry.compound(row.compound_id))}
              >
                {row.smiles ? (
                  <MoleculeChip smiles={row.smiles} size={thumbPx} />
                ) : (
                  <Box
                    sx={{
                      width: thumbPx,
                      height: thumbPx,
                      bgcolor: 'grey.100',
                      borderRadius: 0.5,
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                    }}
                  >
                    <Typography variant="caption" color="text.secondary">
                      —
                    </Typography>
                  </Box>
                )}
              </Box>
              <Typography
                variant="body2"
                sx={{
                  fontFamily: 'monospace',
                  fontWeight: 500,
                  cursor: 'pointer',
                  '&:hover': { textDecoration: 'underline' },
                }}
                onClick={() => router.push(routes.registry.compound(row.compound_id))}
              >
                {row.formatted_id}
              </Typography>
              {orderedIndices.map((origIdx, i) => (
                <BulletCell
                  key={i}
                  evaluation={evals[origIdx]}
                  protocols={protocols}
                  concentrationDisplay={concentrationDisplay}
                />
              ))}
            </Fragment>
          );
        })}
      </Box>
    </Paper>
  );
}

// ---------------------------------------------------------------------------

export const GREEN_HUE = 140;
export const RED_HUE = 5;

export function hueAt(t: number): number {
  return RED_HUE + (GREEN_HUE - RED_HUE) * t;
}

/** Background gradient strip (poor → mid → excellent) used as the
 *  "track" behind a bullet's coloured fill. Same in BulletsView and the
 *  compact card body — pulled into a helper so they cannot drift. */
export function bulletTrackGradient(): string {
  return `linear-gradient(to right, hsl(${RED_HUE}, 50%, 92%), hsl(${(RED_HUE + GREEN_HUE) / 2}, 50%, 92%), hsl(${GREEN_HUE}, 50%, 92%))`;
}

/**
 * One bar in the small-multiples grid. Faint background gradient (poor→excellent)
 * with a coloured fill from 0 to t. Tooltip carries the raw value + tier so the
 * cell stays compact while remaining inspectable. Cell text honours the
 * concentration-display mode (natural/nM/uM/mM/pConc) for protocol/worst_of
 * axes; ratios show as ×N (unitless), Lipinski as integer counts.
 */
function BulletCell({
  evaluation,
  protocols,
  concentrationDisplay,
}: {
  evaluation: AxisEvaluation;
  protocols: ProtocolInfo[];
  concentrationDisplay: ConcentrationDisplayMode;
}) {
  const { axis, value, t } = evaluation;

  if (t == null) {
    return (
      <Tooltip title="No data" arrow>
        <Box
          sx={{
            height: 22,
            border: '1px dashed',
            borderColor: 'divider',
            borderRadius: 0.5,
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
          }}
        >
          <Typography variant="caption" color="text.disabled" sx={{ fontSize: '0.65rem' }}>
            —
          </Typography>
        </Box>
      </Tooltip>
    );
  }

  const fillColour = `hsl(${hueAt(t).toFixed(1)}, 65%, 62%)`;
  const bgGradient = bulletTrackGradient();

  const display = formatAxisValueForBullet(axis, value, protocols, concentrationDisplay);
  const tooltip = `${display}  ${tierLabel(t)}`;

  return (
    <Tooltip title={tooltip} arrow>
      <Box
        sx={{
          position: 'relative',
          height: 22,
          background: bgGradient,
          border: 1,
          borderColor: 'divider',
          borderRadius: 0.5,
          overflow: 'hidden',
        }}
      >
        <Box
          sx={{
            position: 'absolute',
            left: 0,
            top: 0,
            bottom: 0,
            width: `${(t * 100).toFixed(1)}%`,
            backgroundColor: fillColour,
          }}
        />
        <Typography
          variant="caption"
          sx={{
            position: 'relative',
            display: 'block',
            textAlign: 'center',
            lineHeight: '22px',
            fontFamily: 'monospace',
            fontSize: '0.7rem',
            fontWeight: 500,
            color: t > 0.5 ? 'rgba(0,0,0,0.85)' : 'text.secondary',
          }}
        >
          {display}
        </Typography>
      </Box>
    </Tooltip>
  );
}

export function formatAxisValueForBullet(
  axis: ScorecardAxis,
  value: number | null,
  protocols: ProtocolInfo[],
  concentrationDisplay: ConcentrationDisplayMode,
): string {
  if (value == null || !Number.isFinite(value)) return '—';
  switch (axis.kind) {
    case 'protocol': {
      const unit = protocols.find((p) => p.id === axis.protocol_id)?.kpi_unit;
      return formatWithUnit(value, unit, concentrationDisplay);
    }
    case 'worst_of': {
      const firstId = axis.protocol_ids?.[0];
      const unit = firstId ? protocols.find((p) => p.id === firstId)?.kpi_unit : null;
      return formatWithUnit(value, unit, concentrationDisplay);
    }
    case 'ratio':
      return `${formatBare(value)}×`;
    case 'lipinski':
      return Number.isInteger(value) ? String(value) : value.toFixed(1);
  }
}

function formatWithUnit(
  value: number,
  unit: string | null | undefined,
  concentrationDisplay: ConcentrationDisplayMode,
): string {
  if (!unit) return formatBare(value);
  const { displayValue, displayUnit } = formatConcentrationValue(value, unit, concentrationDisplay);
  return displayUnit ? `${displayValue} ${displayUnit}` : displayValue;
}

function formatBare(v: number): string {
  if (!Number.isFinite(v)) return '—';
  if (Math.abs(v) >= 100 || Math.abs(v) < 0.1) return v.toPrecision(3);
  return v.toFixed(2);
}

export function tierLabel(t: number): string {
  if (t >= 1) return '(excellent)';
  if (t >= 2 / 3) return '(good)';
  if (t >= 1 / 3) return '(mid)';
  if (t > 0) return '(poor)';
  return '(failing)';
}
