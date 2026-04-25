'use client';

import { useMemo } from 'react';
import { Box, Typography } from '@mui/material';
import {
  evaluateScorecard,
  groupAxesBySector,
  sectorColour,
  type AxisEvaluation,
} from '@/lib/compounds/scorecard';
import { formatConcentrationValue } from '@/lib/compounds/aggregation-api';
import type {
  CompactRow,
  ConcentrationDisplayMode,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import type { ScorecardAxis, ScorecardConfig } from '@/types/compounds/models';

interface Props {
  config: ScorecardConfig;
  compound: CompactRow;
  /** Used to resolve KPI units per axis (per-protocol axes inherit their
   *  protocol's unit; ratios cancel; lipinski is unitless). */
  protocols?: ProtocolInfo[];
  /** Honour the current concentration display mode (natural / nM / uM /
   *  mM / pConc) for protocol-kind / worst_of axes. */
  concentrationDisplay?: ConcentrationDisplayMode;
  /** Optional caption rendered above the table. */
  caption?: string;
  /** Compact mode tightens vertical spacing for use in narrow card cells. */
  dense?: boolean;
}

/**
 * Per-axis values for the chosen compound, in the same sector-grouped
 * order as the spider so a chemist scanning the table maps directly to
 * a wedge. Used both as the live preview in the scorecard editor and
 * as the key+value section in the Cards aggregation view when a
 * scorecard is configured.
 */
export function ScorecardValuesTable({
  config,
  compound,
  protocols = [],
  concentrationDisplay = 'natural',
  caption,
  dense = false,
}: Props) {
  const evals = useMemo(
    () => evaluateScorecard(config, compound),
    [config, compound],
  );
  const orderedEvals = useMemo<AxisEvaluation[]>(() => {
    if (evals.length === 0) return [];
    const groups = groupAxesBySector(evals.map((e) => e.axis));
    return groups.flatMap((g) => g.items).map(({ index }) => evals[index]);
  }, [evals]);

  if (orderedEvals.length === 0) return null;

  return (
    <Box>
      {caption && (
        <Typography
          variant="caption"
          color="text.secondary"
          sx={{ display: 'block', mb: 0.5 }}
        >
          {caption}
        </Typography>
      )}
      {/* Native <table> rather than CSS Grid: html2canvas has only partial
          Grid support and was mis-positioning text in copied PNGs (words
          overlapping). Tables are html2canvas-friendly and render
          identically in the browser. */}
      <Box
        component="table"
        sx={{
          width: '100%',
          borderCollapse: 'collapse',
          fontSize: '0.8rem',
          '& td': {
            verticalAlign: 'middle',
            paddingTop: dense ? '1px' : '2px',
            paddingBottom: dense ? '1px' : '2px',
          },
        }}
      >
        <tbody>
          {orderedEvals.map((ev, i) => (
            <AxisRow
              key={i}
              evaluation={ev}
              protocols={protocols}
              concentrationDisplay={concentrationDisplay}
            />
          ))}
        </tbody>
      </Box>
    </Box>
  );
}

function AxisRow({
  evaluation,
  protocols,
  concentrationDisplay,
}: {
  evaluation: AxisEvaluation;
  protocols: ProtocolInfo[];
  concentrationDisplay: ConcentrationDisplayMode;
}) {
  const { axis, value, t } = evaluation;
  const dotColour = axis.sector ? sectorColour(axis.sector) : 'transparent';

  const display = useMemo(
    () => formatAxisValue(axis, value, protocols, concentrationDisplay),
    [axis, value, protocols, concentrationDisplay],
  );

  return (
    <tr>
      <td style={{ width: 16, paddingRight: 8 }}>
        <span
          style={{
            display: 'inline-block',
            width: 10,
            height: 10,
            borderRadius: '50%',
            background: dotColour,
            border: axis.sector ? 'none' : '1px dashed rgba(0,0,0,0.2)',
          }}
        />
      </td>
      <td style={{ paddingRight: 12 }}>
        {axis.label || <em style={{ color: '#999' }}>(unnamed)</em>}
      </td>
      <td style={{ fontFamily: 'monospace', textAlign: 'right', paddingRight: 12, whiteSpace: 'nowrap' }}>
        {display}
      </td>
      <td style={{ color: tierColour(t), textAlign: 'right', minWidth: 64, whiteSpace: 'nowrap' }}>
        {tierLabel(t)}
      </td>
    </tr>
  );
}

/**
 * Format a raw axis value for display, honouring units where they exist:
 *  - protocol axis → use that protocol's kpi_unit + concentrationDisplay
 *  - ratio axis    → unitless (cancels when both protocols share units)
 *  - worst_of axis → use the first protocol's kpi_unit
 *  - lipinski      → integer count (no unit)
 */
function formatAxisValue(
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
      const unit = firstId
        ? protocols.find((p) => p.id === firstId)?.kpi_unit
        : null;
      return formatWithUnit(value, unit, concentrationDisplay);
    }
    case 'ratio': {
      // Ratio is unitless when numerator/denominator share units. Either
      // way we show it as a bare number with an explicit "×" suffix to
      // emphasise it's a fold-difference, not a concentration.
      return `${formatBare(value)}×`;
    }
    case 'lipinski': {
      // Integer count out of 4 — drop trailing decimals if value is integral.
      return Number.isInteger(value) ? String(value) : value.toFixed(1);
    }
  }
}

function formatWithUnit(
  value: number,
  unit: string | null | undefined,
  concentrationDisplay: ConcentrationDisplayMode,
): string {
  if (!unit) return formatBare(value);
  const { displayValue, displayUnit } = formatConcentrationValue(
    value,
    unit,
    concentrationDisplay,
  );
  return displayUnit ? `${displayValue} ${displayUnit}` : displayValue;
}

function formatBare(v: number): string {
  if (!Number.isFinite(v)) return '—';
  if (Math.abs(v) >= 100 || Math.abs(v) < 0.1) return v.toPrecision(3);
  return v.toFixed(2);
}

function tierLabel(t: number | null): string {
  if (t == null) return 'no data';
  if (t >= 1) return 'excellent';
  if (t >= 2 / 3) return 'good';
  if (t >= 1 / 3) return 'mid';
  if (t > 0) return 'poor';
  return 'failing';
}

function tierColour(t: number | null): string {
  if (t == null) return 'rgba(0, 0, 0, 0.4)';
  // Same hue ramp as the spider/bullets fill — green at top, red at bottom.
  const hue = 5 + (140 - 5) * t;
  return `hsl(${hue.toFixed(0)}, 55%, 38%)`;
}
