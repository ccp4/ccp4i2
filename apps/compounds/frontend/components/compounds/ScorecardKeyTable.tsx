'use client';

import { Fragment, useMemo } from 'react';
import { Box, Typography } from '@mui/material';
import { groupAxesBySector, sectorColour } from '@/lib/compounds/scorecard';
import { formatConcentrationValue } from '@/lib/compounds/aggregation-api';
import type {
  ConcentrationDisplayMode,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import type { ScorecardAxis, ScorecardConfig } from '@/types/compounds/models';

interface Props {
  config: ScorecardConfig;
  /** Used to resolve protocol names referenced by axis IDs and inherit
   *  per-axis kpi_units for protocol-kind / worst_of axes. */
  protocols?: ProtocolInfo[];
  /** Honour the global concentration display mode for thresholds that
   *  carry units (protocol/worst_of axes). */
  concentrationDisplay?: ConcentrationDisplayMode;
}

/**
 * Compact "key" panel for the scorecard, intended for embedding in
 * SAR slides alongside 4–5 spider cards. Sector-grouped, one row per
 * axis, with a plain-language formula and direction-aware thresholds
 * (excellent → poor). Uses the same sector palette and ordering as
 * the spider so there's no second mental mapping.
 */
export function ScorecardKeyTable({
  config,
  protocols = [],
  concentrationDisplay = 'natural',
}: Props) {
  const groups = useMemo(
    () => groupAxesBySector(config.axes),
    [config.axes],
  );

  if (config.axes.length === 0) {
    return (
      <Typography variant="caption" color="text.secondary">
        No axes configured.
      </Typography>
    );
  }

  return (
    // Native <table> rather than CSS Grid: html2canvas has only partial
    // Grid support and was mis-positioning text in copied PNGs (words
    // overlapping). Tables are html2canvas-friendly and render
    // identically in the browser.
    <Box
      component="table"
      sx={{
        width: '100%',
        borderCollapse: 'collapse',
        fontSize: '0.78rem',
        lineHeight: 1.35,
      }}
    >
      <tbody>
        {groups.map((group, gi) => (
          <Fragment key={gi}>
            {/* Sector banner row spanning all columns */}
            <tr>
              <td
                colSpan={3}
                style={{
                  borderBottom: '1px solid #e0e0e0',
                  paddingTop: gi === 0 ? 0 : 8,
                  paddingBottom: 1,
                  whiteSpace: 'nowrap',
                }}
              >
                <span
                  style={{
                    display: 'inline-block',
                    width: 8,
                    height: 8,
                    borderRadius: '50%',
                    background: group.sector ? sectorColour(group.sector) : 'transparent',
                    border: group.sector ? 'none' : '1px dashed rgba(0,0,0,0.3)',
                    marginRight: 6,
                    verticalAlign: 'middle',
                  }}
                />
                <span
                  style={{
                    fontWeight: 700,
                    textTransform: 'uppercase',
                    fontSize: '0.7rem',
                    color: '#757575',
                  }}
                >
                  {group.sector ?? '—'}
                </span>
              </td>
            </tr>
            {group.items.map(({ item: axis }, i) => (
              <AxisRow
                key={i}
                axis={axis}
                protocols={protocols}
                concentrationDisplay={concentrationDisplay}
              />
            ))}
          </Fragment>
        ))}
      </tbody>
    </Box>
  );
}

function AxisRow({
  axis,
  protocols,
  concentrationDisplay,
}: {
  axis: ScorecardAxis;
  protocols: ProtocolInfo[];
  concentrationDisplay: ConcentrationDisplayMode;
}) {
  return (
    <tr>
      <td style={{ paddingLeft: 12, paddingRight: 10, paddingTop: 1, paddingBottom: 1, fontWeight: 600, verticalAlign: 'baseline', whiteSpace: 'nowrap' }}>
        {axis.label || <em style={{ color: '#999' }}>(unnamed)</em>}
      </td>
      <td
        style={{
          paddingRight: 10,
          paddingTop: 1,
          paddingBottom: 1,
          color: '#757575',
          fontStyle: 'italic',
          // Long formulas wrap within the cell rather than truncating.
          wordBreak: 'break-word',
          verticalAlign: 'baseline',
        }}
        title={formulaText(axis, protocols)}
      >
        {formulaText(axis, protocols)}
      </td>
      <td
        style={{
          fontFamily: 'monospace',
          whiteSpace: 'nowrap',
          textAlign: 'right',
          verticalAlign: 'baseline',
          paddingTop: 1,
          paddingBottom: 1,
        }}
      >
        {thresholdText(axis, protocols, concentrationDisplay)}
      </td>
    </tr>
  );
}

// ---------------------------------------------------------------------------
// Plain-language formula and threshold text per axis kind.
// ---------------------------------------------------------------------------

function formulaText(axis: ScorecardAxis, protocols: ProtocolInfo[]): string {
  const nameOf = (id: string | null | undefined) =>
    (id && protocols.find((p) => p.id === id)?.name) || '?';

  switch (axis.kind) {
    case 'protocol':
      return `${nameOf(axis.protocol_id)} (geomean)`;
    case 'ratio':
      return `${nameOf(axis.numerator_id)} ÷ ${nameOf(axis.denominator_id)}`;
    case 'worst_of': {
      const names = (axis.protocol_ids ?? []).map(nameOf);
      if (names.length === 0) return '(no protocols)';
      if (names.length === 1) return names[0];
      return `worst of: ${names.join(', ')}`;
    }
    case 'lipinski':
      return 'rule-of-5 count (MW≤500, cLogP≤5, HBD≤5, HBA≤10)';
  }
}

function thresholdText(
  axis: ScorecardAxis,
  protocols: ProtocolInfo[],
  concentrationDisplay: ConcentrationDisplayMode,
): string {
  const target = axis.target_value;
  const poor = axis.poor_value;
  if (target == null || poor == null) {
    return '— uncoloured —';
  }
  const fmt = (v: number) => formatAxisValue(axis, v, protocols, concentrationDisplay);
  // Direction-aware: target < poor → lower-better → "≤ target → ≥ poor".
  // target > poor → higher-better → "≥ target → ≤ poor".
  if (target < poor) {
    return `★ ≤ ${fmt(target)}   ✗ ≥ ${fmt(poor)}`;
  }
  return `★ ≥ ${fmt(target)}   ✗ ≤ ${fmt(poor)}`;
}

function formatAxisValue(
  axis: ScorecardAxis,
  value: number,
  protocols: ProtocolInfo[],
  concentrationDisplay: ConcentrationDisplayMode,
): string {
  if (!Number.isFinite(value)) return '—';
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
