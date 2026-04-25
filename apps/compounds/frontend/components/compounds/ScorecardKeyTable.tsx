'use client';

import { useMemo } from 'react';
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
    <Box
      sx={{
        display: 'flex',
        flexDirection: 'column',
        gap: 0.75,
        // Designed for slide-embedded use: fixed-ish width, compact line
        // height. Bump font sizes only via the parent if you need it
        // larger for a particular paste.
        fontSize: '0.78rem',
        '& .MuiTypography-root': { fontSize: '0.78rem', lineHeight: 1.35 },
      }}
    >
      {groups.map((group, gi) => (
        <Box key={gi}>
          <Box
            sx={{
              display: 'flex',
              alignItems: 'center',
              gap: 0.75,
              borderBottom: '1px solid',
              borderColor: 'divider',
              mb: 0.25,
              pb: 0.1,
            }}
          >
            <Box
              sx={{
                width: 8,
                height: 8,
                borderRadius: '50%',
                bgcolor: group.sector ? sectorColour(group.sector) : 'transparent',
                border: group.sector ? 'none' : '1px dashed rgba(0,0,0,0.3)',
              }}
            />
            <Typography
              sx={{
                fontWeight: 700,
                letterSpacing: '0.06em',
                textTransform: 'uppercase',
                fontSize: '0.7rem',
                color: 'text.secondary',
              }}
            >
              {group.sector ?? '—'}
            </Typography>
          </Box>
          <Box
            sx={{
              display: 'grid',
              gridTemplateColumns: 'minmax(110px, max-content) 1fr auto',
              columnGap: 1.25,
              rowGap: 0.15,
              alignItems: 'baseline',
              pl: 1.5,
            }}
          >
            {group.items.map(({ item: axis }, i) => (
              <AxisRow
                key={i}
                axis={axis}
                protocols={protocols}
                concentrationDisplay={concentrationDisplay}
              />
            ))}
          </Box>
        </Box>
      ))}
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
    <>
      <Typography sx={{ fontWeight: 600 }}>
        {axis.label || <em style={{ color: '#999' }}>(unnamed)</em>}
      </Typography>
      <Typography
        sx={{
          color: 'text.secondary',
          fontStyle: 'italic',
          whiteSpace: 'nowrap',
          overflow: 'hidden',
          textOverflow: 'ellipsis',
        }}
        title={formulaText(axis, protocols)}
      >
        {formulaText(axis, protocols)}
      </Typography>
      <Typography sx={{ fontFamily: 'monospace', whiteSpace: 'nowrap' }}>
        {thresholdText(axis, protocols, concentrationDisplay)}
      </Typography>
    </>
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
