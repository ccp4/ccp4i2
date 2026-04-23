/**
 * Evaluate a target's scorecard config against a single compound's aggregation
 * row, producing per-axis raw values and 0-1 normalised scores suitable for
 * feeding to a radar chart.
 *
 * Pure function — no React, no side effects. Mirrors the two-anchor
 * (target_value, poor_value) normalisation used by the protocol colour
 * helper. Direction-agnostic: the sign of (target - poor) flips the sense
 * of "better" naturally through the arithmetic.
 */

import type { CompactRow } from '@/types/compounds/aggregation';
import type {
  ScorecardAxis,
  ScorecardConfig,
  ScorecardProtocolAxis,
  ScorecardRatioAxis,
  ScorecardWorstOfAxis,
} from '@/types/compounds/models';

export interface AxisEvaluation {
  axis: ScorecardAxis;
  /** Raw axis value in native units. null when inputs are missing. */
  value: number | null;
  /** Normalised position: 1 = at or beyond target, 0 = at or beyond poor,
   *  clamped. null when either threshold is unconfigured or value is missing. */
  t: number | null;
}

export function evaluateScorecard(
  config: ScorecardConfig | null | undefined,
  compound: CompactRow,
): AxisEvaluation[] {
  if (!config?.axes) return [];
  return config.axes.map((axis) => evaluateAxis(axis, compound));
}

function evaluateAxis(axis: ScorecardAxis, compound: CompactRow): AxisEvaluation {
  const value = computeAxisValue(axis, compound);
  const t = value == null ? null : normaliseAxis(axis, value);
  return { axis, value, t };
}

function computeAxisValue(axis: ScorecardAxis, compound: CompactRow): number | null {
  switch (axis.kind) {
    case 'protocol':
      return compound.protocols[(axis as ScorecardProtocolAxis).protocol_id]?.geomean ?? null;

    case 'ratio': {
      const { numerator_id, denominator_id } = axis as ScorecardRatioAxis;
      const num = compound.protocols[numerator_id]?.geomean;
      const den = compound.protocols[denominator_id]?.geomean;
      if (num == null || den == null || den === 0) return null;
      return num / den;
    }

    case 'worst_of': {
      const { protocol_ids } = axis as ScorecardWorstOfAxis;
      const values: number[] = [];
      for (const id of protocol_ids) {
        const g = compound.protocols[id]?.geomean;
        if (g != null && Number.isFinite(g)) values.push(g);
      }
      if (values.length === 0) return null;
      // Direction from thresholds: target < poor → lower-better, worst = max.
      const { target_value: target, poor_value: poor } = axis;
      if (target == null || poor == null) return null;
      return target < poor ? Math.max(...values) : Math.min(...values);
    }

    case 'lipinski': {
      const props = compound.properties;
      if (!props) return null;
      let pass = 0;
      if (props.molecular_weight != null && props.molecular_weight <= 500) pass++;
      if (props.clogp != null && props.clogp <= 5) pass++;
      if (props.hbd != null && props.hbd <= 5) pass++;
      if (props.hba != null && props.hba <= 10) pass++;
      return pass;
    }
  }
}

function normaliseAxis(axis: ScorecardAxis, value: number): number | null {
  const { target_value: target, poor_value: poor } = axis;
  if (target == null || poor == null || target === poor) return null;

  const scale = axis.threshold_scale ?? (axis.kind === 'lipinski' ? 'linear' : 'log');
  const canLog = scale === 'log' && value > 0 && target > 0 && poor > 0;

  const v = canLog ? Math.log10(value) : value;
  const t = canLog ? Math.log10(target) : target;
  const p = canLog ? Math.log10(poor) : poor;

  const raw = (v - p) / (t - p);
  return Math.max(0, Math.min(1, raw));
}
