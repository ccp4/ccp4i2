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

import type {
  CompactRow,
  MolecularPropertyName,
} from '@/types/compounds/aggregation';
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
  // Every branch is defensive against saved configs that predate the current
  // validator, and against geomeans that arrive as strings from Django decimal
  // fields. Any unexpected shape collapses to null rather than throwing.
  const protocols = compound.protocols ?? {};

  try {
    switch (axis.kind) {
      case 'protocol':
        return asFiniteNumber(
          protocols[(axis as ScorecardProtocolAxis).protocol_id]?.geomean,
        );

      case 'ratio': {
        const { numerator_id, denominator_id } = axis as ScorecardRatioAxis;
        const num = asFiniteNumber(protocols[numerator_id]?.geomean);
        const den = asFiniteNumber(protocols[denominator_id]?.geomean);
        if (num == null || den == null || den === 0) return null;
        return num / den;
      }

      case 'worst_of': {
        const { protocol_ids } = axis as ScorecardWorstOfAxis;
        if (!Array.isArray(protocol_ids) || protocol_ids.length === 0) return null;
        const values: number[] = [];
        for (const id of protocol_ids) {
          const g = asFiniteNumber(protocols[id]?.geomean);
          if (g != null) values.push(g);
        }
        if (values.length === 0) return null;
        const { target_value: target, poor_value: poor } = axis;
        if (target == null || poor == null) return null;
        return target < poor ? Math.max(...values) : Math.min(...values);
      }

      case 'lipinski': {
        const props = compound.properties;
        if (!props) return null;
        let pass = 0;
        const mw = asFiniteNumber(props.molecular_weight);
        const clogp = asFiniteNumber(props.clogp);
        const hbd = asFiniteNumber(props.hbd);
        const hba = asFiniteNumber(props.hba);
        if (mw != null && mw <= 500) pass++;
        if (clogp != null && clogp <= 5) pass++;
        if (hbd != null && hbd <= 5) pass++;
        if (hba != null && hba <= 10) pass++;
        return pass;
      }
    }
  } catch {
    // Unexpected config shape — fail soft so one bad axis doesn't crash the
    // whole card grid.
    return null;
  }
}

/** Coerce number-ish inputs (Decimal strings, null, undefined, NaN) to a
 *  finite number or null. */
function asFiniteNumber(v: unknown): number | null {
  if (v == null) return null;
  const n = typeof v === 'number' ? v : Number(v);
  return Number.isFinite(n) ? n : null;
}

// ---------------------------------------------------------------------------
// Sector helpers — used by the spider, bullets, and editor to give axes a
// consistent visual grouping by drug-discovery theme (potency / pk / etc.).
// ---------------------------------------------------------------------------

/**
 * Hand-picked categorical palette tuned for colour-blindness (Okabe-Ito).
 * Sectors are assigned colours deterministically by name, so the same
 * sector reads the same colour across every spider, every card, and every
 * bullets banner — chemists learn the mapping once.
 */
const SECTOR_PALETTE: ReadonlyArray<string> = [
  '#0072B2', // blue
  '#D55E00', // vermillion
  '#009E73', // green
  '#CC79A7', // pink
  '#F0E442', // yellow
  '#56B4E9', // sky blue
  '#E69F00', // orange
];

const NEUTRAL_SECTOR_COLOUR = '#9e9e9e';

/** Common sector names suggested by the editor's autocomplete. Free-string
 *  field, so projects can invent their own — these are just hints. */
export const CANONICAL_SECTORS: ReadonlyArray<string> = [
  'potency',
  'selectivity',
  'cellular',
  'pk',
  'phys-props',
  'safety',
];

function normaliseSector(sector: string | null | undefined): string | null {
  if (sector == null) return null;
  const trimmed = sector.trim().toLowerCase();
  return trimmed === '' ? null : trimmed;
}

/**
 * Deterministic colour assignment for a sector name. Runs the normalised
 * name through a simple hash and indexes into the palette, so two configs
 * that use the same sector name get the same colour without coordination.
 */
export function sectorColour(sector: string | null | undefined): string {
  const name = normaliseSector(sector);
  if (!name) return NEUTRAL_SECTOR_COLOUR;
  let h = 0;
  for (let i = 0; i < name.length; i++) {
    h = (h * 31 + name.charCodeAt(i)) >>> 0;
  }
  return SECTOR_PALETTE[h % SECTOR_PALETTE.length];
}

/**
 * Return the axes ordered so that axes sharing a sector are adjacent.
 * Within a sector, axes keep their relative order from the source array.
 * Axes with no sector go to the end as an "(unsectored)" trailing group.
 *
 * Returned indices preserve a back-reference to the original positions so
 * that callers can map evaluations / values without re-keying.
 */
export interface SectorGroup<T> {
  sector: string | null;
  /** Pairs of [original index, item] in the source array. */
  items: Array<{ index: number; item: T }>;
}

export function groupAxesBySector<T extends { sector?: string | null }>(
  axes: ReadonlyArray<T>,
): SectorGroup<T>[] {
  const order: string[] = [];
  const buckets = new Map<string, Array<{ index: number; item: T }>>();
  const NONE = '__none__';

  axes.forEach((item, index) => {
    const key = normaliseSector(item.sector) ?? NONE;
    if (!buckets.has(key)) {
      buckets.set(key, []);
      order.push(key);
    }
    buckets.get(key)!.push({ index, item });
  });

  // Stable sort: encountered-order, with the unsectored bucket pushed last.
  const sectoredKeys = order.filter((k) => k !== NONE);
  const tail = order.includes(NONE) ? [NONE] : [];
  return [...sectoredKeys, ...tail].map((key) => ({
    sector: key === NONE ? null : key,
    items: buckets.get(key)!,
  }));
}

/**
 * Flatten the sector groups back into a single array in display order
 * (sectors grouped, axes within sectors keeping source order).
 */
export function sortAxesBySector<T extends { sector?: string | null }>(
  axes: ReadonlyArray<T>,
): Array<{ index: number; item: T }> {
  return groupAxesBySector(axes).flatMap((g) => g.items);
}

/**
 * Aggregation-request inputs a scorecard depends on. The aggregation page
 * merges these with the user's explicit selections so chemists don't have
 * to manually tick boxes for each protocol referenced by the scorecard.
 */
export interface ScorecardDataNeeds {
  protocolIds: string[];
  properties: MolecularPropertyName[];
}

const LIPINSKI_PROPERTIES: MolecularPropertyName[] = [
  'molecular_weight',
  'clogp',
  'hbd',
  'hba',
];

export function scorecardDataNeeds(
  config: ScorecardConfig | null | undefined,
): ScorecardDataNeeds {
  const protocolIds = new Set<string>();
  let needsLipinski = false;

  for (const axis of config?.axes ?? []) {
    switch (axis.kind) {
      case 'protocol':
        if (axis.protocol_id) protocolIds.add(axis.protocol_id);
        break;
      case 'ratio':
        if (axis.numerator_id) protocolIds.add(axis.numerator_id);
        if (axis.denominator_id) protocolIds.add(axis.denominator_id);
        break;
      case 'worst_of':
        for (const id of axis.protocol_ids ?? []) {
          if (id) protocolIds.add(id);
        }
        break;
      case 'lipinski':
        needsLipinski = true;
        break;
    }
  }

  return {
    protocolIds: Array.from(protocolIds),
    properties: needsLipinski ? LIPINSKI_PROPERTIES : [],
  };
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
