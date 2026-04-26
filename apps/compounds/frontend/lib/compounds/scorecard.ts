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
 * Canonical sectors get fixed colours (no hash collisions possible on the
 * names most projects use); non-canonical sectors hash into a fallback
 * palette that excludes the canonical colours.
 */
const CANONICAL_SECTOR_COLOURS: Record<string, string> = {
  potency: '#0072B2',       // blue
  selectivity: '#D55E00',   // vermillion
  cellular: '#009E73',      // green
  pk: '#CC79A7',            // pink
  'phys-props': '#E69F00',  // orange
  safety: '#56B4E9',        // sky blue
};

/** Fallback palette for non-canonical sectors. Disjoint from the canonical
 *  colours above so a custom sector name never reads as one of the built-ins. */
const FALLBACK_PALETTE: ReadonlyArray<string> = [
  '#F0E442', // yellow
  '#999999', // dark grey
  '#882255', // burgundy
  '#44AA99', // teal
  '#AA4499', // purple
  '#117733', // dark green
];

const NEUTRAL_SECTOR_COLOUR = '#9e9e9e';

/** Common sector names suggested by the editor's autocomplete and given
 *  fixed colours in the palette. Also drives the canonical display order
 *  on radar / bullets (lowest rank first). */
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
 * Deterministic colour for a sector. Canonical sectors get fixed values
 * (no collision risk); non-canonical names hash into FALLBACK_PALETTE.
 */
export function sectorColour(sector: string | null | undefined): string {
  const name = normaliseSector(sector);
  if (!name) return NEUTRAL_SECTOR_COLOUR;
  if (name in CANONICAL_SECTOR_COLOURS) return CANONICAL_SECTOR_COLOURS[name];
  let h = 0;
  for (let i = 0; i < name.length; i++) {
    h = (h * 31 + name.charCodeAt(i)) >>> 0;
  }
  return FALLBACK_PALETTE[h % FALLBACK_PALETTE.length];
}

/**
 * Rank for ordering sectors on the radar and bullets: canonical sectors
 * first, in the order defined above; non-canonical sectors after, ranked
 * alphabetically; unsectored axes last.
 */
function sectorRank(sector: string | null | undefined): number {
  const name = normaliseSector(sector);
  if (!name) return 1_000_000; // trailing
  const canonicalIdx = CANONICAL_SECTORS.indexOf(name);
  if (canonicalIdx >= 0) return canonicalIdx;
  // Non-canonical sectors slot in after canonicals, ordered by name hash
  // (stable across renders, not alphabetical because that would require
  // passing a full sorted sector list through — marginal cost, ignore).
  return CANONICAL_SECTORS.length + name.charCodeAt(0) * 1000 + (name.charCodeAt(1) ?? 0);
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

// ---------------------------------------------------------------------------
// Numeric formatting for scorecard values shown in the bullets / values /
// key tables. Centralised here so both surfaces use the same scientific
// notation rendering — JS's default `toPrecision` returns "1.28e+7", which
// reads as a strikethrough in many fonts (the '+' baseline collides with
// the 'e' descender). Switch to "1.28×10⁷" using Unicode superscripts:
// real codepoints, present in any font, html2canvas captures cleanly.
// ---------------------------------------------------------------------------

const SUPERSCRIPT_DIGIT: Record<string, string> = {
  '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
  '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
};

function digitsToSuperscript(s: string): string {
  let out = '';
  for (const c of s) out += SUPERSCRIPT_DIGIT[c] ?? c;
  return out;
}

export function formatScientific(v: number): string {
  // 3 sig figs (1 + 2 decimals) — matches the precision of the legacy
  // `toPrecision(3)` path it replaces.
  const exp = v.toExponential(2);
  const [mantissa, exponentRaw] = exp.split('e');
  const sign = exponentRaw.startsWith('-') ? '⁻' : '';
  const digits = exponentRaw.replace(/^[+-]/, '');
  return `${mantissa}×10${sign}${digitsToSuperscript(digits)}`;
}

export function formatBareScalar(v: number): string {
  if (!Number.isFinite(v)) return '—';
  if (Math.abs(v) >= 100 || Math.abs(v) < 0.1) return formatScientific(v);
  return v.toFixed(2);
}
