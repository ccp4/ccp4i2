/**
 * Colour-code an assay KPI value against a protocol's absolute thresholds.
 *
 * Thresholds live on the Protocol (target_value / poor_value / threshold_scale).
 * Direction of "better" is implied by the ordering of target_value vs poor_value:
 *   target_value < poor_value  -> lower is better (IC50, EC50, clearance)
 *   target_value > poor_value  -> higher is better (solubility, permeability)
 *
 * When either threshold is null (or both equal), the cell renders uncoloured —
 * this is an honest signal that the protocol owner has not yet curated thresholds,
 * rather than inventing a project-relative scale that would mislead readers.
 */

export interface ProtocolThresholds {
  target_value?: number | null;
  poor_value?: number | null;
  threshold_scale?: 'log' | 'linear' | null;
}

export interface ProtocolColour {
  /** CSS HSL background string, or null if thresholds are unconfigured */
  background: string | null;
  /** Normalised position: 1 = at or beyond target, 0 = at or beyond poor, clamped; null if unconfigured */
  t: number | null;
}

const GREEN_HUE = 140; // slightly bluer than pure 120 — reads as "good" without shouting
const RED_HUE = 5;    // slightly warmer than pure 0
const SATURATION = 65;
const LIGHTNESS = 82; // high lightness keeps text readable on top of the tint

/**
 * Decide the background colour for a single KPI cell.
 * Pure function — safe to call from any render.
 */
export function protocolColour(
  value: number | null | undefined,
  protocol: ProtocolThresholds | null | undefined,
): ProtocolColour {
  if (value == null || !Number.isFinite(value)) {
    return { background: null, t: null };
  }
  if (!protocol) {
    return { background: null, t: null };
  }

  const target = protocol.target_value;
  const poor = protocol.poor_value;
  if (target == null || poor == null || target === poor) {
    return { background: null, t: null };
  }

  const scale = protocol.threshold_scale ?? 'log';

  // Log scale requires all three values to be strictly positive.
  // If any are not, silently fall back to linear rather than returning nothing —
  // unconfigured and invalid should look different.
  const canLog =
    scale === 'log' && value > 0 && target > 0 && poor > 0;

  const v = canLog ? Math.log10(value) : value;
  const t = canLog ? Math.log10(target) : target;
  const p = canLog ? Math.log10(poor) : poor;

  const raw = (v - p) / (t - p);
  const clamped = Math.max(0, Math.min(1, raw));

  const hue = RED_HUE + (GREEN_HUE - RED_HUE) * clamped;

  return {
    background: `hsl(${hue.toFixed(1)}, ${SATURATION}%, ${LIGHTNESS}%)`,
    t: clamped,
  };
}

/**
 * Coarse-grained band label for a normalised position.
 * Useful for legends, filters, and scorecards.
 */
export type ProtocolTier =
  | 'excellent' // at or beyond target
  | 'good'      // top third between target and poor
  | 'mid'       // middle third
  | 'poor'      // bottom third
  | 'failing'   // at or beyond poor
  | null;

export function protocolTier(t: number | null): ProtocolTier {
  if (t == null) return null;
  if (t >= 1) return 'excellent';
  if (t >= 2 / 3) return 'good';
  if (t >= 1 / 3) return 'mid';
  if (t > 0) return 'poor';
  return 'failing';
}
