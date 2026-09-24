/**
 * Flattening Moorhen's map-contour redux slice for the scene lifter.
 *
 * `mapContourSettings` stores one array per attribute, each keyed by molNo;
 * the lifter wants one object per map. Shared by every wrapper that can
 * capture a scene, so a Moorhen release that moves an attribute is one edit.
 */
import type { MapRenderState } from "./moorhen-scene-lifter";

interface ContourEntry {
  molNo: number;
  contourLevel?: number;
  radius?: number;
  alpha?: number;
  style?: "lines" | "solid" | "lit-lines";
  rgb?: { r: number; g: number; b: number };
}

/**
 * Flatten Moorhen's mapContourSettings slice — which stores per-attribute
 * arrays keyed by molNo — into a single MapRenderState per molNo for the
 * lifter to consume. Defensive: every field is optional, and we tolerate
 * the slice being absent (older Moorhen versions, or no maps loaded).
 */
export function collectMapRenderState(
  state: unknown,
): Record<number, MapRenderState> {
  const out: Record<number, MapRenderState> = {};
  const slice = (state as { mapContourSettings?: Record<string, unknown> })
    .mapContourSettings;
  if (!slice) return out;
  const ensure = (molNo: number): MapRenderState => {
    if (!out[molNo]) out[molNo] = {};
    return out[molNo];
  };
  const merge = (
    rows: unknown,
    setter: (s: MapRenderState, v: ContourEntry) => void,
  ) => {
    if (!Array.isArray(rows)) return;
    for (const row of rows as ContourEntry[]) {
      if (typeof row?.molNo === "number") setter(ensure(row.molNo), row);
    }
  };
  merge(slice.contourLevels, (s, v) => { if (v.contourLevel !== undefined) s.contourLevel = v.contourLevel; });
  merge(slice.mapRadii, (s, v) => { if (v.radius !== undefined) s.radius = v.radius; });
  merge(slice.mapAlpha, (s, v) => { if (v.alpha !== undefined) s.alpha = v.alpha; });
  merge(slice.mapStyles, (s, v) => { if (v.style) s.style = v.style; });
  merge(slice.mapColours, (s, v) => { if (v.rgb) s.colour = rgb01ToHex(v.rgb); });
  merge(slice.positiveMapColours, (s, v) => { if (v.rgb) s.positiveColour = rgb01ToHex(v.rgb); });
  merge(slice.negativeMapColours, (s, v) => { if (v.rgb) s.negativeColour = rgb01ToHex(v.rgb); });
  const visible = (slice as { visibleMaps?: unknown }).visibleMaps;
  if (Array.isArray(visible)) {
    const visibleSet = new Set<number>(
      (visible as unknown[]).filter((n): n is number => typeof n === "number"),
    );
    // visibleMaps lists the *visible* molNos; map state already includes
    // every molNo we've seen, but make sure we don't drop a map that has
    // no other contour settings touched.
    for (const molNo of visibleSet) ensure(molNo);
    for (const molNo of Object.keys(out).map(Number)) {
      out[molNo].visible = visibleSet.has(molNo);
    }
  }
  return out;
}

/** Convert Moorhen's {r,g,b} 0-1 shape to a 6-hex string. */
function rgb01ToHex(rgb: { r: number; g: number; b: number }): string {
  const to8 = (v: number) =>
    Math.max(0, Math.min(255, Math.round(v * 255)))
      .toString(16)
      .padStart(2, "0");
  return `#${to8(rgb.r)}${to8(rgb.g)}${to8(rgb.b)}`;
}
