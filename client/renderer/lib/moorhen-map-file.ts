/**
 * Helpers for rendering real-space CCP4 map files (application/CCP4-map) in
 * Moorhen — including masks, which are the same FileType distinguished only by
 * File.sub_type (CMapDataFile.SUBTYPE_MASK).
 *
 * Real-space maps load via MoorhenMap.loadToCootFromMapData (not the MTZ
 * loadToCootFromMtzData path). A mask reads best as a translucent solid
 * surface so it shows the *region* it covers over the model; these defaults
 * are shared by the job/file/campaign viewers and the scene resolver so a
 * mask looks the same however it's loaded.
 */
import type { Dispatch } from "redux";
import {
  setMapStyle,
  setMapAlpha,
  setContourLevel,
  setMapColours,
} from "moorhen/react-lib";

/** File.sub_type marking a CCP4-map file as a mask (CMapDataFile.SUBTYPE_MASK). */
export const MASK_SUBTYPE = 4;

/** Default mask look: a translucent solid surface in a soft blue. */
export const MASK_STYLE = "solid" as const;
export const MASK_ALPHA = 0.4;
export const MASK_CONTOUR_LEVEL = 0.5;
export const MASK_COLOUR = "#7e9cd8";
export const MASK_COLOUR_RGB = { r: 0.494, g: 0.612, b: 0.847 };

/** True if a DB file's sub_type marks it as a mask. */
export function isMaskSubType(subType: number | null | undefined): boolean {
  return subType === MASK_SUBTYPE;
}

/**
 * Tag a freshly-loaded MoorhenMap as a mask. Sets:
 *  - `isCcp4Mask` — read by the lifter (with `isCcp4MapFile`) to emit isMask;
 *  - `mapSubType` — so the contour-slider label reads "Mask" (not "2Fo-Fc");
 *  - `isEM` — so Moorhen contours it in absolute (EM-style) mode rather than
 *    rmsd-relative: a mask holds 0..1 region values, not crystallographic
 *    density, so the rmsd path gives a degenerate histogram / contour range.
 */
export function markMaskMap(map: unknown): void {
  const m = map as {
    isCcp4Mask?: boolean;
    mapSubType?: number;
    isEM?: boolean;
  };
  m.isCcp4Mask = true;
  m.mapSubType = MASK_SUBTYPE;
  m.isEM = true;
}

/**
 * Apply the default mask look to a freshly-loaded Moorhen map. Dispatched as
 * separate actions (mirroring the resolver's applyMapState) because the
 * Moorhen map-setting action typings have drifted across versions — the casts
 * keep this compatible with the installed shape.
 */
export function applyMaskDefaults(dispatch: Dispatch, molNo: number): void {
  dispatch(setMapStyle({ molNo, style: MASK_STYLE } as never));
  dispatch(setMapAlpha({ molNo, alpha: MASK_ALPHA } as never));
  dispatch(setContourLevel({ molNo, contourLevel: MASK_CONTOUR_LEVEL } as never));
  dispatch(setMapColours({ molNo, rgb: MASK_COLOUR_RGB } as never));
}
