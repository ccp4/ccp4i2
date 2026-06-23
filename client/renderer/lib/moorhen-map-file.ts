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
  setMapRadius,
  setContourLevel,
  setMapColours,
} from "moorhen/react-lib";

/** File.sub_type marking a CCP4-map file as a mask (CMapDataFile.SUBTYPE_MASK). */
export const MASK_SUBTYPE = 4;

// --------------------------------------------------------------------------
// CCP4 map mode-0 -> float conversion
//
// dm's masks are written as mode-0 (signed int8) CCP4 maps — the natural mask
// grid, and what dm reads back as NCSIN. coot/Moorhen, though, read mode-0
// maps with degenerate statistics (zero range), which makes the map card crash
// (toFixed on a bad digit count) and contour pathologically. Rather than change
// the files dm consumes, we convert mode-0 -> mode-2 (float32) in the browser
// before handing the bytes to Moorhen: a 0.0/1.0 float map with real header
// stats renders cleanly. Other modes (incl. mode-2 float) pass through
// untouched.
// --------------------------------------------------------------------------

// CCP4 map header layout (256 4-byte words). Byte offsets we touch:
const OFF_NC = 0, OFF_NR = 4, OFF_NS = 8, OFF_MODE = 12;
const OFF_DMIN = 76, OFF_DMAX = 80, OFF_DMEAN = 84, OFF_NSYMBT = 92;
const OFF_MACHST = 212, OFF_RMS = 216;
const CCP4_HEADER_BYTES = 1024;

/** CCP4 map MODE word (0=int8, 1=int16, 2=float32, ...), or null if too short. */
export function ccp4MapMode(buffer: ArrayBuffer, littleEndian = true): number | null {
  if (buffer.byteLength < CCP4_HEADER_BYTES) return null;
  return new DataView(buffer).getInt32(OFF_MODE, littleEndian);
}

/**
 * Convert a mode-0 (signed int8) CCP4 map to mode-2 (float32), recomputing the
 * DMIN/DMAX/DMEAN/RMS header stats. Returns a new ArrayBuffer. Any non-mode-0
 * input (or a too-short buffer) is returned unchanged.
 */
export function ccp4Mode0ToFloat(buffer: ArrayBuffer): ArrayBuffer {
  if (buffer.byteLength < CCP4_HEADER_BYTES) return buffer;
  const src = new DataView(buffer);
  // Endianness from the machine stamp (BE stamp starts 0x11); default little.
  const le = src.getUint8(OFF_MACHST) !== 0x11;
  if (src.getInt32(OFF_MODE, le) !== 0) return buffer; // only mode-0

  const nc = src.getInt32(OFF_NC, le);
  const nr = src.getInt32(OFF_NR, le);
  const ns = src.getInt32(OFF_NS, le);
  const nsymbt = src.getInt32(OFF_NSYMBT, le);
  const n = nc * nr * ns;
  const dataStart = CCP4_HEADER_BYTES + nsymbt;
  if (n <= 0 || dataStart + n > buffer.byteLength) return buffer; // malformed

  // Read int8 values -> float32, accumulating stats.
  const i8 = new Int8Array(buffer, dataStart, n);
  const floats = new Float32Array(n);
  let min = Infinity, max = -Infinity, sum = 0, sumSq = 0;
  for (let i = 0; i < n; i++) {
    const v = i8[i];
    floats[i] = v;
    if (v < min) min = v;
    if (v > max) max = v;
    sum += v;
    sumSq += v * v;
  }
  const mean = n ? sum / n : 0;
  const rms = n ? Math.sqrt(Math.max(0, sumSq / n - mean * mean)) : 0;

  // New buffer: header + symmetry records (copied) + float32 data.
  const out = new ArrayBuffer(dataStart + n * 4);
  const outBytes = new Uint8Array(out);
  outBytes.set(new Uint8Array(buffer, 0, dataStart)); // header + symmetry verbatim
  const dst = new DataView(out);
  dst.setInt32(OFF_MODE, 2, le);
  dst.setFloat32(OFF_DMIN, min === Infinity ? 0 : min, le);
  dst.setFloat32(OFF_DMAX, max === -Infinity ? 0 : max, le);
  dst.setFloat32(OFF_DMEAN, mean, le);
  dst.setFloat32(OFF_RMS, rms, le);
  // Write float data in the file's byte order.
  for (let i = 0; i < n; i++) {
    dst.setFloat32(dataStart + i * 4, floats[i], le);
  }
  return out;
}

/** Default mask look: a translucent solid surface in a soft blue. NB Moorhen's
 *  setMapColours expects 0-255 components (it divides by 255 internally) — 0-1
 *  floats render as a near-black mesh. */
export const MASK_STYLE = "solid" as const;
export const MASK_ALPHA = 0.45;
export const MASK_COLOUR = "#7e9cd8";
export const MASK_COLOUR_RGB = { r: 126, g: 156, b: 216 };

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
  };
  m.isCcp4Mask = true;
  m.mapSubType = MASK_SUBTYPE;
}

/**
 * Apply the default mask look to a freshly-loaded Moorhen map.
 *
 * Crucially the contour LEVEL and RADIUS are taken from Coot's own
 * getSuggestedSettings(), not hardcoded: Coot reads contour level in sigma
 * units, and a mask (mostly 0, with 1s) has a tiny sigma, so a "sensible"
 * absolute level clears nothing -> empty mesh -> undefined bounding sphere ->
 * crash. We override only cosmetics (style / alpha / colour), then force a
 * contour + colour flush so the GL buffers actually pick them up.
 *
 * (Lesson set distilled from the FragVol sub-volume work.)
 */
export async function applyMaskDefaults(
  dispatch: Dispatch,
  map: { molNo: number; suggestedContourLevel?: number; suggestedRadius?: number;
         getSuggestedSettings?: () => Promise<void>;
         drawMapContour?: () => Promise<void>;
         fetchColourAndRedraw?: () => Promise<void>; },
): Promise<void> {
  const molNo = map.molNo;
  try {
    await map.getSuggestedSettings?.();
  } catch {
    /* fall back to whatever level Coot already has */
  }
  if (typeof map.suggestedContourLevel === "number") {
    dispatch(setContourLevel({ molNo, contourLevel: map.suggestedContourLevel } as never));
  }
  if (typeof map.suggestedRadius === "number") {
    dispatch(setMapRadius({ molNo, radius: map.suggestedRadius } as never));
  }
  // Cosmetics only.
  dispatch(setMapStyle({ molNo, style: MASK_STYLE } as never));
  dispatch(setMapAlpha({ molNo, alpha: MASK_ALPHA } as never));
  dispatch(setMapColours({ molNo, rgb: MASK_COLOUR_RGB } as never));
  // The colour/level don't reach the GL buffers until a contour runs.
  try {
    await map.drawMapContour?.();
    await map.fetchColourAndRedraw?.();
  } catch {
    /* best-effort flush */
  }
}
