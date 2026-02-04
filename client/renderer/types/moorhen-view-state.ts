/**
 * Moorhen View State types for URL encoding/decoding.
 *
 * These types support capturing and restoring the Moorhen viewer state
 * via URL query parameters, enabling shareable view links.
 *
 * Keys are abbreviated to minimize URL length when base64 encoded.
 */

/**
 * Compact view state for URL encoding.
 */
export interface MoorhenViewState {
  // Camera state (always present)
  o: [number, number, number]; // origin
  q: [number, number, number, number]; // quaternion
  z: number; // zoom

  // Clip/fog planes (optional - only include if non-default)
  cs?: number; // clipStart (default: 0)
  ce?: number; // clipEnd (default: 1000)
  fs?: number; // fogStart (default: 250)
  fe?: number; // fogEnd (default: 1250)

  // Visibility maps - keyed by fileId (stable ID, not runtime molNo)
  // Only include entries that differ from default (loaded = visible)
  m?: Record<number, MoleculeVisibility>; // molecules
  p?: Record<number, MapVisibility>; // maps
}

export interface MoleculeVisibility {
  v: boolean; // visible
}

export interface MapVisibility {
  v: boolean; // visible
}

/**
 * Default values from Moorhen's glRefSlice.ts initialState.
 * Used to determine which values need to be stored (non-defaults only).
 */
export const MOORHEN_DEFAULTS = {
  clipStart: 0,
  clipEnd: 1000,
  fogStart: 250,
  fogEnd: 1250,
} as const;
