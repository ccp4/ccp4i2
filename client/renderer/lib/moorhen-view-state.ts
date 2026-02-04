/**
 * Moorhen View State utilities for URL encoding/decoding.
 *
 * Provides functions to capture, encode, decode, and apply Moorhen viewer
 * state via URL query parameters.
 */

import {
  MoorhenViewState,
  MOORHEN_DEFAULTS,
} from "../types/moorhen-view-state";
import { moorhen } from "moorhen/types/moorhen";

/**
 * Extract file ID from a molecule/map's uniqueId URL.
 * Handles both /api/proxy/files/123/download/ and /api/proxy/ccp4i2/files/123/download/
 */
export function extractFileIdFromUniqueId(uniqueId: string): number | null {
  // Match any path containing /files/123/download/
  const match = uniqueId.match(/\/files\/(\d+)\/download\//);
  return match ? parseInt(match[1], 10) : null;
}

/**
 * Encode view state to a compact base64url string for URL.
 * Uses base64url encoding which is URL-safe (no + / = characters).
 */
export function encodeViewState(state: MoorhenViewState): string {
  // Validate input state
  if (!state || !Array.isArray(state.o) || !Array.isArray(state.q)) {
    console.error("encodeViewState - Invalid input state:", state);
    throw new Error("Invalid view state: missing origin or quaternion");
  }

  // Round numbers to reduce URL size while maintaining precision
  const compacted: MoorhenViewState = {
    o: state.o.map((v) => Math.round(v * 1000) / 1000) as [
      number,
      number,
      number
    ],
    q: state.q.map((v) => Math.round(v * 10000) / 10000) as [
      number,
      number,
      number,
      number
    ],
    z: Math.round(state.z * 1000) / 1000,
  };

  // Only include non-default clip/fog values
  if (state.cs !== undefined) compacted.cs = state.cs;
  if (state.ce !== undefined) compacted.ce = state.ce;
  if (state.fs !== undefined) compacted.fs = state.fs;
  if (state.fe !== undefined) compacted.fe = state.fe;

  // Include visibility maps if present
  if (state.m && Object.keys(state.m).length > 0) compacted.m = state.m;
  if (state.p && Object.keys(state.p).length > 0) compacted.p = state.p;

  const json = JSON.stringify(compacted);
  console.log("encodeViewState - JSON to encode:", json);

  // Use base64url encoding (URL-safe: replace + with -, / with _, remove padding)
  const encoded = btoa(json).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
  console.log("encodeViewState - Encoded result:", encoded);

  return encoded;
}

/**
 * Decode view state from a URL parameter.
 * Returns null if decoding fails or validation fails.
 */
export function decodeViewState(encoded: string): MoorhenViewState | null {
  try {
    console.log("decodeViewState - Input (first 100 chars):", encoded.substring(0, 100));

    // Restore base64 padding and characters
    let base64 = encoded.replace(/-/g, "+").replace(/_/g, "/");
    while (base64.length % 4) base64 += "=";

    console.log("decodeViewState - Restored base64 (first 100 chars):", base64.substring(0, 100));

    const json = atob(base64);
    console.log("decodeViewState - Decoded JSON:", json);

    const parsed = JSON.parse(json);
    console.log("decodeViewState - Parsed object:", parsed);

    // Validate required fields
    if (!Array.isArray(parsed.o) || parsed.o.length !== 3) {
      console.error("decodeViewState - Invalid origin:", parsed.o);
      return null;
    }
    if (!Array.isArray(parsed.q) || parsed.q.length !== 4) {
      console.error("decodeViewState - Invalid quaternion:", parsed.q);
      return null;
    }
    if (typeof parsed.z !== "number") {
      console.error("decodeViewState - Invalid zoom:", parsed.z);
      return null;
    }

    // Validate array elements are numbers
    if (!parsed.o.every((v: unknown) => typeof v === "number")) {
      console.error("decodeViewState - Origin contains non-numbers:", parsed.o);
      return null;
    }
    if (!parsed.q.every((v: unknown) => typeof v === "number")) {
      console.error("decodeViewState - Quaternion contains non-numbers:", parsed.q);
      return null;
    }

    console.log("decodeViewState - Success, returning:", parsed);
    return parsed as MoorhenViewState;
  } catch (e) {
    console.error("decodeViewState - Exception:", e);
    return null;
  }
}

/**
 * Capture current view state from Moorhen Redux store.
 *
 * @param store - Redux store instance
 * @param molecules - Array of loaded molecules
 * @param maps - Array of loaded maps
 * @param visibleMolecules - Array of visible molecule molNo values
 * @param visibleMaps - Array of visible map molNo values
 */
export function captureViewState(
  store: { getState: () => moorhen.State },
  molecules: moorhen.Molecule[],
  maps: moorhen.Map[],
  visibleMolecules: number[],
  visibleMaps: number[]
): MoorhenViewState {
  const state = store.getState();
  // glRef is at the root level of the Redux store
  const glRef = (state as unknown as { glRef: {
    origin: number[];
    quat: number[];
    zoom: number;
    clipStart: number;
    clipEnd: number;
    fogStart: number;
    fogEnd: number;
  }}).glRef;

  // Debug logging
  console.log("captureViewState - glRef:", {
    origin: glRef?.origin,
    quat: glRef?.quat,
    zoom: glRef?.zoom,
    fogStart: glRef?.fogStart,
    fogEnd: glRef?.fogEnd,
    clipStart: glRef?.clipStart,
    clipEnd: glRef?.clipEnd,
  });

  // Helper to check for array-like objects (includes Float32Array and other typed arrays)
  // Note: Array.isArray(Float32Array) returns false, so we check for length property
  const isArrayLike = (obj: unknown): boolean =>
    obj !== null &&
    typeof obj === "object" &&
    typeof (obj as { length?: unknown }).length === "number";

  // Validate glRef data - accept both regular arrays and typed arrays
  if (!glRef || !isArrayLike(glRef.origin) || !isArrayLike(glRef.quat)) {
    console.error("captureViewState - Invalid glRef state:", glRef);
    // Return default state if glRef is invalid
    return {
      o: [0, 0, 0],
      q: [0, 0, 0, -1],
      z: 1.0,
    };
  }

  // Convert typed arrays (Float32Array) to plain arrays for JSON serialization
  // Typed arrays serialize as objects {"0": val, "1": val} instead of arrays
  const viewState: MoorhenViewState = {
    o: Array.from(glRef.origin).slice(0, 3) as [number, number, number],
    q: Array.from(glRef.quat).slice(0, 4) as [number, number, number, number],
    z: glRef.zoom,
  };

  // Only include clip/fog if defined and different from defaults
  if (glRef.clipStart !== undefined && glRef.clipStart !== MOORHEN_DEFAULTS.clipStart) {
    viewState.cs = glRef.clipStart;
  }
  if (glRef.clipEnd !== undefined && glRef.clipEnd !== MOORHEN_DEFAULTS.clipEnd) {
    viewState.ce = glRef.clipEnd;
  }
  if (glRef.fogStart !== undefined && glRef.fogStart !== MOORHEN_DEFAULTS.fogStart) {
    viewState.fs = glRef.fogStart;
  }
  if (glRef.fogEnd !== undefined && glRef.fogEnd !== MOORHEN_DEFAULTS.fogEnd) {
    viewState.fe = glRef.fogEnd;
  }

  // Capture molecule visibility (only hidden items - default is visible)
  const moleculeVisibility: Record<number, { v: boolean }> = {};
  for (const mol of molecules) {
    const fileId = extractFileIdFromUniqueId(mol.uniqueId || "");
    if (fileId !== null) {
      const isVisible = visibleMolecules.includes(mol.molNo);
      // Default is visible (loaded = visible), so only record if hidden
      if (!isVisible) {
        moleculeVisibility[fileId] = { v: false };
      }
    }
  }
  if (Object.keys(moleculeVisibility).length > 0) {
    viewState.m = moleculeVisibility;
  }

  // Capture map visibility (only hidden items)
  const mapVisibility: Record<number, { v: boolean }> = {};
  for (const map of maps) {
    const fileId = extractFileIdFromUniqueId(map.uniqueId || "");
    if (fileId !== null) {
      const isVisible = visibleMaps.includes(map.molNo);
      if (!isVisible) {
        mapVisibility[fileId] = { v: false };
      }
    }
  }
  if (Object.keys(mapVisibility).length > 0) {
    viewState.p = mapVisibility;
  }

  return viewState;
}

/**
 * Build the full URL with view state parameter.
 * Uses the current window location as the base.
 */
export function buildViewStateUrl(viewState: MoorhenViewState): string {
  const encoded = encodeViewState(viewState);
  const url = new URL(window.location.href);
  url.searchParams.set("view", encoded);
  return url.toString();
}
