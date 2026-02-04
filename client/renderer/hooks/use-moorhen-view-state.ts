/**
 * Hook for managing Moorhen view state via URL parameters.
 *
 * Provides functionality to:
 * - Apply view state from URL ?view= parameter when files are loaded
 * - Capture current view state for generating shareable URLs
 */

import { useEffect, useRef, useCallback } from "react";
import { useDispatch, useSelector, useStore } from "react-redux";
import { moorhen } from "moorhen/types/moorhen";
import {
  setOrigin,
  setQuat,
  setZoom,
  setClipStart,
  setClipEnd,
  setFogStart,
  setFogEnd,
  setRequestDrawScene,
  hideMolecule,
  hideMap,
} from "moorhen";
import { MoorhenViewState, MOORHEN_DEFAULTS } from "../types/moorhen-view-state";
import {
  decodeViewState,
  captureViewState,
  buildViewStateUrl,
  extractFileIdFromUniqueId,
} from "../lib/moorhen-view-state";

interface UseMoorhenViewStateProps {
  viewParam: string | null;
  onViewRestored?: () => void;
}

interface UseMoorhenViewStateReturn {
  getViewUrl: () => string;
}

export function useMoorhenViewState({
  viewParam,
  onViewRestored,
}: UseMoorhenViewStateProps): UseMoorhenViewStateReturn {
  const dispatch = useDispatch();
  const store = useStore();
  const hasAppliedViewState = useRef(false);

  const cootInitialized = useSelector(
    (state: moorhen.State) => state.generalStates.cootInitialized
  );
  const molecules = useSelector(
    (state: moorhen.State) => state.molecules.moleculeList
  );
  const visibleMolecules = useSelector(
    (state: moorhen.State) => state.molecules.visibleMolecules
  );
  const visibleMaps = useSelector(
    (state: moorhen.State) => state.mapContourSettings.visibleMaps
  );

  // Get all maps from the store
  // Maps are stored at state.maps as an array (mapsSlice)
  const allMaps = useSelector(
    (state: moorhen.State) => (state as unknown as { maps: moorhen.Map[] }).maps || []
  );

  // Apply view state after files are loaded
  useEffect(() => {
    // Skip if no view param, not initialized, or already applied
    if (!viewParam || viewParam.length < 10 || !cootInitialized || hasAppliedViewState.current) {
      return;
    }

    // Wait for at least one molecule or map to be loaded
    // This indicates file loading has progressed
    if (molecules.length === 0 && allMaps.length === 0) {
      return;
    }

    const viewState = decodeViewState(viewParam);
    if (!viewState) {
      console.warn("Failed to decode view state from URL parameter:", viewParam.substring(0, 50));
      // Mark as applied to prevent repeated warnings
      hasAppliedViewState.current = true;
      return;
    }

    // Mark as applied to prevent re-application
    hasAppliedViewState.current = true;

    // Use setTimeout to ensure this runs after the current render cycle
    // and after Moorhen has finished its initial setup
    setTimeout(() => {
      applyViewState(viewState);
      onViewRestored?.();
    }, 500);
  }, [viewParam, cootInitialized, molecules.length, allMaps.length]);

  const applyViewState = useCallback(
    (viewState: MoorhenViewState) => {
      console.log("Applying view state:", viewState);

      // Apply camera state
      dispatch(setOrigin(viewState.o));
      dispatch(setQuat(viewState.q));
      dispatch(setZoom(viewState.z));

      // Apply clip/fog planes (use defaults if not specified)
      dispatch(setClipStart(viewState.cs ?? MOORHEN_DEFAULTS.clipStart));
      dispatch(setClipEnd(viewState.ce ?? MOORHEN_DEFAULTS.clipEnd));
      dispatch(setFogStart(viewState.fs ?? MOORHEN_DEFAULTS.fogStart));
      dispatch(setFogEnd(viewState.fe ?? MOORHEN_DEFAULTS.fogEnd));

      // Apply molecule visibility
      if (viewState.m) {
        for (const mol of molecules) {
          const fileId = extractFileIdFromUniqueId(mol.uniqueId || "");
          if (fileId !== null && viewState.m[fileId]) {
            if (!viewState.m[fileId].v) {
              dispatch(hideMolecule(mol));
            }
          }
        }
      }

      // Apply map visibility
      if (viewState.p) {
        for (const map of allMaps) {
          const fileId = extractFileIdFromUniqueId(map.uniqueId || "");
          if (fileId !== null && viewState.p[fileId]) {
            if (!viewState.p[fileId].v) {
              dispatch(hideMap(map));
            }
          }
        }
      }

      // Trigger redraw
      dispatch(setRequestDrawScene(true));
    },
    [dispatch, molecules, allMaps]
  );

  // Function to capture and return current view URL
  const getViewUrl = useCallback(() => {
    const viewState = captureViewState(
      store as { getState: () => moorhen.State },
      molecules,
      allMaps,
      visibleMolecules,
      visibleMaps
    );
    return buildViewStateUrl(viewState);
  }, [store, molecules, allMaps, visibleMolecules, visibleMaps]);

  return { getViewUrl };
}
