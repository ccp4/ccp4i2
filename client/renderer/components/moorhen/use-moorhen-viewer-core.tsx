"use client";

/**
 * What every CCP4i2 Moorhen viewer page has in common, beneath the drawer:
 * the refs MoorhenContainer fills, the container sizing, the Redux wiring,
 * theme and side-panel effects, the asset prefix and WASM prefetch, the
 * menu system, and the viewer frame itself.
 *
 * The project/job/file viewer (moorhen-wrapper.tsx) and the campaign viewer
 * (campaign-moorhen-wrapper.tsx) were two copies of all of this that had
 * drifted; each now takes it from here and keeps only what is its own. The
 * loaders and scene handlers that sit on top are in use-moorhen-loaders.ts.
 */
import {
  MoorhenContainer,
  MoorhenInstanceProvider,
  MoorhenMenuSystem,
  setBackgroundColor,
  setShownSidePanel,
  setTheme,
} from "moorhen/react-lib";
// @ts-ignore - moorhen 0.23 type may lack .d.ts depending on build
import type { MoorhenInstance, MoorhenPanel } from "moorhen/react-lib";
import { RefObject, useCallback, useEffect, useMemo, useRef, useState } from "react";
import { moorhen } from "moorhen/types/moorhen";
import { useDispatch, useSelector, useStore } from "react-redux";
import { webGL } from "moorhen/types/mgWebGL";
import { useTheme } from "../../theme/theme-provider";
import { isElectronWindow, moorhenUrlPrefix } from "../../lib/moorhen-asset-path";
import { prefetchMoorhenWasm } from "../../lib/moorhen-wasm-prefetch";
import { readCameraState } from "../../lib/moorhen-view-state";
import {
  MoorhenErrorBoundary,
  MoorhenFallback,
  isSafariBrowser,
  useMoorhenCapabilities,
} from "./moorhen-capability-check";
import { usePopcorn } from "../../providers/popcorn-provider";

export const MONOMER_LIBRARY_PATH =
  "https://raw.githubusercontent.com/MonomerLibrary/monomers/master/";

export interface MoorhenViewerCoreOptions {
  /** The key of the one side panel the page registers, opened on mount.
   *  (Moorhen's extension API shows only one tab however many are given;
   *  sub-tabs live inside the panel's own content.) */
  sidePanelKey: string;
}

export type MoorhenViewerCore = ReturnType<typeof useMoorhenViewerCore>;

export function useMoorhenViewerCore({ sidePanelKey }: MoorhenViewerCoreOptions) {
  const capabilities = useMoorhenCapabilities();
  const [isSafari] = useState(() => isSafariBrowser());
  const { setMessage } = usePopcorn();
  const dispatch = useDispatch();
  const theme = useTheme();

  // Container ref for measuring available height below the AppBar.
  // Moorhen 1.0 no longer accepts a setMoorhenDimensions callback; it takes a
  // static `size` prop and otherwise defaults to the full window.innerHeight,
  // which overflows the viewport because the viewer renders under a toolbar.
  // So we measure the container's offset and feed Moorhen an explicit size,
  // keeping it in sync on window resize.
  const moorhenContainerRef = useRef<HTMLDivElement>(null);
  const [size, setSize] = useState<[number, number]>(() =>
    typeof window === "undefined" ? [0, 0] : [window.innerWidth, window.innerHeight],
  );
  useEffect(() => {
    const measure = () => {
      const top = moorhenContainerRef.current?.getBoundingClientRect().top ?? 0;
      setSize([window.innerWidth, window.innerHeight - top]);
    };
    measure();
    window.addEventListener("resize", measure);
    return () => window.removeEventListener("resize", measure);
  }, []);

  const glRef: RefObject<webGL.MGWebGL | null> = useRef(null);
  const commandCentre = useRef<null | moorhen.CommandCentre>(null);
  // Filled by MoorhenContainer on mount. Molecules and maps are built from it.
  const moorhenInstanceRef = useRef<null | MoorhenInstance>(null);
  const moleculesRef = useRef<null | moorhen.Molecule[]>(null);
  const mapsRef = useRef<null | moorhen.Map[]>(null);
  const activeMapRef = useRef<moorhen.Map>(null);
  const lastHoveredAtom = useRef<null | moorhen.HoveredAtom>(null);
  const prevActiveMoleculeRef = useRef<null | moorhen.Molecule>(null);
  const timeCapsuleRef = useRef(null);

  const cootInitialized = useSelector(
    (state: moorhen.State) => state.generalStates.cootInitialized,
  );
  const molecules = useSelector((state: moorhen.State) => state.molecules.moleculeList);
  const maps = useSelector(
    (state: moorhen.State) => (state as unknown as { maps: moorhen.Map[] }).maps || [],
  );
  const backgroundColor = useSelector(
    (state: moorhen.State) => state.sceneSettings.backgroundColor,
  );
  const defaultBondSmoothness = useSelector(
    (state: moorhen.State) => state.sceneSettings.defaultBondSmoothness,
  );
  const store = useStore();

  useEffect(() => {
    dispatch(setBackgroundColor(theme.mode === "light" ? [1, 1, 1, 1] : [0, 0, 0, 1]));
    dispatch(setTheme(theme.mode === "light" ? "flatly" : "darkly"));
  }, [theme.mode, dispatch]);

  useEffect(() => {
    dispatch(setShownSidePanel(sidePanelKey));
  }, [dispatch, sidePanelKey]);

  // URL prefix for Moorhen to load its resources (CSS, pixmaps, monomers, etc.)
  // In web browsers, use API route for CORP headers (COEP compatibility), with
  // the Moorhen version in the path so the route's immutable cache is honest
  // across upgrades. In Electron, serve directly from public/MoorhenAssets.
  const isElectron = isElectronWindow();
  const urlPrefix = moorhenUrlPrefix(isElectron);

  // Start the WASM download now, rather than after the data archives.
  //
  // MoorhenCommandCentre.init() awaits three .tar.gz fetches before posting
  // CootInitialize, and only that message makes the worker fetch the WASM --
  // the biggest file on the page, last in the queue. Warming the cache here
  // lets the two run together; the worker's own fetch then joins it. No effect
  // in Electron, which reads these files from disk.
  useEffect(() => {
    if (isElectron) return;
    return prefetchMoorhenWasm(urlPrefix);
  }, [isElectron, urlPrefix]);

  // Note: Don't subscribe to the camera state here - it changes every frame during rotation
  // and would cause constant re-renders. Access origin directly from store when needed.
  const getOrigin = useCallback(() => {
    return readCameraState(store.getState() as moorhen.State).origin;
  }, [store]);

  // Moorhen 1.0 requires the InstanceProvider to be seeded with a menu system
  // (it builds the per-instance MoorhenInstance from it). One per viewer.
  const menuSystem = useMemo(() => new MoorhenMenuSystem(), []);

  // Show Safari advisory as a non-blocking snackbar
  const isElectronEnv = typeof window !== "undefined" && !!(window as any).electronAPI;
  useEffect(() => {
    if (isSafari && !isElectronEnv) {
      setMessage(
        "Safari may have compatibility issues with the Moorhen viewer. For the best experience, consider using Chrome, Edge, or Firefox.",
        "warning",
      );
    }
  }, [isSafari, isElectronEnv, setMessage]);

  return {
    capabilities,
    setMessage,
    dispatch,
    store,
    moorhenContainerRef,
    size,
    glRef,
    commandCentre,
    moorhenInstanceRef,
    moleculesRef,
    mapsRef,
    activeMapRef,
    lastHoveredAtom,
    prevActiveMoleculeRef,
    timeCapsuleRef,
    cootInitialized,
    molecules,
    maps,
    backgroundColor,
    defaultBondSmoothness,
    monomerLibraryPath: MONOMER_LIBRARY_PATH,
    urlPrefix,
    getOrigin,
    menuSystem,
  };
}

/** The viewer itself: Moorhen inside its error boundary, with the page's
 *  side panel. The props MoorhenContainer takes are memoised on what can
 *  change (the panel, the size), so the container does not re-mount. */
export function MoorhenViewerFrame({
  core,
  extraSidePanels,
}: {
  core: MoorhenViewerCore;
  extraSidePanels: Record<string, MoorhenPanel>;
}) {
  const {
    capabilities,
    store,
    menuSystem,
    moorhenContainerRef,
    size,
    glRef,
    timeCapsuleRef,
    commandCentre,
    moorhenInstanceRef,
    moleculesRef,
    mapsRef,
    activeMapRef,
    lastHoveredAtom,
    prevActiveMoleculeRef,
    monomerLibraryPath,
    urlPrefix,
  } = core;
  const collectedProps = useMemo(
    () => ({
      glRef,
      timeCapsuleRef,
      commandCentre,
      moorhenInstanceRef,
      moleculesRef,
      mapsRef,
      activeMapRef: activeMapRef as React.RefObject<moorhen.Map>,
      lastHoveredAtom,
      prevActiveMoleculeRef,
      monomerLibraryPath,
      urlPrefix,
      store,
      viewOnly: false,
      extraSidePanels,
      size,
    }),
    // The refs are stable for the component's life; only these change.
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [urlPrefix, store, extraSidePanels, size],
  );
  return (
    <div ref={moorhenContainerRef}>
      <MoorhenErrorBoundary
        fallback={<MoorhenFallback reason="runtime_error" capabilities={capabilities} />}
      >
        {store && (
          <MoorhenInstanceProvider menuSystem={menuSystem}>
            <MoorhenContainer {...collectedProps} />
          </MoorhenInstanceProvider>
        )}
      </MoorhenErrorBoundary>
    </div>
  );
}
