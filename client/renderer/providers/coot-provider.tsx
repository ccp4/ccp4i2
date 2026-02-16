import { PropsWithChildren, useEffect } from "react";
import { checkMoorhenCapabilities } from "../components/moorhen/moorhen-capability-check";

// Check if running in Electron (no COEP restrictions) vs web browser
const isElectron = typeof window !== "undefined" && !!(window as any).electronAPI;

/**
 * CootProvider - Moorhen 0.23+
 *
 * In Moorhen 0.23, the MoorhenContainer handles its own WASM module loading
 * internally via windowCootCCP4Loader and CommandCentre worker initialization.
 * This provider only performs capability checking; WASM loading is no longer
 * needed here.
 */
export const CootProvider: React.FC<PropsWithChildren> = (props) => {
  // Log capability check on mount
  useEffect(() => {
    const capabilities = checkMoorhenCapabilities();
    console.log(`[Moorhen] Browser capabilities:`, capabilities);

    // Warn if capabilities are missing (but don't block - let it try)
    if (!capabilities.isSupported && !isElectron) {
      console.warn(
        `[Moorhen] ${capabilities.message}. Attempting to load anyway...`
      );
    }
  }, []);

  return <>{props.children}</>;
};
