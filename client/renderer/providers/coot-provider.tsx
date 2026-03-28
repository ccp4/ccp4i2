/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
    // Capability check performed on mount

    // Warn if capabilities are missing (but don't block - let it try)
    if (!capabilities.isSupported && !isElectron) {
      console.warn(
        `[Moorhen] ${capabilities.message}. Attempting to load anyway...`
      );
    }
  }, []);

  return <>{props.children}</>;
};
