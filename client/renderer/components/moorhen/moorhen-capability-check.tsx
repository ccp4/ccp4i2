"use client";

import React, { useState, useEffect, useCallback } from "react";

/**
 * Browser capability check results for Moorhen/WebAssembly threading.
 */
export interface MoorhenCapabilities {
  /** Page has Cross-Origin-Embedder-Policy set correctly */
  crossOriginIsolated: boolean;
  /** SharedArrayBuffer is available (required for WASM threading) */
  sharedArrayBufferAvailable: boolean;
  /** Overall: can Moorhen likely work? */
  isSupported: boolean;
  /** Descriptive message about the capability status */
  message: string;
}

/**
 * Check browser capabilities required for Moorhen's WebAssembly threading.
 *
 * Instead of blocking based on user-agent, we check actual capabilities:
 * - crossOriginIsolated: Required for SharedArrayBuffer in secure contexts
 * - SharedArrayBuffer: Required for WASM threading
 *
 * This allows browsers that have improved their support to work, while still
 * providing helpful feedback when things won't work.
 */
export function checkMoorhenCapabilities(): MoorhenCapabilities {
  if (typeof window === "undefined") {
    return {
      crossOriginIsolated: false,
      sharedArrayBufferAvailable: false,
      isSupported: false,
      message: "Server-side rendering - capabilities unknown",
    };
  }

  const crossOriginIsolated = window.crossOriginIsolated === true;
  const sharedArrayBufferAvailable = typeof SharedArrayBuffer !== "undefined";

  // Both are required for Moorhen to work
  const isSupported = crossOriginIsolated && sharedArrayBufferAvailable;

  let message: string;
  if (isSupported) {
    message = "Browser supports WebAssembly threading";
  } else if (!crossOriginIsolated && !sharedArrayBufferAvailable) {
    message = "Cross-origin isolation and SharedArrayBuffer are not available";
  } else if (!crossOriginIsolated) {
    message = "Cross-origin isolation is not enabled for this page";
  } else {
    message = "SharedArrayBuffer is not available in this browser";
  }

  return {
    crossOriginIsolated,
    sharedArrayBufferAvailable,
    isSupported,
    message,
  };
}

/**
 * Detect Safari browser (known to have WASM threading issues with Moorhen).
 */
export const isSafariBrowser = (): boolean => {
  if (typeof navigator === "undefined") return false;
  const ua = navigator.userAgent;
  const isSafari = /^((?!chrome|android).)*safari/i.test(ua);
  const isIOS = /iPad|iPhone|iPod/.test(ua);
  return isSafari || isIOS;
};

/**
 * Safari experimental warning - shown when Safari is detected.
 * Capabilities may look OK but WASM threading is known to crash.
 */
export interface SafariExperimentalWarningProps {
  onProceed: () => void;
}

export const SafariExperimentalWarning: React.FC<SafariExperimentalWarningProps> = ({
  onProceed,
}) => (
  <div
    style={{
      display: "flex",
      flexDirection: "column",
      alignItems: "center",
      justifyContent: "center",
      height: "100%",
      padding: "40px",
      textAlign: "center",
      backgroundColor: "#fff3cd",
      border: "1px solid #ffc107",
      borderRadius: "8px",
      margin: "20px",
    }}
  >
    <h2 style={{ color: "#856404", marginBottom: "16px" }}>
      Safari Has Known Compatibility Issues
    </h2>
    <div style={{ color: "#856404", maxWidth: "600px", lineHeight: "1.6" }}>
      <p style={{ marginBottom: "12px" }}>
        The Moorhen molecular viewer uses WebAssembly threading which has known
        issues in Safari. While your browser reports the required features as
        available, Safari&apos;s implementation often crashes during use.
      </p>
      <p style={{ marginBottom: "16px" }}>
        For the best experience, we recommend using <strong>Google Chrome</strong>,{" "}
        <strong>Microsoft Edge</strong>, or <strong>Firefox</strong>.
      </p>
      <div style={{ display: "flex", gap: "12px", justifyContent: "center", flexWrap: "wrap" }}>
        <button
          onClick={onProceed}
          style={{
            padding: "10px 20px",
            backgroundColor: "#856404",
            color: "white",
            border: "none",
            borderRadius: "4px",
            cursor: "pointer",
            fontSize: "14px",
          }}
        >
          Try Anyway (Experimental)
        </button>
      </div>
      <p style={{ marginTop: "16px", fontSize: "12px", opacity: 0.8 }}>
        If you proceed, you may experience crashes or errors when loading molecules.
      </p>
    </div>
  </div>
);

/**
 * Props for the MoorhenFallback component.
 */
export interface MoorhenFallbackProps {
  /** Why Moorhen couldn't load */
  reason: "capabilities" | "load_error" | "runtime_error";
  /** Additional error details */
  error?: Error | string;
  /** Capability check results (for capabilities reason) */
  capabilities?: MoorhenCapabilities;
}

/**
 * Fallback UI shown when Moorhen cannot load.
 *
 * Provides helpful information about why Moorhen isn't available
 * and what the user can do about it.
 */
export const MoorhenFallback: React.FC<MoorhenFallbackProps> = ({
  reason,
  error,
  capabilities,
}) => {
  const getTitle = () => {
    switch (reason) {
      case "capabilities":
        return "Browser Compatibility Issue";
      case "load_error":
        return "Failed to Load Molecular Viewer";
      case "runtime_error":
        return "Molecular Viewer Error";
    }
  };

  const getMessage = () => {
    switch (reason) {
      case "capabilities":
        return (
          <>
            <p style={{ marginBottom: "12px" }}>
              The Moorhen molecular viewer requires WebAssembly threading features
              that are not available in your current browser configuration.
            </p>
            {capabilities && (
              <ul style={{ textAlign: "left", marginBottom: "12px", paddingLeft: "20px" }}>
                <li>
                  Cross-Origin Isolation:{" "}
                  {capabilities.crossOriginIsolated ? "✓ Enabled" : "✗ Not enabled"}
                </li>
                <li>
                  SharedArrayBuffer:{" "}
                  {capabilities.sharedArrayBufferAvailable ? "✓ Available" : "✗ Not available"}
                </li>
              </ul>
            )}
          </>
        );
      case "load_error":
        return (
          <p style={{ marginBottom: "12px" }}>
            The molecular viewer failed to initialize. This may be due to browser
            compatibility issues or network problems loading required resources.
          </p>
        );
      case "runtime_error":
        return (
          <p style={{ marginBottom: "12px" }}>
            An error occurred while running the molecular viewer.
          </p>
        );
    }
  };

  const errorDetails = error
    ? typeof error === "string"
      ? error
      : error.message
    : null;

  return (
    <div
      style={{
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        justifyContent: "center",
        height: "100%",
        padding: "40px",
        textAlign: "center",
        backgroundColor: "#fff3cd",
        border: "1px solid #ffc107",
        borderRadius: "8px",
        margin: "20px",
      }}
    >
      <h2 style={{ color: "#856404", marginBottom: "16px" }}>{getTitle()}</h2>
      <div style={{ color: "#856404", maxWidth: "600px", lineHeight: "1.6" }}>
        {getMessage()}

        {errorDetails && (
          <details style={{ marginTop: "12px", textAlign: "left" }}>
            <summary style={{ cursor: "pointer" }}>Technical details</summary>
            <pre
              style={{
                marginTop: "8px",
                padding: "8px",
                backgroundColor: "#fff8e1",
                borderRadius: "4px",
                fontSize: "12px",
                whiteSpace: "pre-wrap",
                wordBreak: "break-word",
              }}
            >
              {errorDetails}
            </pre>
          </details>
        )}

        <div style={{ marginTop: "16px" }}>
          <p style={{ fontWeight: "bold", marginBottom: "8px" }}>Try one of these options:</p>
          <ul style={{ textAlign: "left", paddingLeft: "20px" }}>
            <li>Use <strong>Google Chrome</strong> or <strong>Microsoft Edge</strong></li>
            <li>Use <strong>Firefox</strong> (version 79+)</li>
            <li>Use the CCP4i2 desktop application (Electron)</li>
          </ul>
        </div>
      </div>
    </div>
  );
};

/**
 * Error boundary for catching Moorhen rendering errors.
 */
interface MoorhenErrorBoundaryProps {
  children: React.ReactNode;
  fallback?: React.ReactNode;
}

interface MoorhenErrorBoundaryState {
  hasError: boolean;
  error: Error | null;
}

export class MoorhenErrorBoundary extends React.Component<
  MoorhenErrorBoundaryProps,
  MoorhenErrorBoundaryState
> {
  constructor(props: MoorhenErrorBoundaryProps) {
    super(props);
    this.state = { hasError: false, error: null };
  }

  static getDerivedStateFromError(error: Error): MoorhenErrorBoundaryState {
    return { hasError: true, error };
  }

  componentDidCatch(error: Error, errorInfo: React.ErrorInfo) {
    console.error("[MoorhenErrorBoundary] Caught error:", error, errorInfo);
  }

  render() {
    if (this.state.hasError) {
      if (this.props.fallback) {
        return this.props.fallback;
      }
      return (
        <MoorhenFallback
          reason="runtime_error"
          error={this.state.error || undefined}
        />
      );
    }

    return this.props.children;
  }
}

/**
 * Hook for checking Moorhen capabilities with client-side only execution.
 *
 * Returns undefined during SSR, then the actual capabilities after hydration.
 */
export function useMoorhenCapabilities(): MoorhenCapabilities | undefined {
  const [capabilities, setCapabilities] = useState<MoorhenCapabilities | undefined>(
    undefined
  );

  useEffect(() => {
    setCapabilities(checkMoorhenCapabilities());
  }, []);

  return capabilities;
}
