"use client";

import React, { PropsWithChildren, useEffect } from "react";

/**
 * Browser-style back/forward keyboard shortcuts:
 *   Mac:           Cmd+[  /  Cmd+]
 *   Windows/Linux: Alt+Left  /  Alt+Right
 *
 * Mounted once, in app/ccp4i2/layout.tsx, so it covers every /ccp4i2 route.
 * The listener is reference-counted at module level, so if a second copy is
 * ever mounted inside the first (or React StrictMode replays the effect) there
 * is still exactly one keydown handler — a nested pair must not call
 * history.back() twice per keystroke.
 */
function handleKeyDown(e: KeyboardEvent) {
  if (e.metaKey && e.key === "[") {
    e.preventDefault();
    window.history.back();
  } else if (e.metaKey && e.key === "]") {
    e.preventDefault();
    window.history.forward();
  } else if (e.altKey && e.key === "ArrowLeft") {
    e.preventDefault();
    window.history.back();
  } else if (e.altKey && e.key === "ArrowRight") {
    e.preventDefault();
    window.history.forward();
  }
}

let mountedProviders = 0;

export function NavigationShortcutsProvider({ children }: PropsWithChildren) {
  useEffect(() => {
    if (mountedProviders++ === 0) {
      window.addEventListener("keydown", handleKeyDown);
    }
    return () => {
      if (--mountedProviders === 0) {
        window.removeEventListener("keydown", handleKeyDown);
      }
    };
  }, []);

  return <>{children}</>;
}
