"use client";

import React, { PropsWithChildren, useEffect } from "react";

export function NavigationShortcutsProvider({ children }: PropsWithChildren) {
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      // Mac: Cmd+[ or Cmd+]
      if (e.metaKey && e.key === "[") {
        e.preventDefault();
        window.history.back();
      }
      if (e.metaKey && e.key === "]") {
        e.preventDefault();
        window.history.forward();
      }
      // Windows/Linux: Alt+Left or Alt+Right
      if (e.altKey && e.key === "ArrowLeft") {
        e.preventDefault();
        window.history.back();
      }
      if (e.altKey && e.key === "ArrowRight") {
        e.preventDefault();
        window.history.forward();
      }
    };
    window.addEventListener("keydown", handleKeyDown);
    return () => window.removeEventListener("keydown", handleKeyDown);
  }, []);

  return <>{children}</>;
}
