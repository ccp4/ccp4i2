/*
 * Copyright (C) 2025 Newcastle University
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
