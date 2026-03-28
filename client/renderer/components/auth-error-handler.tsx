/*
 * Copyright (C) 2026 Newcastle University
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

import { useEffect, useRef } from "react";
import { usePopcorn } from "../providers/popcorn-provider";
import { AUTH_ERROR_EVENT, AuthErrorDetail } from "../api-fetch";
import { logout } from "../utils/auth-token";

/**
 * Listens for 401/403 auth errors from the API layer and surfaces them
 * to the user via snackbar notifications.
 *
 * - 401: Session expired — shows error snackbar, then triggers re-auth
 * - 403: Forbidden — shows warning snackbar (no redirect, user may need
 *   to contact an admin)
 *
 * Debounces rapid-fire errors (e.g. many SWR hooks failing at once)
 * so the user sees a single notification rather than a flood.
 *
 * Must be rendered inside PopcornProvider.
 */
export const AuthErrorHandler: React.FC = () => {
  const { setMessage } = usePopcorn();
  const lastNotifiedAt = useRef(0);

  useEffect(() => {
    const DEBOUNCE_MS = 5000;

    const handleAuthError = (event: Event) => {
      const { status, message } = (event as CustomEvent<AuthErrorDetail>).detail;

      // Debounce: only show one notification per window
      const now = Date.now();
      if (now - lastNotifiedAt.current < DEBOUNCE_MS) return;
      lastNotifiedAt.current = now;

      if (status === 401) {
        setMessage(message, "error");
        // Give the user a moment to read the snackbar before redirecting
        setTimeout(() => logout(), 2000);
      } else {
        // 403 — don't redirect, just inform
        setMessage(message, "warning");
      }
    };

    window.addEventListener(AUTH_ERROR_EVENT, handleAuthError);
    return () => window.removeEventListener(AUTH_ERROR_EVENT, handleAuthError);
  }, [setMessage]);

  return null;
};
