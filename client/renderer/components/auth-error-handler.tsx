"use client";

import { useEffect, useRef } from "react";
import { usePopcorn } from "../providers/popcorn-provider";
import { AUTH_ERROR_EVENT, AuthErrorDetail } from "../api-fetch";
import { logout } from "@ccp4/ccp4i2-api";
import { reauthenticate } from "../utils/reauth";

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
        // Show snackbar with an explicit "Sign in" action so the user
        // can recover deliberately. No auto-logout: yanking the user
        // away after 2s gives them no time to read the message and no
        // agency to wait/cancel. The snackbar stays put (popcorn's
        // action-snackbars don't auto-hide) until they click.
        //
        // The action signs the user back IN. It used to call logout(),
        // which tore down the AAD session that would have made the return
        // trip silent and then sent them through credentials and MFA to
        // land on "/" rather than the page they were reading. Signing out
        // belongs on the Sign-out menu item, not on a recovery prompt.
        // Falling back to logout() keeps the button honest where no re-auth
        // is possible -- there is nothing else it could usefully do.
        setMessage(message, "error", {
          label: "Sign in",
          onClick: async () => {
            switch (await reauthenticate()) {
              case "started":
                // The page is on its way to AAD; say nothing over it.
                return;
              case "in-progress":
                setMessage("Signing in — this may take a moment.", "info");
                return;
              case "failed":
                // It did not start, but the session is exactly as it was, so
                // signing out would cost the user more than it saves. Offer
                // it as their choice rather than doing it to them.
                setMessage("Could not start sign-in.", "error", {
                  label: "Sign out",
                  onClick: () => logout(),
                });
                return;
              case "unavailable":
                // Nothing can re-authenticate here -- the desktop session,
                // where logout is itself a no-op. Keeps the button honest.
                logout();
                return;
            }
          },
        });
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
