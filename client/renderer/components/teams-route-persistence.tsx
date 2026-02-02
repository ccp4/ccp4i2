"use client";

import { useEffect } from "react";
import { usePathname, useSearchParams } from "next/navigation";

const TEAMS_ROUTE_KEY = "ccp4i2-teams-last-route";

/**
 * Check if running in an iframe (Teams context).
 */
function isRunningInIframe(): boolean {
  if (typeof window === "undefined") return false;
  try {
    return window.self !== window.top;
  } catch {
    return true;
  }
}

/**
 * Save the last visited route to sessionStorage.
 * Only saves non-root routes and only in Teams context.
 *
 * When navigating to "/" (app-selector), clears the saved route.
 * This allows deliberate navigation to the app-selector without bounce-back.
 */
export function saveTeamsRoute(route: string): void {
  if (typeof window === "undefined") return;
  if (!isRunningInIframe()) return;

  // When navigating to root, clear saved route (deliberate navigation to app-selector)
  if (route === "/") {
    clearSavedTeamsRoute();
    return;
  }

  // Don't save auth pages or API routes
  if (route.startsWith("/auth/") || route.startsWith("/api/")) {
    return;
  }

  try {
    sessionStorage.setItem(TEAMS_ROUTE_KEY, route);
  } catch (e) {
    // Ignore storage errors
  }
}

/**
 * Get the saved route from sessionStorage.
 * Returns null if not in Teams context or no route saved.
 */
export function getSavedTeamsRoute(): string | null {
  if (typeof window === "undefined") return null;
  if (!isRunningInIframe()) return null;

  try {
    return sessionStorage.getItem(TEAMS_ROUTE_KEY);
  } catch {
    return null;
  }
}

/**
 * Clear the saved route (called after restoring).
 */
export function clearSavedTeamsRoute(): void {
  if (typeof window === "undefined") return;

  try {
    sessionStorage.removeItem(TEAMS_ROUTE_KEY);
  } catch {
    // Ignore
  }
}

/**
 * Component that tracks route changes and saves them for Teams context.
 * Add this to the root layout to enable route persistence.
 */
export function TeamsRoutePersistence() {
  const pathname = usePathname();
  const searchParams = useSearchParams();

  useEffect(() => {
    // Build full path including search params
    const search = searchParams.toString();
    const fullPath = search ? `${pathname}?${search}` : pathname;

    saveTeamsRoute(fullPath);
  }, [pathname, searchParams]);

  // This component renders nothing
  return null;
}
