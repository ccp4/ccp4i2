"use client";

import { useEffect, useRef } from "react";
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
 * @param route The route to save
 * @param isDeliberateHomeNavigation If true, navigating to "/" clears saved route
 */
export function saveTeamsRoute(route: string, isDeliberateHomeNavigation: boolean = false): void {
  if (typeof window === "undefined") return;
  if (!isRunningInIframe()) return;

  // When deliberately navigating to root (not fresh load), clear saved route
  if (route === "/" && isDeliberateHomeNavigation) {
    clearSavedTeamsRoute();
    return;
  }

  // Don't save root, auth pages, or API routes
  if (route === "/" || route.startsWith("/auth/") || route.startsWith("/api/")) {
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
 *
 * Tracks whether user has visited a non-root route to distinguish:
 * - Fresh Teams reload to "/" -> should restore saved route
 * - Deliberate navigation to "/" -> should clear and stay
 */
export function TeamsRoutePersistence() {
  const pathname = usePathname();
  const searchParams = useSearchParams();
  // Track if we've visited a non-root route (meaning navigation to "/" is deliberate)
  const hasVisitedNonRootRoute = useRef(false);

  useEffect(() => {
    // Build full path including search params
    const search = searchParams.toString();
    const fullPath = search ? `${pathname}?${search}` : pathname;

    // Check if this is deliberate navigation to home (we've been to another route first)
    const isDeliberateHomeNavigation = pathname === "/" && hasVisitedNonRootRoute.current;

    // Track if we've visited a non-root route
    if (pathname !== "/" && !pathname.startsWith("/auth/")) {
      hasVisitedNonRootRoute.current = true;
    }

    saveTeamsRoute(fullPath, isDeliberateHomeNavigation);
  }, [pathname, searchParams]);

  // This component renders nothing
  return null;
}
