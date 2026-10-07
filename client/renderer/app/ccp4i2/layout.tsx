"use client";
import { PropsWithChildren } from "react";
import { SWRConfig } from "swr";
import { ApiError } from "@ccp4/ccp4i2-api";
import { NavigationShortcutsProvider } from "@/providers/navigation-shortcuts-provider";

/**
 * /ccp4i2 segment layout.
 *
 * The real provider stack (RequireAuth, CCP4i2App, FindInPageProvider,
 * RecentlyStartedJobsProvider, DeleteDialogProvider, PopcornProvider,
 * AuthErrorHandler) lives one level deeper at app/ccp4i2/(authed)/layout.tsx.
 * Everything except /ccp4i2/config sits inside that route group.
 *
 * Config stays there, outside the (authed) group, because it must remain
 * reachable when Django isn't running yet (initial setup, server
 * reconfiguration). RequireAuth on a non-running backend just spins;
 * AuthErrorHandler on a non-running backend would surface "session
 * expired" snackbars in response to connection failures.
 *
 * The one thing that belongs at this level is what every /ccp4i2 route needs
 * regardless of auth: the back/forward keyboard shortcuts. Mounting them here
 * once (rather than per page) is what makes them work on every route — the
 * projects list, campaigns and graph viewer used to be missed. The SWR
 * defaults belong here for the same reason.
 */
// SWR retries failed fetches indefinitely by default. A 4xx (a deleted job,
// a missing file) will fail the same way every time, so only retry the rest.
const shouldRetryOnError = (error: Error) =>
  !(error instanceof ApiError && error.status >= 400 && error.status < 500);

export default function CCP4i2Layout(props: PropsWithChildren) {
  return (
    <SWRConfig value={{ shouldRetryOnError }}>
      <NavigationShortcutsProvider>{props.children}</NavigationShortcutsProvider>
    </SWRConfig>
  );
}
