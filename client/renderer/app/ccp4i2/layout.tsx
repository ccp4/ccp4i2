"use client";
import { PropsWithChildren } from "react";

/**
 * /ccp4i2 segment layout — passthrough.
 *
 * The real provider stack (RequireAuth, CCP4i2App, FindInPageProvider,
 * RecentlyStartedJobsProvider, DeleteDialogProvider, PopcornProvider,
 * AuthErrorHandler) lives one level deeper at app/ccp4i2/(authed)/layout.tsx.
 * Everything except /ccp4i2/config sits inside that route group.
 *
 * Config stays here, outside the (authed) group, because it must remain
 * reachable when Django isn't running yet (initial setup, server
 * reconfiguration). RequireAuth on a non-running backend just spins;
 * AuthErrorHandler on a non-running backend would surface "session
 * expired" snackbars in response to connection failures.
 *
 * Kept as an explicit file (rather than deleted) so this rationale lives
 * with the layout boundary it defends.
 */
export default function CCP4i2Layout(props: PropsWithChildren) {
  return <>{props.children}</>;
}
