"use client";
import { PropsWithChildren } from "react";
import { DeleteDialogProvider } from "@/providers/delete-dialog";
import { FindInPageProvider } from "@/providers/find-in-page-provider";
import { RecentlyStartedJobsProvider } from "@/providers/recently-started-jobs-context";
import { CCP4i2App } from "@/providers/ccp4i2-app";
import { PopcornProvider } from "@/providers/popcorn-provider";
import { AuthErrorHandler } from "@/components/auth-error-handler";
import RequireAuth from "@/components/require-auth";

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

/**
 * (authed) route-group layout — wraps every real ccp4i2 route (project,
 * job, campaigns, moorhen-page, import/new-project, graph-viewer, plus
 * the /ccp4i2 landing page).
 *
 * /ccp4i2/config/* deliberately sits one level up, outside this group:
 * config can be reached when the Django server hasn't started yet (or
 * is being reconfigured), and we don't want AuthErrorHandler firing a
 * "Your session has expired" snackbar + auto-logout in response to
 * connection failures from a backend the user is mid-way through
 * setting up.
 *
 * PopcornProvider + AuthErrorHandler live here (not inside CCP4i2App)
 * so every (authed) route — and any modal/dialog opened from one — sees
 * the snackbar surface and the 401-→-logout handler. Previously they
 * were nested inside CCP4i2App, which meant components mounted outside
 * the app shell (e.g. landing-page dialogs) had no popcorn context and
 * 401s went to the void.
 */
export default function AuthedLayout(props: PropsWithChildren) {
  return (
    <PopcornProvider>
      <AuthErrorHandler />
      <FindInPageProvider>
        <RecentlyStartedJobsProvider>
          <DeleteDialogProvider>
            {REQUIRE_AUTH ? (
              <RequireAuth>
                <CCP4i2App>{props.children}</CCP4i2App>
              </RequireAuth>
            ) : (
              <CCP4i2App>{props.children}</CCP4i2App>
            )}
          </DeleteDialogProvider>
        </RecentlyStartedJobsProvider>
      </FindInPageProvider>
    </PopcornProvider>
  );
}
