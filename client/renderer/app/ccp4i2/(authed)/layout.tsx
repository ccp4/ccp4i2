"use client";
import { PropsWithChildren } from "react";
import { usePathname } from "next/navigation";
import { Stack } from "@mui/material";
import { DeleteDialogProvider } from "@/providers/delete-dialog";
import { ImportProvenanceProvider } from "@/providers/import-provenance-provider";
import { FindInPageProvider } from "@/providers/find-in-page-provider";
import { RecentlyStartedJobsProvider } from "@/providers/recently-started-jobs-context";
import { RunningProcessesProvider } from "@/providers/running-processes";
import { CCP4i2App } from "@/providers/ccp4i2-app";
import { PopcornProvider } from "@/providers/popcorn-provider";
import { TopBarProvider } from "@/providers/top-bar-context";
import { AuthErrorHandler } from "@/components/auth-error-handler";
import CCP4i2AppBar from "@/components/ccp4i2-app-bar";
import RequireAuth from "@/components/require-auth";

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

/**
 * Routes that bring their own chrome. Both open as separate OS windows with
 * no route back into i2, so the shared bar's menus would only offer to
 * navigate them somewhere they cannot return from; Moorhen also needs the
 * vertical space.
 */
const OWN_CHROME = ["/ccp4i2/moorhen-page", "/ccp4i2/graph-viewer"];

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
 *
 * The app bar is mounted here, once, so that every route gets the same one
 * (pages contribute their name through useTopBar). RunningProcessesProvider
 * has to sit above it, because the Utilities menu in the bar opens that
 * dialog and the dialog is about jobs across all projects, not this one.
 */
export default function AuthedLayout(props: PropsWithChildren) {
  const pathname = usePathname() ?? "";
  const showAppBar = !OWN_CHROME.some((prefix) => pathname.startsWith(prefix));

  const shell = (
    <RunningProcessesProvider>
      <TopBarProvider>
        {showAppBar ? (
          <Stack sx={{ height: "100svh", width: "100%", overflow: "hidden" }}>
            <CCP4i2AppBar />
            <Stack sx={{ flex: 1, minHeight: 0 }}>{props.children}</Stack>
          </Stack>
        ) : (
          props.children
        )}
      </TopBarProvider>
    </RunningProcessesProvider>
  );

  return (
    <PopcornProvider>
      <AuthErrorHandler />
      <FindInPageProvider>
        <RecentlyStartedJobsProvider>
          <DeleteDialogProvider>
            <ImportProvenanceProvider>
              {REQUIRE_AUTH ? (
                <RequireAuth>
                  <CCP4i2App>{shell}</CCP4i2App>
                </RequireAuth>
              ) : (
                <CCP4i2App>{shell}</CCP4i2App>
              )}
            </ImportProvenanceProvider>
          </DeleteDialogProvider>
        </RecentlyStartedJobsProvider>
      </FindInPageProvider>
    </PopcornProvider>
  );
}
