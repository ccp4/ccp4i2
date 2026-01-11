"use client";
import { PropsWithChildren } from "react";
import { DeleteDialogProvider } from "../../providers/delete-dialog";
import { RecentlyStartedJobsProvider } from "../../providers/recently-started-jobs-context";
import { CCP4i2App } from "../../providers/ccp4i2-app";
import RequireAuth from "../../components/require-auth";

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

/**
 * CCP4i2 app layout - wraps all ccp4i2 pages with:
 * - Authentication (when required) - AuthProvider is in root layout
 * - Job tracking context
 * - Delete dialog provider
 * - CCP4i2 app shell (header, navigation)
 */
export default function CCP4i2Layout(props: PropsWithChildren) {
  if (REQUIRE_AUTH) {
    return (
      <RecentlyStartedJobsProvider>
        <DeleteDialogProvider>
          <RequireAuth>
            <CCP4i2App>{props.children}</CCP4i2App>
          </RequireAuth>
        </DeleteDialogProvider>
      </RecentlyStartedJobsProvider>
    );
  }

  return (
    <RecentlyStartedJobsProvider>
      <DeleteDialogProvider>
        <CCP4i2App>{props.children}</CCP4i2App>
      </DeleteDialogProvider>
    </RecentlyStartedJobsProvider>
  );
}
