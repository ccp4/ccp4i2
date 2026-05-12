"use client";
import React, { PropsWithChildren, useState, useCallback, useMemo } from "react";
import { CCP4i2Context } from "../app-context";
import { CssBaseline } from "@mui/material";
import { File, Job } from "../types/models";
import { usePopcorn } from "./popcorn-provider";
import { RunCheckProvider } from "./run-check-provider";
import { useStalledJobWarnings } from "./recently-started-jobs-context";
// PopcornProvider + AuthErrorHandler are mounted by the
// app/ccp4i2/(authed)/layout.tsx route group, so all real ccp4i2 routes
// share one snackbar surface and one 401 handler — including modals and
// dialogs that mount outside the CCP4i2App shell. /ccp4i2/config lives
// outside (authed) and so does NOT get them, by design.

/**
 * Component to display stalled job warnings via Popcorn.
 * Must be rendered inside both PopcornProvider and RecentlyStartedJobsProvider.
 */
const StalledJobWarningsHandler: React.FC<PropsWithChildren> = ({ children }) => {
  const { setMessage } = usePopcorn();

  // Memoized callback to display warnings - uses setMessage directly
  const showWarning = useCallback(
    (message: string) => {
      setMessage(message, "warning");
    },
    [setMessage]
  );

  // Hook that periodically checks for stalled job warnings and displays them
  useStalledJobWarnings(showWarning);

  return <>{children}</>;
};

export const CCP4i2App = (props: PropsWithChildren) => {
  const [projectId, setProjectId] = useState<number | null>(null);
  const [jobId, setJobId] = useState<number | null>(null);
  const [cootModule, setCootModule] = useState<any | null>(null);
  const [cootModuleError, setCootModuleError] = useState<Error | null>(null);
  const [devMode, setDevMode] = useState<boolean>(true);
  const [activeDragItem, setActiveDragItem] = useState<Job | File | null>(null);

  // Memoize context value to prevent unnecessary re-renders of consumers
  // Note: RDKit module is provided separately via RDKitProvider (useRDKit hook)
  const contextValue = useMemo(
    () => ({
      projectId,
      setProjectId,
      jobId,
      setJobId,
      cootModule,
      setCootModule,
      cootModuleError,
      setCootModuleError,
      devMode,
      setDevMode,
      activeDragItem,
      setActiveDragItem,
    }),
    [projectId, jobId, cootModule, cootModuleError, devMode, activeDragItem]
  );

  return (
    <CCP4i2Context.Provider value={contextValue}>
      <CssBaseline />
      <RunCheckProvider>
        <StalledJobWarningsHandler>
          {props.children}
        </StalledJobWarningsHandler>
      </RunCheckProvider>
    </CCP4i2Context.Provider>
  );
};
