"use client";
import React, { PropsWithChildren, useEffect, useState, useCallback } from "react";
import { CCP4i2Context } from "../app-context";
import { CssBaseline } from "@mui/material";
import { File, Job } from "../types/models";
import { PopcornProvider, usePopcorn } from "./popcorn-provider";
import { RunCheckProvider } from "./run-check-provider";
import { useStalledJobWarnings } from "./recently-started-jobs-context";

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
  const [rdkitModule, setRdkitModule] = useState<any | null>(null);
  const [jobPanelSize, setJobPanelSize] = useState<number>(70);
  const [devMode, setDevMode] = useState<boolean>(true);
  const [activeDragItem, setActiveDragItem] = useState<Job | File | null>(null);

  return (
    <PopcornProvider>
      <CCP4i2Context.Provider
        value={{
          projectId,
          setProjectId,
          jobId,
          setJobId,
          cootModule,
          setCootModule,
          rdkitModule,
          setRdkitModule,
          jobPanelSize,
          setJobPanelSize,
          devMode,
          setDevMode,
          activeDragItem,
          setActiveDragItem,
        }}
      >
        <CssBaseline />
        <RunCheckProvider>
          <StalledJobWarningsHandler>
            {props.children}
          </StalledJobWarningsHandler>
        </RunCheckProvider>
      </CCP4i2Context.Provider>
    </PopcornProvider>
  );
};
