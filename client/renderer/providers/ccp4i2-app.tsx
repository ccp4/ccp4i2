"use client";
import { PropsWithChildren, useEffect, useState } from "react";
import { CCP4i2Context } from "../app-context";
import { CssBaseline } from "@mui/material";
import { File, Job } from "../types/models";
import { PopcornProvider } from "./popcorn-provider";
import { RunCheckProvider } from "./run-check-provider";

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
        <RunCheckProvider>{props.children}</RunCheckProvider>
      </CCP4i2Context.Provider>
    </PopcornProvider>
  );
};
