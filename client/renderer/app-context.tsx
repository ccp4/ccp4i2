import React, { createContext, useContext } from "react";
import { File, Job } from "./types/models";

interface CCP4i2Context {
  projectId?: number | null;
  setProjectId?: (projectId: number) => void | null;
  jobId?: number | null;
  setJobId?: (jobId: number | null) => void | null;
  jobPanelSize?: number;
  setJobPanelSize?: (size: number) => void;
  cootModule?: any | null;
  setCootModule?: (module: any | null) => void;
  rdkitModule?: any | null;
  setRdkitModule?: (module: any | null) => void;
  devMode: boolean;
  setDevMode: (devMode: boolean) => void;
  activeDragItem: Job | File | null;
  setActiveDragItem: (item: Job | File | null) => void;
}
export const CCP4i2Context = createContext<CCP4i2Context>({
  projectId: null,
  setProjectId: () => {},
  jobId: null,
  setJobId: () => {},
  jobPanelSize: 7,
  setJobPanelSize: () => {},
  cootModule: null,
  setCootModule: () => {},
  rdkitModule: null,
  setRdkitModule: () => {},
  devMode: true,
  setDevMode: () => {},
  activeDragItem: null,
  setActiveDragItem: () => {},
});

export const useCCP4i2Window = () => {
  const context = useContext(CCP4i2Context);
  if (!context) {
    throw new Error("useCCP4i2Window must be used within a CCP4i2Provider");
  }
  return context;
};
