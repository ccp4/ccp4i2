/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
  cootModuleError?: Error | null;
  setCootModuleError?: (error: Error | null) => void;
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
  cootModuleError: null,
  setCootModuleError: () => {},
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
