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
import React, { createContext, useContext, useState, useMemo } from "react";

interface JobTabContextType {
  jobTabValue: number;
  setJobTabValue: (value: number) => void;
}

const JobTabContext = createContext<JobTabContextType | undefined>(undefined);

export const JobTabProvider: React.FC<{ children: React.ReactNode }> = ({
  children,
}) => {
  const [jobTabValue, setJobTabValue] = useState<number>(0);

  const contextValue = useMemo(
    () => ({ jobTabValue, setJobTabValue }),
    [jobTabValue]
  );

  return (
    <JobTabContext.Provider value={contextValue}>
      {children}
    </JobTabContext.Provider>
  );
};

export function useJobTab() {
  const context = useContext(JobTabContext);
  if (!context) {
    throw new Error("useJobTab must be used within a JobTabProvider");
  }
  return context;
}
