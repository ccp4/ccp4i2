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
