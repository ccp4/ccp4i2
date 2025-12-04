import React, { createContext, useContext, useState } from "react";

interface JobTabContextType {
  jobTabValue: number;
  setJobTabValue: (value: number) => void;
}

const JobTabContext = createContext<JobTabContextType | undefined>(undefined);

export const JobTabProvider: React.FC<{ children: React.ReactNode }> = ({
  children,
}) => {
  const [jobTabValue, setJobTabValue] = useState<number>(0);

  return (
    <JobTabContext.Provider value={{ jobTabValue, setJobTabValue }}>
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
