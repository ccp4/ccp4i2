import React, {
  createContext,
  PropsWithChildren,
  useContext,
  useMemo,
} from "react";
import type { Job } from "../../../types/models";
import { useJob } from "../../../utils";

interface TaskInterfaceContextValue {
  job: Job;
  useTaskItem: ReturnType<typeof useJob>["useTaskItem"];
}

const TaskInterfaceContext = createContext<TaskInterfaceContextValue | null>(
  null
);

export const TaskInterfaceProvider: React.FC<
  PropsWithChildren<{ job: Job }>
> = ({ job, children }) => {
  const { useTaskItem } = useJob(job.id);
  const value = useMemo(
    () => ({ job, useTaskItem }),
    [job, useTaskItem]
  );
  return (
    <TaskInterfaceContext.Provider value={value}>
      {children}
    </TaskInterfaceContext.Provider>
  );
};

export function useTaskInterface(): TaskInterfaceContextValue {
  const ctx = useContext(TaskInterfaceContext);
  if (!ctx) {
    throw new Error(
      "useTaskInterface must be used inside <TaskInterfaceProvider>"
    );
  }
  return ctx;
}
