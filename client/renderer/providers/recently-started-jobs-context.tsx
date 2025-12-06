"use client";

import React, {
  createContext,
  useContext,
  useState,
  useCallback,
  useRef,
  useEffect,
} from "react";

// ============================================================================
// Types
// ============================================================================

interface StartedJobInfo {
  jobId: number;
  projectId: number;
  jobNumber: string;
  jobTitle: string;
  startTime: number;
  /** Last observed status from polling (updated by components) */
  lastObservedStatus?: number;
  /** Whether we've already warned about this job being stalled */
  warnedAboutStall: boolean;
}

interface RecentlyStartedJobsContextValue {
  /**
   * Mark a job as recently started. This triggers fast polling for the grace period
   * even if the database status hasn't caught up yet.
   */
  markJobAsStarting: (
    jobId: number,
    projectId: number,
    jobNumber: string,
    jobTitle: string
  ) => void;

  /**
   * Update the observed status of a tracked job.
   * Called by polling hooks when they fetch fresh job data.
   */
  updateJobObservedStatus: (jobId: number, status: number) => void;

  /**
   * Check if a job was recently started and is still within the grace period.
   */
  isJobRecentlyStarted: (jobId: number) => boolean;

  /**
   * Check if ANY job in the given project was recently started.
   * Useful for project-level polling decisions.
   */
  hasRecentlyStartedJobsInProject: (projectId: number) => boolean;

  /**
   * Get any stalled job warnings that need to be shown.
   * Returns and clears the pending warnings.
   */
  consumeStalledJobWarnings: () => StalledJobWarning[];
}

interface StalledJobWarning {
  jobId: number;
  jobNumber: string;
  jobTitle: string;
  message: string;
}

// ============================================================================
// Constants
// ============================================================================

/** Grace period after a job is started before we stop treating it as "active" */
const GRACE_PERIOD_MS = 20000; // 20 seconds

/** How often to clean up expired entries and check for stalls */
const CLEANUP_INTERVAL_MS = 5000;

/** Active job statuses: Queued (2), Running (3), Running remotely (7) */
const ACTIVE_JOB_STATUSES = [2, 3, 7];

// ============================================================================
// Context
// ============================================================================

const RecentlyStartedJobsContext =
  createContext<RecentlyStartedJobsContextValue | null>(null);

// ============================================================================
// Provider
// ============================================================================

interface RecentlyStartedJobsProviderProps {
  children: React.ReactNode;
}

export const RecentlyStartedJobsProvider: React.FC<
  RecentlyStartedJobsProviderProps
> = ({ children }) => {
  // Map of jobId -> job info
  const [startedJobs, setStartedJobs] = useState<Map<number, StartedJobInfo>>(
    new Map()
  );

  // Pending stalled job warnings
  const [stalledWarnings, setStalledWarnings] = useState<StalledJobWarning[]>(
    []
  );

  // Mark a job as recently started
  const markJobAsStarting = useCallback(
    (
      jobId: number,
      projectId: number,
      jobNumber: string,
      jobTitle: string
    ) => {
      console.log(
        `[RecentlyStartedJobs] Marking job ${jobNumber} (${jobId}) as starting`
      );
      setStartedJobs((prev) => {
        const next = new Map(prev);
        next.set(jobId, {
          jobId,
          projectId,
          jobNumber,
          jobTitle,
          startTime: Date.now(),
          lastObservedStatus: undefined,
          warnedAboutStall: false,
        });
        return next;
      });
    },
    []
  );

  // Update observed status for a tracked job
  const updateJobObservedStatus = useCallback(
    (jobId: number, status: number) => {
      setStartedJobs((prev) => {
        const info = prev.get(jobId);
        if (!info) return prev;

        // If status changed, update it
        if (info.lastObservedStatus !== status) {
          const next = new Map(prev);
          next.set(jobId, { ...info, lastObservedStatus: status });
          return next;
        }
        return prev;
      });
    },
    []
  );

  // Check if a specific job was recently started
  const isJobRecentlyStarted = useCallback(
    (jobId: number): boolean => {
      const info = startedJobs.get(jobId);
      if (!info) return false;
      return Date.now() - info.startTime < GRACE_PERIOD_MS;
    },
    [startedJobs]
  );

  // Check if any job in a project was recently started
  const hasRecentlyStartedJobsInProject = useCallback(
    (projectId: number): boolean => {
      const now = Date.now();
      for (const [, info] of startedJobs.entries()) {
        if (now - info.startTime >= GRACE_PERIOD_MS) continue;
        if (info.projectId === projectId) return true;
      }
      return false;
    },
    [startedJobs]
  );

  // Consume and clear pending stalled job warnings
  const consumeStalledJobWarnings = useCallback((): StalledJobWarning[] => {
    const warnings = [...stalledWarnings];
    if (warnings.length > 0) {
      setStalledWarnings([]);
    }
    return warnings;
  }, [stalledWarnings]);

  // Periodic cleanup and stall detection
  useEffect(() => {
    const cleanup = () => {
      const now = Date.now();
      const newWarnings: StalledJobWarning[] = [];

      setStartedJobs((prev) => {
        let hasChanges = false;
        const next = new Map<number, StartedJobInfo>();

        for (const [jobId, info] of prev.entries()) {
          const elapsed = now - info.startTime;

          if (elapsed < GRACE_PERIOD_MS) {
            // Still within grace period, keep tracking
            next.set(jobId, info);
          } else {
            // Grace period expired
            hasChanges = true;

            // Check if job became active
            const becameActive =
              info.lastObservedStatus !== undefined &&
              ACTIVE_JOB_STATUSES.includes(info.lastObservedStatus);

            if (!becameActive && !info.warnedAboutStall) {
              // Job didn't become active within grace period - likely stalled
              console.warn(
                `[RecentlyStartedJobs] Job ${info.jobNumber} (${jobId}) appears stalled - ` +
                  `last observed status: ${info.lastObservedStatus ?? "unknown"}`
              );
              newWarnings.push({
                jobId: info.jobId,
                jobNumber: info.jobNumber,
                jobTitle: info.jobTitle,
                message: `Job ${info.jobNumber} "${info.jobTitle}" may not have started - check the server logs`,
              });
            } else if (becameActive) {
              console.log(
                `[RecentlyStartedJobs] Job ${info.jobNumber} successfully started (status: ${info.lastObservedStatus})`
              );
            }
            // Don't add to next map - remove from tracking
          }
        }

        if (!hasChanges && newWarnings.length === 0) return prev;
        return next;
      });

      // Add any new warnings
      if (newWarnings.length > 0) {
        setStalledWarnings((prev) => [...prev, ...newWarnings]);
      }
    };

    const intervalId = setInterval(cleanup, CLEANUP_INTERVAL_MS);
    return () => clearInterval(intervalId);
  }, []);

  const value: RecentlyStartedJobsContextValue = {
    markJobAsStarting,
    updateJobObservedStatus,
    isJobRecentlyStarted,
    hasRecentlyStartedJobsInProject,
    consumeStalledJobWarnings,
  };

  return (
    <RecentlyStartedJobsContext.Provider value={value}>
      {children}
    </RecentlyStartedJobsContext.Provider>
  );
};

// ============================================================================
// Hooks
// ============================================================================

export const useRecentlyStartedJobs = (): RecentlyStartedJobsContextValue => {
  const context = useContext(RecentlyStartedJobsContext);
  if (!context) {
    throw new Error(
      "useRecentlyStartedJobs must be used within a RecentlyStartedJobsProvider"
    );
  }
  return context;
};

/**
 * Convenience hook that combines the context check with actual job status.
 * Returns true if the job is either:
 * 1. In an active database status (2=Queued, 3=Running, 7=Running remotely)
 * 2. Was recently started and is still within the grace period
 *
 * Also updates the observed status in the context for stall detection.
 */
export const useIsJobEffectivelyActive = (
  jobId: number | undefined | null,
  jobStatus: number | undefined
): boolean => {
  const context = useContext(RecentlyStartedJobsContext);

  // Update observed status when it changes
  useEffect(() => {
    if (context && jobId && jobStatus !== undefined) {
      context.updateJobObservedStatus(jobId, jobStatus);
    }
  }, [context, jobId, jobStatus]);

  // If no context (provider not mounted), fall back to status-only check
  if (!context) {
    return jobStatus !== undefined && ACTIVE_JOB_STATUSES.includes(jobStatus);
  }

  const { isJobRecentlyStarted } = context;

  // Check database status
  const isStatusActive =
    jobStatus !== undefined && ACTIVE_JOB_STATUSES.includes(jobStatus);

  // Check grace period
  const isRecentlyStarted = jobId ? isJobRecentlyStarted(jobId) : false;

  return isStatusActive || isRecentlyStarted;
};

/**
 * Hook to display stalled job warnings via Popcorn.
 * Should be used once at a high level in the app.
 */
export const useStalledJobWarnings = (
  showWarning: (message: string) => void
) => {
  const context = useContext(RecentlyStartedJobsContext);

  useEffect(() => {
    if (!context) return;

    const checkWarnings = () => {
      const warnings = context.consumeStalledJobWarnings();
      for (const warning of warnings) {
        showWarning(warning.message);
      }
    };

    // Check for warnings periodically
    const intervalId = setInterval(checkWarnings, 1000);
    return () => clearInterval(intervalId);
  }, [context, showWarning]);
};
