import React, {
  createContext,
  useCallback,
  useContext,
  useState,
  useMemo,
  ReactNode,
  useEffect,
  useRef,
} from "react";
import {
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
} from "@mui/material";
import { Button } from "@mui/material";
import { useCCP4i2Window } from "../app-context";
import { useApi } from "../api";
import { useJob, useProject } from "../utils";
import { File as DjangoFile, Project } from "../types/models";
import { InlineTaskModal } from "../components/task/task-elements/inline-task-modal";

/**
 * An error report keyed by parameter objectPath.
 * Each entry contains a maxSeverity (number) and an array of messages (strings).
 *
 * Severity mapping (from server XML):
 *   1 = WARNING / UNDEFINED_ERROR — shown in dialog, does not block
 *   2 = ERROR — shown in dialog, blocks Confirm button
 */
export interface CCP4i2ErrorReport {
  [key: string]: {
    maxSeverity: number;
    messages: string[];
  };
}

interface RunCheckContextType {
  runTaskRequested: number | null;
  setRunTaskRequested: (taskId: number | null) => void;
  confirmTaskRun: (taskId: number) => Promise<boolean>;
}

export const RunCheckContext = createContext<RunCheckContextType>({
  runTaskRequested: null,
  setRunTaskRequested: () => {},
  confirmTaskRun: () => Promise.resolve(false),
});

interface RunCheckProviderProps {
  children: ReactNode;
}

export const RunCheckProvider: React.FC<RunCheckProviderProps> = ({
  children,
}) => {
  const [runTaskRequested, setRunTaskRequested] = useState<number | null>(null);
  const [pendingResolve, setPendingResolve] = useState<
    ((value: boolean) => void) | null
  >(null);

  const confirmTaskRun = useCallback((taskId: number): Promise<boolean> => {
    return new Promise((resolve) => {
      setRunTaskRequested(taskId);
      setPendingResolve(() => resolve);
    });
  }, []);

  const handleConfirm = () => {
    if (pendingResolve) {
      pendingResolve(true);
      setPendingResolve(null);
    }
    setRunTaskRequested(null);
  };

  const handleCancel = () => {
    if (pendingResolve) {
      pendingResolve(false);
      setPendingResolve(null);
    }
    setRunTaskRequested(null);
  };

  const contextValue = useMemo(
    () => ({
      runTaskRequested,
      setRunTaskRequested,
      confirmTaskRun,
    }),
    [runTaskRequested, confirmTaskRun]
  );

  return (
    <RunCheckContext.Provider value={contextValue}>
      {children}
      <ErrorAwareRunDialog
        runTaskRequested={runTaskRequested}
        handleConfirm={handleConfirm}
        handleCancel={handleCancel}
      />
    </RunCheckContext.Provider>
  );
};

/**
 * Converts a raw object path like
 *   "prosmart_refmac.container.controlParameters.RIGID_BODY_SELECTION[2].chainId"
 * into a user-friendly label like
 *   "RIGID_BODY_SELECTION[2]"
 */
const formatErrorPath = (path: string): string => {
  const parts = path.split(".");
  const skip = new Set(["container", "inputData", "controlParameters", "outputData"]);
  const meaningful = parts.filter((p, i) => {
    if (i === 0) return false; // task name
    return !skip.has(p);
  });
  if (meaningful.length <= 1) return meaningful[0] || path;
  return meaningful.slice(0, -1).join(" \u2192 ");
};

/**
 * Quick-action hints: when a warning objectPath matches one of these
 * suffixes, show a button that creates the helper task and navigates to it.
 */
const QUICK_ACTIONS: {
  suffix: string;
  taskName: string;
  label: string;
}[] = [
  { suffix: ".FREERFLAG", taskName: "freerflag", label: "Generate Free-R Flags" },
  { suffix: ".SEQIN",     taskName: "ProvideAsuContents", label: "Create ASU Contents" },
  { suffix: ".ASUFILE",   taskName: "ProvideAsuContents", label: "Create ASU Contents" },
  { suffix: ".ASUIN",     taskName: "ProvideAsuContents", label: "Create ASU Contents" },
];

interface ErrorAwareRunDialogProps {
  runTaskRequested: number | null;
  handleConfirm: () => void;
  handleCancel: () => void;
}

const ErrorAwareRunDialog: React.FC<ErrorAwareRunDialogProps> = ({
  runTaskRequested,
  handleConfirm,
  handleCancel,
}) => {
  const autoSubmitTimer = useRef<NodeJS.Timeout | null>(null);
  const { jobId } = useCCP4i2Window();
  const {
    validation, job, createPeerTask, mutateValidation,
    setParameter, fileItemToParameterArg, mutateContainer,
  } = useJob(jobId);
  const { jobs: projectJobs } = useProject(job?.project);
  const api = useApi();
  const { data: projects } = api.get<Project[]>(runTaskRequested !== null ? "projects" : null);

  // Fetch heavier run-time validation (monomer coverage etc.) only
  // when the run dialog is actually open.
  const { data: runTimeValidation, mutate: mutateRunTimeValidation } = api.get_validation(
    runTaskRequested !== null
      ? { type: "jobs", id: jobId, endpoint: "run_time_validation" }
      : null
  );

  // State for inline helper task modal (e.g. "Generate Free-R Flags")
  const [inlineTask, setInlineTask] = useState<{
    taskName: string;
    title: string;
    objectPath: string;
    helperJob: any;
  } | null>(null);

  // When the inline task completes, auto-select its output in the parent field
  const handleOutputFileReady = useCallback(
    async (outputFile: DjangoFile) => {
      if (!inlineTask?.objectPath || !projectJobs) return;
      try {
        const paramArg = fileItemToParameterArg(
          outputFile,
          inlineTask.objectPath,
          projectJobs,
          projects || []
        );
        await setParameter(paramArg);
        await mutateContainer();
      } catch (error) {
        console.error("Error auto-selecting output file:", error);
      }
      setInlineTask(null);
      await Promise.all([mutateValidation(), mutateRunTimeValidation()]);
    },
    [inlineTask?.objectPath, projectJobs, projects, fileItemToParameterArg,
     setParameter, mutateContainer, mutateValidation, mutateRunTimeValidation]
  );

  // Merge server validation sources: cheap polling + expensive run-time.
  // Run-time validation wins on key collisions (spread order).
  const receivedErrors: CCP4i2ErrorReport = {
    ...(validation || {}),
    ...(runTimeValidation || {}),
  };

  // seriousIssues: anything worth showing (WARNING or ERROR)
  const seriousIssues = Object.fromEntries(
    Object.entries(receivedErrors).filter(
      ([_, value]) => value.maxSeverity >= 1
    )
  );

  // blockingIssues: ERROR only — disables Confirm button
  const blockingIssues = Object.fromEntries(
    Object.entries(receivedErrors).filter(
      ([_, value]) => value.maxSeverity >= 2
    )
  );

  const hasSeriousIssues = Object.keys(seriousIssues).length > 0;
  const hasBlockingIssues = Object.keys(blockingIssues).length > 0;

  // Auto-submit only after run-time validation has loaded and found no issues.
  const runTimeValidationLoaded =
    runTaskRequested !== null && runTimeValidation !== undefined;

  useEffect(() => {
    if (autoSubmitTimer.current) {
      clearTimeout(autoSubmitTimer.current);
      autoSubmitTimer.current = null;
    }
    if (
      runTimeValidationLoaded &&
      !hasSeriousIssues &&
      runTaskRequested !== null &&
      jobId !== null &&
      jobId === runTaskRequested
    ) {
      autoSubmitTimer.current = setTimeout(() => {
        handleConfirm();
      }, 200);
    }
    return () => {
      if (autoSubmitTimer.current) {
        clearTimeout(autoSubmitTimer.current);
        autoSubmitTimer.current = null;
      }
    };
  }, [seriousIssues, runTaskRequested, jobId, runTimeValidationLoaded]);

  return (
    <>
    <Dialog
      open={runTaskRequested !== null && !inlineTask}
      onClose={() => handleCancel()}
      maxWidth="md"
      fullWidth
      slotProps={{ paper: { style: { minWidth: 600 } } }}
    >
      <DialogContent>
        <DialogTitle>Confirm Task Execution</DialogTitle>
        {hasSeriousIssues && (
          <div style={{ margin: "0.5rem 0" }}>
            {Object.entries(seriousIssues).map(([key, issueSet]) => {
              const quickAction = QUICK_ACTIONS.find((a) =>
                key.endsWith(a.suffix)
              );
              return issueSet.messages.map((issue, issueIndex) => (
                <div
                  key={`${key}_${issueIndex}`}
                  style={{
                    display: "flex",
                    alignItems: "center",
                    gap: "0.5rem",
                    padding: "0.35rem 0.75rem",
                    margin: "0.25rem 0",
                    borderLeft: `3px solid ${
                      issueSet.maxSeverity >= 2 ? "red" : "orange"
                    }`,
                    color: issueSet.maxSeverity >= 2 ? "red" : "inherit",
                    fontSize: "0.9rem",
                  }}
                >
                  <span style={{ flex: 1 }}>
                    <span style={{ fontWeight: 500 }}>
                      {formatErrorPath(key)}
                    </span>
                    {": "}
                    {issue}
                  </span>
                  {quickAction && issueIndex === 0 && (
                    <Button
                      size="small"
                      variant="outlined"
                      sx={{ textTransform: "none", whiteSpace: "nowrap", flexShrink: 0 }}
                      onClick={async () => {
                        const created = await createPeerTask(quickAction.taskName);
                        if (created) {
                          setInlineTask({
                            taskName: quickAction.taskName,
                            title: quickAction.label,
                            objectPath: key,
                            helperJob: created,
                          });
                        }
                      }}
                    >
                      {quickAction.label}
                    </Button>
                  )}
                </div>
              ));
            })}
          </div>
        )}
        <DialogActions>
          <Button onClick={handleCancel}>Cancel</Button>
          <Button onClick={handleConfirm} disabled={hasBlockingIssues}>
            Confirm
          </Button>
        </DialogActions>
      </DialogContent>
    </Dialog>

    {inlineTask && job && (
      <InlineTaskModal
        open={true}
        onClose={() => {
          setInlineTask(null);
          mutateValidation();
        }}
        taskName={inlineTask.taskName}
        parentJob={job}
        helperJobOverride={inlineTask.helperJob}
        onOutputFileReady={handleOutputFileReady}
        title={inlineTask.title}
      />
    )}
    </>
  );
};

export const useRunCheck = (): RunCheckContextType => {
  const context = useContext(RunCheckContext);
  if (!context) {
    throw new Error("useRunCheck must be used within a RunCheckProvider");
  }
  return context;
};
