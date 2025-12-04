import React, {
  createContext,
  useContext,
  useState,
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
import { useJob } from "../utils";

/**
 * An error report keyed by parameter name or error type.
 * Each entry contains a maxSeverity (number) and an array of messages (strings).
 */
export interface CCP4i2ErrorReport {
  [key: string]: {
    maxSeverity: number;
    messages: string[];
    // You can add more fields here if needed
  };
}

/**
 * An object whose keys are strings and values are ReactNode actions.
 */
export interface CCP4i2RunActions {
  [key: string]: ReactNode;
}

interface RunCheckContextType {
  runTaskRequested: number | null;
  setRunTaskRequested: (taskId: number | null) => void;
  confirmTaskRun: (taskId: number) => Promise<boolean>;
  extraDialogActions: CCP4i2RunActions | null;
  setExtraDialogActions: (actions: CCP4i2RunActions | null) => void;
  processedErrors: CCP4i2ErrorReport | null;
  setProcessedErrors: (errors: CCP4i2ErrorReport | null) => void;
}

export const RunCheckContext = createContext<RunCheckContextType>({
  runTaskRequested: null,
  setRunTaskRequested: () => {},
  confirmTaskRun: () => Promise.resolve(false),
  extraDialogActions: null,
  setExtraDialogActions: () => {},
  processedErrors: null,
  setProcessedErrors: () => {},
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
  const [extraDialogActions, setExtraDialogActions] =
    useState<CCP4i2RunActions | null>(null);
  const [processedErrors, setProcessedErrors] =
    useState<CCP4i2ErrorReport | null>(null);

  const confirmTaskRun = (taskId: number): Promise<boolean> => {
    return new Promise((resolve) => {
      setRunTaskRequested(taskId);
      setPendingResolve(() => resolve);
    });
  };

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

  return (
    <RunCheckContext.Provider
      value={
        {
          runTaskRequested,
          setRunTaskRequested,
          confirmTaskRun,
          extraDialogActions,
          setExtraDialogActions,
          processedErrors,
          setProcessedErrors,
        } as RunCheckContextType
      }
    >
      {children}
      <ErrorAwareRunDialog
        runTaskRequested={runTaskRequested}
        handleConfirm={handleConfirm}
        handleCancel={handleCancel}
      />
    </RunCheckContext.Provider>
  );
};

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
  const { extraDialogActions, processedErrors } = useRunCheck();
  const autoSubmitTimer = useRef<NodeJS.Timeout | null>(null);
  const { jobId } = useCCP4i2Window();
  const { validation } = useJob(jobId);

  // ...inside ErrorAwareRunDialog...
  const receivedErrors: CCP4i2ErrorReport | null =
    processedErrors || validation || {};

  const seriousIssues: CCP4i2ErrorReport | null = receivedErrors
    ? Object.fromEntries(
        Object.entries(receivedErrors).filter(
          ([_, value]) => value.maxSeverity === 2 || value.maxSeverity === 3
        )
      )
    : null;

  const blockingIssues: CCP4i2ErrorReport | null = receivedErrors
    ? Object.fromEntries(
        Object.entries(receivedErrors).filter(
          ([_, value]) => value.maxSeverity === 2
        )
      )
    : null;

  const hasSeriousIssues = seriousIssues
    ? Object.keys(seriousIssues).length > 0
    : false;

  useEffect(() => {
    if (autoSubmitTimer.current) {
      clearTimeout(autoSubmitTimer.current);
      autoSubmitTimer.current;
    }
    if (
      !hasSeriousIssues &&
      runTaskRequested !== null &&
      jobId !== null &&
      jobId === runTaskRequested
    ) {
      autoSubmitTimer.current = setTimeout(() => {
        handleConfirm();
      }, 200); // Auto-submit after 200 milliseconds if no serious issues
    }
    return () => {
      if (autoSubmitTimer.current) {
        clearTimeout(autoSubmitTimer.current);
        autoSubmitTimer.current = null;
      }
    };
  }, [seriousIssues, runTaskRequested, jobId]);

  return (
    <Dialog
      open={runTaskRequested !== null}
      onClose={() => handleCancel()}
      maxWidth="md"
      fullWidth
      slotProps={{
        paper: {
          style: { minWidth: 600 },
        },
      }}
    >
      <DialogContent>
        <DialogTitle>Confirm Task Execution</DialogTitle>
        {seriousIssues && Object.keys(seriousIssues).length > 0 && (
          <pre style={{ color: "red" }}>
            {Object.entries(seriousIssues).map(([key, issueSet], index) =>
              issueSet.messages.map((issue, issueIndex) => (
                <div key={`${key}_${issueIndex}`}>{key}: {issue}</div>
              ))
            )}
          </pre>
        )}
        <DialogActions>
          {extraDialogActions &&
            Object.keys(extraDialogActions)?.map((actionName, index) => (
              <React.Fragment key={index}>
                {extraDialogActions[actionName]}
              </React.Fragment> // Ensure each action is wrapped in a fragment
            ))}
          <Button onClick={handleCancel}>Cancel</Button>
          <Button
            onClick={handleConfirm}
            disabled={
              blockingIssues !== null && Object.keys(blockingIssues).length > 0
            }
          >
            Confirm
          </Button>
        </DialogActions>
      </DialogContent>
    </Dialog>
  );
};

export const useRunCheck = (): RunCheckContextType => {
  const context = useContext(RunCheckContext);
  if (!context) {
    throw new Error("useRunCheck must be used within a RunCheckProvider");
  }
  return context;
};
