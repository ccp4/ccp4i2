import React, {
  createContext,
  useCallback,
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
import { useRouter } from "next/navigation";
import { useCCP4i2Window } from "../app-context";
import { useJob } from "../utils";
import { Job } from "../types/models";

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

/**
 * Configuration for a single optional input warning.
 */
export interface OptionalInputWarningConfig {
  /** Unique key for this warning (e.g., "FREERFLAG", "ASUIN", "SEQIN") */
  key: string;
  /** The input value from useTaskItem - checked for dbFileId?.length > 0 */
  inputValue: any;
  /** Error key path (e.g., "prosmart_refmac.container.inputData.FREERFLAG") */
  errorKey: string;
  /** Warning messages to display */
  messages: string[];
  /** Task name to create when button is clicked (e.g., "freerflag", "ProvideAsuContents") */
  createTaskName: string;
  /** Button label (e.g., "Create FreeR task", "Create ASU Contents task") */
  buttonLabel: string;
}

/**
 * Configuration for the optional input warning hook.
 */
interface UseOptionalInputWarningOptions {
  /** The job object */
  job: Job;
  /** Server validation errors */
  validation: CCP4i2ErrorReport | null;
  /** Function to create a peer task */
  createPeerTask: (taskName: string) => Promise<Job | undefined>;
  /** Array of warning configurations for each optional input */
  warnings: OptionalInputWarningConfig[];
  /** Optional function to filter/transform validation errors before processing */
  filterErrors?: (errors: CCP4i2ErrorReport) => CCP4i2ErrorReport;
}

/**
 * Generic hook to manage warnings for optional-but-recommended inputs.
 *
 * This hook:
 * 1. Adds warnings (maxSeverity: 3) for unset optional inputs
 * 2. Adds action buttons to the run confirmation dialog for creating related tasks
 * 3. Removes warnings/buttons when inputs are set
 *
 * @example
 * ```tsx
 * const { value: freeRFlag } = useTaskItem("FREERFLAG");
 * const { value: asuIn } = useTaskItem("ASUIN");
 *
 * useOptionalInputWarning({
 *   job,
 *   validation,
 *   createPeerTask,
 *   warnings: [
 *     {
 *       key: "FREERFLAG",
 *       inputValue: freeRFlag,
 *       errorKey: "parrot.container.inputData.FREERFLAG",
 *       messages: ["Free R flag is recommended for refinement"],
 *       createTaskName: "freerflag",
 *       buttonLabel: "Create FreeR task",
 *     },
 *     {
 *       key: "ASUIN",
 *       inputValue: asuIn,
 *       errorKey: "parrot.container.inputData.ASUIN",
 *       messages: ["ASU content is recommended for accurate solvent estimation"],
 *       createTaskName: "ProvideAsuContents",
 *       buttonLabel: "Create ASU Contents task",
 *     },
 *   ],
 * });
 * ```
 */
export const useOptionalInputWarning = ({
  job,
  validation,
  createPeerTask,
  warnings,
  filterErrors,
}: UseOptionalInputWarningOptions): void => {
  const router = useRouter();
  const {
    processedErrors,
    setProcessedErrors,
    extraDialogActions,
    setExtraDialogActions,
    setRunTaskRequested,
  } = useRunCheck();

  // Create task handler factory
  const createTaskHandler = useCallback(
    (taskName: string) => async () => {
      try {
        const createdJob: Job | undefined = await createPeerTask(taskName);
        if (createdJob) {
          router.push(`/project/${job.project}/job/${createdJob.id}`);
          setRunTaskRequested(null);
        }
      } catch (error) {
        console.error(`Error creating ${taskName} task:`, error);
      }
    },
    [createPeerTask, job.project, router, setRunTaskRequested]
  );

  // Process validation errors - add warnings for unset inputs
  const processErrors = useCallback(() => {
    if (!validation) return;

    // Start with existing validation errors, optionally filtered
    let newProcessedErrors = { ...(validation as CCP4i2ErrorReport) };
    if (filterErrors) {
      newProcessedErrors = filterErrors(newProcessedErrors);
    }

    // Add warnings for each unset input
    for (const warning of warnings) {
      if (!warning.inputValue?.dbFileId?.length) {
        newProcessedErrors[warning.errorKey] = {
          messages: warning.messages,
          maxSeverity: 3, // Allows execution but shows warning
        };
      }
    }

    // Only update if errors have actually changed
    const newErrorsKey = JSON.stringify(newProcessedErrors);
    const currentErrorsKey = JSON.stringify(processedErrors);

    if (newErrorsKey !== currentErrorsKey) {
      setProcessedErrors(newProcessedErrors);
    }
  }, [validation, warnings, processedErrors, setProcessedErrors, filterErrors]);

  // Handle extra dialog actions
  const updateExtraDialogActions = useCallback(() => {
    let newActions = { ...extraDialogActions };
    let hasChanges = false;

    for (const warning of warnings) {
      const isUnset = !warning.inputValue?.dbFileId?.length;
      const hasAction = !!newActions?.[warning.key];

      if (isUnset && !hasAction) {
        // Add action button for unset input
        newActions[warning.key] = (
          <Button
            variant="contained"
            onClick={createTaskHandler(warning.createTaskName)}
          >
            {warning.buttonLabel}
          </Button>
        );
        hasChanges = true;
      } else if (!isUnset && hasAction) {
        // Remove action button when input is set
        delete newActions[warning.key];
        hasChanges = true;
      }
    }

    if (hasChanges) {
      setExtraDialogActions(
        Object.keys(newActions).length > 0 ? newActions : null
      );
    }
  }, [warnings, extraDialogActions, setExtraDialogActions, createTaskHandler]);

  // Effect for error processing
  useEffect(() => {
    const allDefined = warnings.every((w) => w.inputValue !== undefined);
    if (allDefined) {
      processErrors();
    }
  }, [warnings, processErrors]);

  // Effect for handling extra dialog actions
  useEffect(() => {
    updateExtraDialogActions();
  }, [updateExtraDialogActions]);
};

/**
 * Configuration for the FreeR warning hook.
 */
interface UseFreeRWarningOptions {
  /** The job object */
  job: Job;
  /** The task name prefix for the error key (e.g., "prosmart_refmac", "servalcat_pipe") */
  taskName: string;
  /** The FreeR flag value from useTaskItem("FREERFLAG") */
  freeRFlag: any;
  /** Server validation errors */
  validation: CCP4i2ErrorReport | null;
  /** Function to create a peer task */
  createPeerTask: (taskName: string) => Promise<Job | undefined>;
  /** Optional function to filter/transform validation errors before processing */
  filterErrors?: (errors: CCP4i2ErrorReport) => CCP4i2ErrorReport;
}

/**
 * Hook to manage FreeR flag warnings for refinement tasks.
 * This is a convenience wrapper around useOptionalInputWarning.
 *
 * @example
 * ```tsx
 * const { value: freeRFlag } = useTaskItem("FREERFLAG");
 * useFreeRWarning({
 *   job,
 *   taskName: "prosmart_refmac",
 *   freeRFlag,
 *   validation,
 *   createPeerTask,
 * });
 * ```
 */
export const useFreeRWarning = ({
  job,
  taskName,
  freeRFlag,
  validation,
  createPeerTask,
  filterErrors,
}: UseFreeRWarningOptions): void => {
  useOptionalInputWarning({
    job,
    validation,
    createPeerTask,
    filterErrors,
    warnings: [
      {
        key: "FREERFLAG",
        inputValue: freeRFlag,
        errorKey: `${taskName}.container.inputData.FREERFLAG`,
        messages: [
          "Setting the Free R flag file is strongly recommended for refinement",
          "You are advised to select an existing set or create a new one",
        ],
        createTaskName: "freerflag",
        buttonLabel: "Create FreeR task",
      },
    ],
  });
};

/**
 * Configuration for the ASU content warning hook.
 */
interface UseAsuContentWarningOptions {
  /** The job object */
  job: Job;
  /** The task name prefix for the error key (e.g., "parrot", "modelcraft") */
  taskName: string;
  /** The ASU content value from useTaskItem("ASUIN") or similar */
  asuContent: any;
  /** Server validation errors */
  validation: CCP4i2ErrorReport | null;
  /** Function to create a peer task */
  createPeerTask: (taskName: string) => Promise<Job | undefined>;
  /** Optional function to filter/transform validation errors before processing */
  filterErrors?: (errors: CCP4i2ErrorReport) => CCP4i2ErrorReport;
}

/**
 * Hook to manage ASU content warnings for phasing/building tasks.
 * This is a convenience wrapper around useOptionalInputWarning.
 *
 * @example
 * ```tsx
 * const { value: asuIn } = useTaskItem("ASUIN");
 * useAsuContentWarning({
 *   job,
 *   taskName: "parrot",
 *   asuContent: asuIn,
 *   validation,
 *   createPeerTask,
 * });
 * ```
 */
export const useAsuContentWarning = ({
  job,
  taskName,
  asuContent,
  validation,
  createPeerTask,
  filterErrors,
}: UseAsuContentWarningOptions): void => {
  useOptionalInputWarning({
    job,
    validation,
    createPeerTask,
    filterErrors,
    warnings: [
      {
        key: "ASUIN",
        inputValue: asuContent,
        errorKey: `${taskName}.container.inputData.ASUIN`,
        messages: [
          "Providing ASU content information is strongly recommended",
          "It improves solvent estimation and Matthews coefficient calculation",
        ],
        createTaskName: "ProvideAsuContents",
        buttonLabel: "Create ASU Contents task",
      },
    ],
  });
};

/**
 * Configuration for the sequence warning hook.
 */
interface UseSequenceWarningOptions {
  /** The job object */
  job: Job;
  /** The task name prefix for the error key (e.g., "shelx") */
  taskName: string;
  /** The sequence value from useTaskItem("SEQIN") or similar */
  sequence: any;
  /** Server validation errors */
  validation: CCP4i2ErrorReport | null;
  /** Function to create a peer task */
  createPeerTask: (taskName: string) => Promise<Job | undefined>;
  /** Optional function to filter/transform validation errors before processing */
  filterErrors?: (errors: CCP4i2ErrorReport) => CCP4i2ErrorReport;
}

/**
 * Hook to manage sequence input warnings for tasks that benefit from sequence information.
 * This is a convenience wrapper around useOptionalInputWarning.
 *
 * @example
 * ```tsx
 * const { value: seqIn } = useTaskItem("SEQIN");
 * useSequenceWarning({
 *   job,
 *   taskName: "shelx",
 *   sequence: seqIn,
 *   validation,
 *   createPeerTask,
 * });
 * ```
 */
export const useSequenceWarning = ({
  job,
  taskName,
  sequence,
  validation,
  createPeerTask,
  filterErrors,
}: UseSequenceWarningOptions): void => {
  useOptionalInputWarning({
    job,
    validation,
    createPeerTask,
    filterErrors,
    warnings: [
      {
        key: "SEQIN",
        inputValue: sequence,
        errorKey: `${taskName}.container.inputData.SEQIN`,
        messages: [
          "Providing sequence information is strongly recommended",
          "It enables molecular weight calculation and improves model building",
        ],
        createTaskName: "ProvideAsuContents",
        buttonLabel: "Create ASU Contents task",
      },
    ],
  });
};
