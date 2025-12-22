import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import $ from "jquery";
import useSWR, { KeyedMutator, mutate, SWRResponse } from "swr";

import { useApi } from "./api";
import {
  Job,
  Project,
  File as DjangoFile,
} from "./types/models";
import { useRunCheck } from "./providers/run-check-provider";
import { useParameterChangeIntent } from "./providers/parameter-change-intent-provider";
import { apiJson, apiText } from "./api-fetch";
import { useIsJobEffectivelyActive } from "./providers/recently-started-jobs-context";

// ============================================================================
// Types and Interfaces
// ============================================================================

export interface SetParameterArg {
  object_path: string;
  value: any;
}

/**
 * New API response format for set_parameter endpoint.
 * Returns {success: true, data: {updated_item: ...}} on success.
 */
export type SetParameterResponse =
  | {
      success: true;
      data: {
        updated_item: any;
      };
    }
  | {
      success: false;
      error: string;
    };

/**
 * New API response format for create_task endpoint.
 * Returns {success: true, data: {new_job: ...}} on success.
 */
export interface CreateTaskResponse {
  success: boolean;
  data?: {
    new_job: Job;
  };
  error?: string;
}
export interface ValidationError {
  path: string;
  error: {
    maxSeverity: number;
    [key: string]: any;
  };
}

export interface TaskItem {
  item: any;
  value: any;
  update: (value: any) => Promise<boolean | Response>;
  updateNoMutate: (value: any) => Promise<boolean | Response>;
}

export interface ProjectData {
  project: Project | undefined;
  mutateProject: () => void;
  directory: any;
  mutateDirectory: () => void;
  jobs: Job[] | undefined;
  mutateJobs: KeyedMutator<Job[]>;
  mutateJobTree: () => void;
  /** Combined mutation that invalidates both jobs and job_tree endpoints */
  mutateAllJobs: () => void;
  files: DjangoFile[] | undefined;
  mutateFiles: KeyedMutator<DjangoFile[]>;
}

/**
 * API response format for upload_file_param endpoint.
 */
export type UploadFileParamResponse =
  | {
      success: true;
      data: {
        updated_item: any;
      };
    }
  | {
      success: false;
      error: string;
    };

/**
 * Column selection for multi-representation MTZ upload.
 */
export interface ColumnSelectorEntry {
  signature: string;
  columnSelector: string;
  contentFlag: number;
  fileSuffix: string;
  isPrimary: boolean;
}

/**
 * Arguments for uploadFileParam function.
 */
export interface UploadFileParamArg {
  objectPath: string;
  file: Blob;
  fileName: string;
  /** Legacy single column selector (backward compatible) */
  columnSelector?: string;
  /** Enhanced multi-selector format for multiple representations */
  columnSelectors?: ColumnSelectorEntry[];
}

export interface JobData {
  job: Job | undefined;
  mutateJob: () => void;
  container: any;
  mutateContainer: () => void;
  params_xml: any;
  mutateParams_xml: () => void;
  validation: any;
  mutateValidation: () => void;
  diagnostic_xml: any;
  mutateDiagnosticXml: () => void;
  def_xml: any;
  mutateDef_xml: () => void;
  setParameter: (
    arg: SetParameterArg
  ) => Promise<SetParameterResponse | undefined>;
  setParameterNoMutate: (
    arg: SetParameterArg
  ) => Promise<SetParameterResponse | undefined>;
  uploadFileParam: (
    arg: UploadFileParamArg
  ) => Promise<UploadFileParamResponse | undefined>;
  useTaskItem: (paramName: string) => TaskItem;
  createPeerTask: (taskName: string) => Promise<Job | undefined>;
  useFileContent: (paramName: string) => SWRResponse<string, Error>;
  getValidationColor: (item: any) => string;
  getErrors: (item: any) => ValidationError[];
  useFileDigest: (objectPath: string) => SWRResponse<any, Error>;
  fetchDigest: (objectPath: string) => Promise<any | null>;
  fileItemToParameterArg: (
    value: DjangoFile,
    objectPath: string,
    projectJobs: Job[],
    projects: Project[]
  ) => SetParameterArg;
}

// ============================================================================
// Constants
// ============================================================================

const VALIDATION_COLORS = {
  SUCCESS: "success.light",
  WARNING: "warning.light",
  ERROR: "error.light",
} as const;

const SEVERITY_LEVELS = {
  SUCCESS: 0,
  WARNING: 1,
  ERROR: 2,
} as const;

const JOB_STATUS = {
  PENDING: 1,
} as const;

const FILE_DIRECTORIES = {
  JOB_OUTPUT: 1,
  IMPORTED: 2,
} as const;

const FILE_PATHS = {
  IMPORTED_FILES: "CCP4_IMPORTED_FILES",
  JOB_PREFIX: "job_",
  JOBS_DIR: "CCP4_JOBS",
} as const;

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Checks if the given path ends with the specified name, considering dot notation.
 */
const pathMatches = (
  path: string | null | undefined,
  name: string
): boolean => {
  if (!path || !name) return false;

  const normalizedName = `.${name}`.replace(/\.+/g, ".");
  const normalizedPath = `.${path}`.replace(/\.+/g, ".");

  return normalizedPath.endsWith(normalizedName);
};

/**
 * Safely checks if an object has a constructor matching the expected type.
 */
const hasConstructor = (obj: any, constructor: any): boolean => {
  return obj && typeof obj === "object" && obj.constructor === constructor;
};

/**
 * Recursively searches for items within a container that match a specified name.
 */
const findItemsRecursively = (
  name: string,
  container: any,
  multiple: boolean = true,
  accumulator: any[] = []
): any[] => {
  if (!container || !name) return accumulator;

  const initialLength = accumulator.length;

  try {
    // Check if current container matches
    if (pathMatches(container?._objectPath, name)) {
      accumulator.push(container);
      if (!multiple) return accumulator;
    }

    // Handle CList containers
    if (container._baseClass === "CList" && Array.isArray(container._value)) {
      for (const item of container._value) {
        if (pathMatches(item?._objectPath, name)) {
          accumulator.push(item);
          if (!multiple) return accumulator;
        }

        // Recursive search
        findItemsRecursively(name, item, multiple, accumulator);
        if (!multiple && accumulator.length > initialLength) {
          return accumulator;
        }
      }
    }
    // Handle object containers
    else if (hasConstructor(container._value, Object)) {
      for (const key of Object.keys(container._value)) {
        const item = container._value[key];

        if (pathMatches(item?._objectPath, name)) {
          accumulator.push(item);
          if (!multiple) return accumulator;
        }

        // Recursive search
        findItemsRecursively(name, item, multiple, accumulator);
        if (!multiple && accumulator.length > initialLength) {
          return accumulator;
        }
      }
    }
  } catch (error) {
    console.error(`Error searching for items with name "${name}":`, error);
  }

  return accumulator;
};

/**
 * Extracts validation errors for a given item based on the provided validation object.
 */
const extractValidationErrors = (
  item: any,
  validation: any
): ValidationError[] => {
  if (!validation || !item?._objectPath) {
    return [];
  }

  const itemPath = item._objectPath;
  const matchingErrors: ValidationError[] = [];

  try {
    for (const validationPath of Object.keys(validation)) {
      if (
        validationPath === itemPath ||
        validationPath.startsWith(`${itemPath}.`) ||
        validationPath.startsWith(`${itemPath}[`)
      ) {
        matchingErrors.push({
          path: validationPath,
          error: validation[validationPath],
        });
      }
    }
  } catch (error) {
    console.error("Error extracting validation errors:", error);
  }

  return matchingErrors;
};

/**
 * Determines the appropriate validation color based on error severity.
 */
const determineValidationColor = (
  fieldErrors: ValidationError[] | any
): string => {
  if (
    !fieldErrors ||
    !Array.isArray(fieldErrors) ||
    (Array.isArray(fieldErrors) && fieldErrors.length === 0)
  ) {
    return VALIDATION_COLORS.SUCCESS;
  }

  let maxSeverity: number = SEVERITY_LEVELS.SUCCESS;

  try {
    if (Array.isArray(fieldErrors)) {
      maxSeverity = fieldErrors.reduce((highest, error) => {
        const currentSeverity =
          error?.error?.maxSeverity ?? SEVERITY_LEVELS.SUCCESS;
        return Math.max(highest, currentSeverity);
      }, SEVERITY_LEVELS.SUCCESS);
    } else {
      if (
        fieldErrors &&
        typeof fieldErrors === "object" &&
        "maxSeverity" in fieldErrors
      ) {
        maxSeverity =
          (fieldErrors as { maxSeverity: number }).maxSeverity ??
          SEVERITY_LEVELS.SUCCESS;
      } else {
        maxSeverity = SEVERITY_LEVELS.SUCCESS;
      }
    }
  } catch (error) {
    console.error("Error determining validation color:", error);
    return VALIDATION_COLORS.ERROR;
  }

  if (maxSeverity === SEVERITY_LEVELS.SUCCESS) {
    return VALIDATION_COLORS.SUCCESS;
  } else if (maxSeverity === SEVERITY_LEVELS.WARNING) {
    return VALIDATION_COLORS.WARNING;
  } else {
    return VALIDATION_COLORS.ERROR;
  }
};

/**
 * Recursively extracts values from complex data structures.
 */
const extractValueRecursively = (item: any): any => {
  if (!item) return null;

  const { _value } = item;

  // Handle primitive values
  if (
    _value === undefined ||
    _value === null ||
    typeof _value === "string" ||
    typeof _value === "number" ||
    typeof _value === "boolean"
  ) {
    return _value;
  }

  // Handle objects
  if (hasConstructor(_value, Object)) {
    const result: Record<string, any> = {};
    try {
      for (const key of Object.keys(_value)) {
        result[key] = extractValueRecursively(_value[key]);
      }
    } catch (error) {
      console.error("Error extracting object values:", error);
    }
    return result;
  }

  // Handle arrays
  if (Array.isArray(_value)) {
    if (_value.length === 0) return [];

    try {
      return _value.map((value) => extractValueRecursively(value));
    } catch (error) {
      console.error("Error extracting array values:", error);
      return [];
    }
  }

  console.warn("Unknown item type:", _value);
  return _value;
};

/**
 * Creates a safe file reader promise with proper error handling.
 */
const createFileReaderPromise = (
  file: File,
  readAs: "Text" | "ArrayBuffer" | "File" = "Text"
): Promise<string | ArrayBuffer | File> => {
  return new Promise((resolve, reject) => {
    if (readAs === "File") {
      resolve(file);
      return;
    }

    const reader = new FileReader();

    const cleanup = () => {
      reader.onabort = null;
      reader.onerror = null;
      reader.onloadend = null;
    };

    reader.onabort = () => {
      cleanup();
      reject(new Error("File reading was aborted"));
    };

    reader.onerror = () => {
      cleanup();
      reject(new Error("File reading failed"));
    };

    reader.onloadend = () => {
      cleanup();
      if (reader.result !== null) {
        resolve(reader.result);
      } else {
        reject(new Error("File reading returned null"));
      }
    };

    try {
      if (readAs === "Text") {
        reader.readAsText(file);
      } else if (readAs === "ArrayBuffer") {
        reader.readAsArrayBuffer(file);
      }
    } catch (error) {
      cleanup();
      reject(error);
    }
  });
};

/**
 * Prettifies XML using XSLT transformation with error handling.
 */
const prettifyXmlSafely = (sourceXml: Document | Element): string => {
  if (!sourceXml) return "";

  try {
    let targetNode: Document | Element | undefined = sourceXml;

    // Handle jQuery nodes
    if (!targetNode?.nodeName) {
      try {
        const node = $(sourceXml).get(0);
        if (node && node.nodeType === Node.ELEMENT_NODE) {
          targetNode = node as Element;
        } else if (node && node.nodeType === Node.DOCUMENT_NODE) {
          targetNode = node as Document;
        } else {
          console.error(
            "Cannot extract HTML element from jQuery object: unexpected node type"
          );
          return "";
        }
      } catch (error) {
        console.error("Cannot extract HTML element from jQuery object:", error);
        return "";
      }
    }

    if (!targetNode) return "";

    const xsltStylesheet = [
      '<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">',
      '  <xsl:strip-space elements="*"/>',
      '  <xsl:template match="para[content-style][not(text())]">',
      '    <xsl:value-of select="normalize-space(.)"/>',
      "  </xsl:template>",
      '  <xsl:template match="node()|@*">',
      '    <xsl:copy><xsl:apply-templates select="node()|@*"/></xsl:copy>',
      "  </xsl:template>",
      '  <xsl:output indent="yes"/>',
      "</xsl:stylesheet>",
    ].join("\n");

    const xsltDoc = new DOMParser().parseFromString(
      xsltStylesheet,
      "application/xml"
    );
    const xsltProcessor = new XSLTProcessor();
    xsltProcessor.importStylesheet(xsltDoc);

    const resultDoc = xsltProcessor.transformToDocument(targetNode);
    return new XMLSerializer().serializeToString(resultDoc);
  } catch (error) {
    console.error("Error prettifying XML:", error);
    return "";
  }
};

// ============================================================================
// Exported Functions
// ============================================================================

/**
 * Retrieves items from a container that match the specified name.
 */
export const itemsForName = (
  name: string,
  container: any,
  multiple: boolean = true
): any[] => {
  if (!name || !container) return [];
  return findItemsRecursively(name, container, multiple);
};

/**
 * Extracts the value of an item, handling various data types.
 */
export const valueOfItem = extractValueRecursively;

/**
 * Determines the appropriate validation color based on field errors.
 */
export const validationColor = determineValidationColor;

/**
 * Reads file contents with proper error handling.
 */
export const readFilePromise = createFileReaderPromise;

/**
 * Prettifies XML with error handling.
 */
export const prettifyXml = prettifyXmlSafely;

/**
 * Custom hook for async effects with proper cleanup.
 */
export const useAsyncEffect = (
  effect: () => Promise<void>,
  dependencies: any[]
): void => {
  useEffect(() => {
    let cancelled = false;

    const executeEffect = async () => {
      try {
        if (!cancelled) {
          await effect();
        }
      } catch (error) {
        if (!cancelled) {
          console.error("Async effect error:", error);
        }
      }
    };

    executeEffect();

    return () => {
      cancelled = true;
    };
  }, dependencies);
};

// ============================================================================
// Custom Hooks
// ============================================================================

/**
 * Custom hook that returns the previous value of the given input.
 */
export const usePrevious = <T>(value: T): T | undefined => {
  const ref = useRef<T | undefined>(undefined);

  useEffect(() => {
    ref.current = value;
  }, [value]);

  return ref.current;
};

/**
 * Custom hook to fetch and manage project-related data.
 * Accepts undefined/null to skip fetching until projectId is available.
 */
export const useProject = (projectId: number | null | undefined): ProjectData => {
  const api = useApi();

  const { data: project, mutate: mutateProject } = api.get_endpoint<Project>({
    type: "projects",
    id: projectId,
    endpoint: "",
  });

  const { data: directory, mutate: mutateDirectory } = api.get_endpoint<any>({
    type: "projects",
    id: projectId,
    endpoint: "directory",
  });

  const { data: jobs, mutate: mutateJobs } = api.get_endpoint<Job[]>({
    type: "projects",
    id: projectId,
    endpoint: "jobs",
  });

  // Also get mutator for job_tree endpoint (used by classic-jobs-list)
  const { mutate: mutateJobTree } = api.get_endpoint<any>({
    type: "projects",
    id: projectId,
    endpoint: "job_tree",
  });

  // Combined mutation that invalidates both jobs and job_tree endpoints
  const mutateAllJobs = useCallback(() => {
    mutateJobs();
    mutateJobTree();
  }, [mutateJobs, mutateJobTree]);

  const { data: files, mutate: mutateFiles } = api.get_endpoint<DjangoFile[]>({
    type: "projects",
    id: projectId,
    endpoint: "files",
  });

  // Note: KPIs are now embedded in job_tree endpoint, so we no longer need
  // separate job_float_values and job_char_values API calls

  return {
    project,
    mutateProject,
    directory,
    mutateDirectory,
    jobs,
    mutateJobs,
    mutateJobTree,
    mutateAllJobs,
    files,
    mutateFiles,
  };
};

// ============================================================================
// Project Files Hook with Optional Polling
// ============================================================================

/** Polling interval for files when task tab is open */
const FILES_POLL_INTERVAL = 5000;

/**
 * Hook for fetching project files with optional polling.
 *
 * When a task widget is open, files should be polled to pick up newly created
 * output files from completed jobs. Pass `shouldPoll=true` when the task tab
 * is active to enable 5-second polling.
 *
 * @param projectId - The project ID to fetch files for
 * @param shouldPoll - Whether to enable polling (e.g., when task tab is open)
 */
export const useProjectFiles = (
  projectId: number | undefined,
  shouldPoll: boolean = false
) => {
  const api = useApi();

  const pollInterval = shouldPoll ? FILES_POLL_INTERVAL : 0;

  const { data: files, mutate: mutateFiles } = api.get_endpoint<DjangoFile[]>(
    projectId
      ? {
          type: "projects",
          id: projectId,
          endpoint: "files",
        }
      : null,
    pollInterval
  );

  return {
    files,
    mutateFiles,
  };
};

// ============================================================================
// Job-aware Directory Hook
// ============================================================================

/** Polling interval when active job is running */
const DIRECTORY_ACTIVE_POLL_INTERVAL = 5000;
/** Job statuses that indicate active work: Queued (2), Running (3), Running remotely (7) */
const ACTIVE_JOB_STATUSES = [2, 3, 7];

/**
 * Hook for fetching directory data with job-aware polling.
 *
 * - Polls frequently (5s) when the active job is running/queued OR was recently started
 * - Stops polling when job is idle
 * - Forces a refresh when job transitions from running to not-running
 * - Uses grace period tracking to handle DB latency after job submission
 *
 * @param projectId - The project ID to fetch directory for
 * @param activeJob - The currently active job (optional)
 */
export const useJobDirectory = (
  projectId: number | undefined,
  activeJob: Job | undefined | null
) => {
  const api = useApi();

  // Use grace period-aware hook to handle DB latency after job submission
  // This also updates the observed status for stall detection
  const isJobActive = useIsJobEffectivelyActive(activeJob?.id, activeJob?.status);

  // Track previous job status for transition detection
  const previousJob = usePrevious(activeJob);

  // Use adaptive polling - fast when job is running, no polling when idle
  const pollInterval = isJobActive ? DIRECTORY_ACTIVE_POLL_INTERVAL : 0;

  const { data: directory, mutate: mutateDirectory } = api.get_endpoint<any>(
    projectId
      ? {
          type: "projects",
          id: projectId,
          endpoint: "directory",
        }
      : null,
    pollInterval
  );

  // Detect status transitions from running to not-running and force refresh
  useEffect(() => {
    if (!activeJob || !previousJob) return;

    // Check if job transitioned from running to finished
    const wasActive = ACTIVE_JOB_STATUSES.includes(previousJob.status);
    const isNowIdle = !ACTIVE_JOB_STATUSES.includes(activeJob.status);

    // Only refresh if this is the same job (not a job switch) and it just finished
    if (wasActive && isNowIdle && activeJob.id === previousJob.id) {
      console.log(
        `[useJobDirectory] Job ${activeJob.id} finished (${previousJob.status} -> ${activeJob.status}), refreshing directory`
      );
      mutateDirectory();
    }
  }, [activeJob, previousJob, mutateDirectory]);

  return {
    directory,
    mutateDirectory,
    isJobActive,
  };
};

// ============================================================================
// Parameter Setting Queue
// ============================================================================

interface QueuedParameterOperation {
  id: string;
  operation: () => Promise<SetParameterResponse | undefined>;
  resolve: (value: SetParameterResponse | undefined) => void;
  reject: (error: any) => void;
}

class ParameterQueue {
  private queue: QueuedParameterOperation[] = [];
  private isProcessing = false;

  async enqueue(
    operation: () => Promise<SetParameterResponse | undefined>
  ): Promise<SetParameterResponse | undefined> {
    return new Promise((resolve, reject) => {
      const queueItem: QueuedParameterOperation = {
        id: `param_${Date.now()}_${Math.random()}`,
        operation,
        resolve,
        reject,
      };

      this.queue.push(queueItem);
      this.processQueue();
    });
  }

  private async processQueue(): Promise<void> {
    if (this.isProcessing || this.queue.length === 0) {
      return;
    }

    this.isProcessing = true;

    while (this.queue.length > 0) {
      const item = this.queue.shift();
      if (!item) break;

      try {
        console.log(`Processing parameter operation ${item.id}`);
        const result = await item.operation();
        item.resolve(result);
      } catch (error) {
        console.error(`Error in parameter operation ${item.id}:`, error);
        item.reject(error);
      }
    }

    this.isProcessing = false;
  }

  getQueueLength(): number {
    return this.queue.length;
  }

  isQueueProcessing(): boolean {
    return this.isProcessing;
  }
}

// Global parameter queue instance
const parameterQueue = new ParameterQueue();

// ============================================================================
// Modified useJob Hook
// ============================================================================

export const useJob = (jobId: number | null | undefined): JobData => {
  const api = useApi();

  const { data: job, mutate: mutateJob } = api.get_endpoint<Job>(
    {
      type: "jobs",
      id: jobId,
      endpoint: "",
    },
    10000
  );

  const { data: container, mutate: mutateContainer } =
    api.get_wrapped_endpoint_json<any>({
      type: "jobs",
      id: jobId,
      endpoint: "container",
    });

  const { data: params_xml, mutate: mutateParams_xml } =
    api.get_pretty_endpoint_xml({
      type: "jobs",
      id: jobId,
      endpoint: "params_xml",
    });

  const { processedErrors, setProcessedErrors } = useRunCheck();

  // Get mutateValidation from useSWR
  const { data: validation, mutate: mutateValidation } = api.get_validation({
    type: "jobs",
    id: jobId,
    endpoint: "validation",
  });

  // Decorate mutateValidation so it always resets processed errors
  const mutateValidationWithProcessedErrors = useCallback(
    async (...args: any[]) => {
      setProcessedErrors(null);
      return mutateValidation(...args);
    },
    [mutateValidation, setProcessedErrors]
  );

  // Only fetch diagnostic_xml for jobs that have run (status > 1)
  // Status 1 = PENDING, Status > 1 = has run or is running
  const shouldFetchDiagnostic = job?.status && job.status > JOB_STATUS.PENDING;
  const { data: diagnostic_xml, mutate: mutateDiagnosticXml } =
    api.get_pretty_endpoint_xml(
      shouldFetchDiagnostic
        ? {
            type: "jobs",
            id: jobId,
            endpoint: "diagnostic_xml",
          }
        : null
    );

  const { data: def_xml, mutate: mutateDef_xml } = api.get_pretty_endpoint_xml({
    type: "jobs",
    id: jobId,
    endpoint: "def_xml",
  });

  const { mutateJobs } = useProject(job?.project);
  const { setIntent, setIntentForPath, clearIntentForPath } = useParameterChangeIntent();

  // Memoized functions
  const setParameter = useCallback(
    async (
      setParameterArg: SetParameterArg
    ): Promise<SetParameterResponse | undefined> => {
      if (job?.status !== JOB_STATUS.PENDING) {
        console.warn(
          "Attempting to edit interface of task not in pending state"
        );
        return undefined;
      }

      const objectPath = setParameterArg.object_path;

      // Record intent BEFORE making the API call
      // This prevents the container refetch from overwriting local state
      // Get previous value from container lookup if available
      const previousValue = container?.lookup?.[objectPath]?._value;
      setIntentForPath({
        jobId: job.id,
        parameterPath: objectPath,
        reason: "UserEdit",
        previousValue,
      });

      // Enqueue the operation to ensure sequential execution
      return parameterQueue.enqueue(async () => {
        try {
          console.log(
            "Executing setParameter for:",
            objectPath
          );

          const result = await api.post<SetParameterResponse>(
            `jobs/${job.id}/set_parameter`,
            setParameterArg
          );
          setProcessedErrors(null);
          // Update all related data
          await Promise.all([
            mutateValidation(),
            mutateContainer(),
            mutateParams_xml(),
          ]);

          // NOTE: Don't clear intent here - let it expire naturally via auto-cleanup
          // The intent mechanism persists for a short window after changes

          console.log("Parameter set successfully:", result);
          return result;
        } catch (error) {
          // Clear intent on error so future syncs work
          clearIntentForPath(objectPath);
          console.error("Error setting parameter:", error);
          throw error;
        }
      });
    },
    [
      job,
      container,
      mutateContainer,
      mutateValidation,
      mutateParams_xml,
      api,
      setProcessedErrors,
      setIntentForPath,
      clearIntentForPath,
    ]
  );

  const setParameterNoMutate = useCallback(
    async (
      setParameterArg: SetParameterArg
    ): Promise<SetParameterResponse | undefined> => {
      if (job?.status !== JOB_STATUS.PENDING) {
        console.warn(
          "Attempting to edit interface of task not in pending state"
        );
        return undefined;
      }

      const objectPath = setParameterArg.object_path;

      // Record intent even for no-mutate calls
      // This is important for derived updates (e.g., CImportUnmergedElement)
      // where multiple fields are updated before a single mutateContainer()
      const previousValue = container?.lookup?.[objectPath]?._value;
      setIntentForPath({
        jobId: job.id,
        parameterPath: objectPath,
        reason: "UserEdit",
        previousValue,
      });

      // Enqueue the operation to ensure sequential execution
      return parameterQueue.enqueue(async () => {
        try {
          console.log(
            "Executing setParameterNoMutate for:",
            objectPath
          );

          const result = await api.post<SetParameterResponse>(
            `jobs/${job.id}/set_parameter`,
            setParameterArg
          );

          // Note: We don't clear intent here because the caller will typically
          // call mutateContainer() later, and we want the intent to persist
          // until that refetch completes. The auto-cleanup will handle stale intents.

          return result;
        } catch (error) {
          // Clear intent on error
          clearIntentForPath(objectPath);
          console.error("Error setting parameter (no mutate):", error);
          throw error;
        }
      });
    },
    [job, container, mutateParams_xml, mutateValidation, api, setIntentForPath, clearIntentForPath]
  );

  /**
   * Upload a file to a CDataFile parameter with intent tracking.
   * This is the centralized function for all file uploads that should
   * integrate with the intent tracking mechanism.
   */
  const uploadFileParam = useCallback(
    async (
      uploadArg: UploadFileParamArg
    ): Promise<UploadFileParamResponse | undefined> => {
      if (job?.status !== JOB_STATUS.PENDING) {
        console.warn(
          "Attempting to upload file to task not in pending state"
        );
        return undefined;
      }

      const { objectPath, file, fileName, columnSelector, columnSelectors } = uploadArg;

      // Record intent BEFORE making the API call
      // This prevents the container refetch from overwriting local state
      const previousValue = container?.lookup?.[objectPath]?._value;
      setIntentForPath({
        jobId: job.id,
        parameterPath: objectPath,
        reason: "FileUpload",
        previousValue,
      });

      // Enqueue the operation to ensure sequential execution
      return parameterQueue.enqueue(async () => {
        try {

          const formData = new FormData();
          formData.append("objectPath", objectPath);
          formData.append("file", file, fileName);
          if (columnSelector?.trim()) {
            formData.append("column_selector", columnSelector);
          }
          // Enhanced multi-selector format for multiple representations
          if (columnSelectors && columnSelectors.length > 0) {
            formData.append("column_selectors", JSON.stringify(columnSelectors));
          }

          const result = await api.post<UploadFileParamResponse>(
            `jobs/${job.id}/upload_file_param`,
            formData
          );

          setProcessedErrors(null);

          // Update all related data
          await Promise.all([
            mutateValidation(),
            mutateContainer(),
            mutateParams_xml(),
          ]);

          // Invalidate the file digest cache for this objectPath
          // This triggers re-fetch of the digest, which task interfaces can use
          // to extract metadata like wavelength from the uploaded file
          // NOTE: Key must match useFileDigest exactly (no trailing slash)
          const digestKey = `jobs/${job.id}/digest?object_path=${objectPath}`;
          await mutate(digestKey);

          // NOTE: Don't clear intent here - let it expire naturally via auto-cleanup
          // The intent mechanism persists for a short window after changes

          console.log("File uploaded successfully:", result);
          return result;
        } catch (error) {
          // Clear intent on error so future syncs work
          clearIntentForPath(objectPath);
          console.error("Error uploading file:", error);
          throw error;
        }
      }) as Promise<UploadFileParamResponse | undefined>;
    },
    [
      job,
      container,
      mutateContainer,
      mutateValidation,
      mutateParams_xml,
      api,
      setProcessedErrors,
      setIntentForPath,
      clearIntentForPath,
    ]
  );

  const useTaskItem = useMemo(() => {
    return (paramName: string): TaskItem => {
      if (!paramName?.length || !container?.lookup) {
        return {
          item: null,
          value: null,
          update: async () => false,
          updateNoMutate: async () => false,
        };
      }

      const item = container.lookup[paramName];
      const value = valueOfItem(item);

      const update = async (newValue: any): Promise<boolean | Response> => {
        if (!job || job.status !== JOB_STATUS.PENDING) return false;

        // Check if value actually changed
        if (JSON.stringify({ value }) === JSON.stringify({ value: newValue })) {
          return false;
        }

        // Note: Intent is now recorded inside setParameter, no need to call setIntent here

        // Use the queued setParameter instead of direct fetch
        try {
          const result = await setParameter({
            object_path: item._objectPath,
            value: newValue,
          });

          return result?.success ? true : false;
        } catch (error) {
          console.error("Error updating task item:", error);
          return false;
        }
      };

      const updateNoMutate = async (
        newValue: any
      ): Promise<boolean | Response> => {
        if (!job || job.status !== JOB_STATUS.PENDING) return false;

        // Check if value actually changed
        if (JSON.stringify({ value }) === JSON.stringify({ value: newValue })) {
          return false;
        }

        // Use the queued setParameter instead of direct fetch
        try {
          const result = await setParameterNoMutate({
            object_path: item._objectPath,
            value: newValue,
          });

          return result?.success ? true : false;
        } catch (error) {
          console.error("Error updating task item:", error);
          return false;
        }
      };

      return { item, value, update, updateNoMutate };
    };
  }, [container, job, setParameter]);

  const createPeerTask = useCallback(
    async (taskName: string): Promise<Job | undefined> => {
      if (!job || !mutateJobs) {
        console.warn("Cannot create peer task: missing job or mutateJobs");
        return undefined;
      }

      try {
        console.log(`Creating ${taskName} task...`);

        const result = await api.post<CreateTaskResponse>(
          `projects/${job.project}/create_task/`,
          {
            task_name: taskName,
          }
        );

        if (result?.success && result.data?.new_job) {
          const createdJob: Job = result.data.new_job;
          await mutateJobs();
          return createdJob;
        }

        console.warn("Failed to create peer task:", result);
        return undefined;
      } catch (error) {
        console.error("Error creating peer task:", error);
        return undefined;
      }
    },
    [job, api, mutateJobs]
  );

  // Custom hook to fetch file content using SWR
  const useFileContent = (paramName: string): SWRResponse<string, Error> => {
    // Create a unique key for SWR caching
    const item = container?.lookup[paramName];
    const dbFileId = valueOfItem(item)?.dbFileId;

    //console.log("dbFileId", JSON.stringify(dbFileId));
    // Return null key when dbFileId is falsey - this prevents SWR from fetching
    const swrKey = dbFileId ? `files/${dbFileId}/download_by_uuid` : null;

    const fetcher = async (): Promise<string> => {
      if (!swrKey) {
        throw new Error(
          `Parameter "${paramName}" not found or has no dbFileId`
        );
      }
      return apiText(swrKey);
    };

    return useSWR<string, Error>(swrKey, swrKey ? fetcher : null, {
      // SWR options
      revalidateOnFocus: false,
      revalidateOnReconnect: true,
      errorRetryCount: 3,
      errorRetryInterval: 1000,
      // Cache for 5 minutes
      dedupingInterval: 5 * 60 * 1000,
      onError: (error) => {
        console.warn(
          `Error fetching file content for parameter "${paramName}":`,
          error
        );
        // Clear the cache for this key
        mutate(swrKey, null, false); // false = don't revalidate
      },
    });
  };

  // Custom hook to fetch file digest using SWR
  // Note: objectPath should be the full path like "prosmart_refmac.inputData.F_SIGF"
  // Returns unwrapped digest data (extracts .data from API response)
  const useFileDigest = (objectPath: string): SWRResponse<any, Error> => {
    // Create a unique key for SWR caching
    const swrKey = objectPath
      ? `jobs/${job?.id}/digest?object_path=${objectPath}`
      : null;
    const fetcher = async (): Promise<any> => {
      if (!swrKey) {
        throw new Error("Parameter not found");
      }
      // Unwrap API response format: {success: true, data: {...}}
      const result = await apiJson(swrKey);
      return result?.data ?? result;
    };

    return useSWR<any, Error>(swrKey, swrKey ? fetcher : null, {
      // SWR options
      revalidateOnFocus: false,
      revalidateOnReconnect: true,
      errorRetryCount: 3,
      errorRetryInterval: 1000,
      // Cache for 5 minutes
      dedupingInterval: 5 * 60 * 1000,
      onError: () => {
        // Clear the cache for this key
        mutate(swrKey, null, false); // false = don't revalidate
      },
    });
  };

  /**
   * Imperatively fetch the digest for a file parameter.
   * Use this in onChange callbacks for predictable, deterministic behavior.
   *
   * @param objectPath - Full object path like "ProvideAsuContents.inputData.ASUCONTENTIN"
   * @returns Promise resolving to digest data, or null if fetch fails
   *
   * @example
   * const handleFileChange = async () => {
   *   const digest = await fetchDigest(item._objectPath);
   *   if (digest?.seqList) {
   *     await setAsuContent(digest.seqList);
   *   }
   * };
   */
  const fetchDigest = useCallback(
    async (objectPath: string): Promise<any | null> => {
      if (!job?.id || !objectPath) return null;

      const digestKey = `jobs/${job.id}/digest?object_path=${objectPath}`;
      try {
        const result = await apiJson(digestKey);
        // Also update SWR cache so useFileDigest stays in sync
        mutate(digestKey, result, false);
        // Unwrap API response format: {success: true, data: {...}}
        return result?.data ?? result;
      } catch (error) {
        console.error(`Error fetching digest for ${objectPath}:`, error);
        return null;
      }
    },
    [job?.id]
  );

  const getValidationColor = useMemo(() => {
    return (item: any): string => {
      const fieldErrors = extractValidationErrors(
        item,
        processedErrors || validation || {}
      );
      return determineValidationColor(fieldErrors);
    };
  }, [processedErrors, validation]);

  const getErrors = useMemo(() => {
    return (item: any): ValidationError[] => {
      return extractValidationErrors(item, processedErrors || validation || {});
    };
  }, [processedErrors, validation]);

  const fileItemToParameterArg = useCallback(
    (
      value: DjangoFile,
      objectPath: string,
      projectJobs: Job[],
      projects: Project[]
    ): SetParameterArg => {
      // Base parameter structure
      const setParameterArg: SetParameterArg = {
        object_path: objectPath,
        value: {
          dbFileId: value.uuid.replace(/-/g, ""),
          subType: value.sub_type,
          contentFlag: value.content,
          annotation: value.annotation,
          baseName: value.name,
        },
      };

      // Handle different file directory types
      const handleFileDirectory = () => {
        if (value.directory === FILE_DIRECTORIES.IMPORTED) {
          setParameterArg.value.relPath = FILE_PATHS.IMPORTED_FILES;
        } else if (value.directory === FILE_DIRECTORIES.JOB_OUTPUT) {
          const jobOfFile = projectJobs?.find(
            (theJob) => theJob.id === value.job
          );
          if (jobOfFile) {
            const jobDir = jobOfFile.number
              .split(".")
              .map((element) => `${FILE_PATHS.JOB_PREFIX}${element}`)
              .join("/");

            setParameterArg.value.relPath = `${FILE_PATHS.JOBS_DIR}/${jobDir}`;
          }
        }

        // Set project for both directory types (moved outside the if/else)
        const jobOfFile = projectJobs?.find(
          (theJob) => theJob.id === value.job
        );
        if (jobOfFile) {
          const project = projects?.find(
            (theProject) => theProject.id === jobOfFile.project
          );
          if (project) {
            setParameterArg.value.project = project.uuid.replace(/-/g, "");
          }
        }
      };

      // Apply directory-specific handling only if job context exists
      if (job && projectJobs) {
        handleFileDirectory();
      }

      return setParameterArg;
    },
    [job]
  );

  return {
    job,
    mutateJob,
    container,
    mutateContainer,
    params_xml,
    mutateParams_xml,
    validation,
    mutateValidation: mutateValidationWithProcessedErrors, // use the decorated version
    diagnostic_xml,
    mutateDiagnosticXml,
    def_xml,
    mutateDef_xml,
    setParameter,
    setParameterNoMutate,
    uploadFileParam,
    useTaskItem,
    createPeerTask,
    useFileContent,
    getValidationColor,
    getErrors,
    useFileDigest,
    fetchDigest,
    fileItemToParameterArg,
  };
};

// ============================================================================
// Optional: Queue Status Hooks
// ============================================================================

/**
 * Hook to monitor parameter queue status
 */
export const useParameterQueueStatus = () => {
  const [queueLength, setQueueLength] = useState(0);
  const [isProcessing, setIsProcessing] = useState(false);

  useEffect(() => {
    const interval = setInterval(() => {
      setQueueLength(parameterQueue.getQueueLength());
      setIsProcessing(parameterQueue.isQueueProcessing());
    }, 100);

    return () => clearInterval(interval);
  }, []);

  return { queueLength, isProcessing };
};
