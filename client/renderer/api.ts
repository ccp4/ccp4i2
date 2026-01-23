import useSWR, { SWRConfiguration, SWRResponse } from "swr";
import $ from "jquery";
import { prettifyXml } from "./utils";
import {
  apiFetch,
  apiJson,
  apiText,
  apiPost,
  apiPatch,
  apiDelete,
} from "./api-fetch";
import { getAccessToken } from "./utils/auth-token";

// =============================================================================
// Constants
// =============================================================================

/** Standard polling interval for active jobs/dialogs */
export const POLL_INTERVAL = {
  ACTIVE: 5000,
  DISABLED: 0,
} as const;

// =============================================================================
// Types
// =============================================================================

/**
 * Standard API response format from the backend.
 * All endpoints now return this consistent structure.
 */
interface ApiResponse<T = any> {
  success: boolean;
  data?: T;
  error?: string;
  status: number;
}

export interface EndpointFetch {
  type: string;
  id: number | null | undefined;
  endpoint: string;
}

interface ValidationErrors {
  [objectPath: string]: {
    messages: string[];
    maxSeverity: number;
  };
}

// =============================================================================
// URL Helpers
// =============================================================================

/**
 * @deprecated No longer needed - apiGet/apiPost/etc. now auto-prefix with /api/proxy/ccp4i2/
 * Kept for backward compatibility. Just use the endpoint directly with api functions.
 */
export function makeApiUrl(endpoint: string): string {
  let api_path = `/api/proxy/ccp4i2/${endpoint}`;
  if (api_path.charAt(api_path.length - 1) !== "/") api_path += "/";
  return api_path;
}

function endpointToUrl(ef: EndpointFetch): string {
  return `${ef.type}/${ef.id}/${ef.endpoint}`;
}

function isValidEndpoint(ef: EndpointFetch | null | undefined): ef is EndpointFetch {
  // Django model IDs start at 1, so ID 0 is invalid
  return !!(ef?.type && ef.id !== null && ef.id !== undefined && ef.id > 0);
}

function isValidStringEndpoint(endpoint: string | null | undefined): endpoint is string {
  return !!(endpoint && !endpoint.includes("undefined") && !endpoint.includes("null"));
}

// =============================================================================
// Base Fetcher - Single source of truth for API calls
// =============================================================================

/**
 * Base fetcher that handles the standard API response format.
 * All other fetchers compose on top of this.
 */
async function baseFetcher<T>(url: string): Promise<T> {
  const response = await apiJson<ApiResponse<T>>(url);
  if (response?.success && response.data !== undefined) {
    return response.data;
  }
  throw new Error(response?.error || "API request failed");
}

/**
 * Simple JSON fetcher (no unwrapping - for endpoints that return data directly)
 */
const jsonFetcher = <T>(url: string): Promise<T> => apiJson<T>(url);

// =============================================================================
// Transformers - Pure functions that transform data
// =============================================================================

/**
 * Parse XML string to XMLDocument
 */
function parseXml(data: { xml: string } | null): XMLDocument | null {
  if (!data?.xml) return null;
  return $.parseXML(data.xml);
}

/**
 * Parse and prettify XML string
 */
function parseAndPrettifyXml(data: { xml: string } | null): string | null {
  if (!data?.xml) return null;
  return prettifyXml($.parseXML(data.xml));
}

/**
 * Transform validation XML to structured error object
 */
function parseValidationXml(data: { xml: string } | null): ValidationErrors {
  if (!data?.xml) return {};

  const validationXml = $.parseXML(data.xml);
  const objectPaths = $(validationXml).find("objectPath").toArray();
  const results: ValidationErrors = {};

  objectPaths.forEach((errorObjectNode: HTMLElement) => {
    const objectPath = errorObjectNode.textContent?.trim();
    if (!objectPath) return;

    if (!results[objectPath]) {
      results[objectPath] = { messages: [], maxSeverity: 0 };
    }

    const errorNode = $(errorObjectNode).parent();
    if (errorNode) {
      const severity = $(errorNode).find("severity").get(0)?.textContent;
      if (severity?.includes("WARNING") && results[objectPath].maxSeverity < 1) {
        results[objectPath].maxSeverity = 1;
      }
      if (severity?.includes("ERROR") && results[objectPath].maxSeverity < 2) {
        results[objectPath].maxSeverity = 2;
      }

      const description = $(errorNode).find("description").get(0)?.textContent;
      if (description) {
        results[objectPath].messages.push(description);
      }
    }
  });

  return results;
}

/**
 * Build lookup table for container navigation
 */
function buildLookup(container: any, lookup: Record<string, any> = {}): Record<string, any> {
  // Skip items without _objectPath (e.g., top-level container)
  const objectPath = container._objectPath;
  if (objectPath) {
    const pathElements = objectPath.split(".");

    // Add all suffix paths to lookup
    for (let i = 0; i < pathElements.length; i++) {
      const subPath = pathElements.slice(-i).join(".");
      lookup[subPath] = container;
    }
  }

  // Recurse into children
  if (container._baseClass === "CList" && Array.isArray(container._value)) {
    container._value.forEach((item: any) => {
      if (item && typeof item === "object") buildLookup(item, lookup);
    });
  } else if (container._value?.constructor === Object) {
    Object.values(container._value).forEach((item: any) => {
      if (item && typeof item === "object") buildLookup(item, lookup);
    });
  }

  return lookup;
}

/**
 * Transform container response (handles old and new API formats)
 */
function parseContainerResponse(response: any, isContainer: boolean): any {
  let result: any;

  // Handle new format: {"success": true, "data": {"result": ...}}
  if (response?.success && response.data?.result) {
    result = typeof response.data.result === "string"
      ? JSON.parse(response.data.result)
      : response.data.result;
  }
  // Handle old format: {"status": "Success", "result": ...}
  else if (response?.status === "Success" && response.result) {
    result = response.result;
  }
  else {
    throw new Error(response?.error || response?.reason || "Failed to fetch endpoint data");
  }

  // Build lookup for container endpoints
  if (isContainer) {
    return { container: result, lookup: buildLookup(result) };
  }

  return result;
}

// =============================================================================
// Composed Fetchers - Built from base fetcher + transformers
// =============================================================================

/**
 * Fetcher for endpoints returning {xml: string}
 */
function createXmlFetcher(transform: (data: { xml: string } | null) => any) {
  return async (ef: EndpointFetch) => {
    if (!isValidEndpoint(ef)) throw new Error("Invalid endpoint");
    const data = await baseFetcher<{ xml: string }>(endpointToUrl(ef));
    return transform(data);
  };
}

const xmlFetcher = createXmlFetcher(parseXml);
const prettyXmlFetcher = createXmlFetcher(parseAndPrettifyXml);
const validationFetcher = createXmlFetcher(parseValidationXml);

/**
 * Fetcher for wrapped JSON (container, etc.)
 */
async function wrappedJsonFetcher(ef: EndpointFetch): Promise<any> {
  if (!isValidEndpoint(ef)) throw new Error("Invalid endpoint");
  const response = await jsonFetcher<any>(endpointToUrl(ef));
  return parseContainerResponse(response, ef.endpoint === "container");
}

/**
 * Simple endpoint fetcher (returns raw API response)
 */
async function endpointFetcher<T>(ef: EndpointFetch): Promise<T> {
  if (!isValidEndpoint(ef)) throw new Error("Invalid endpoint");
  return jsonFetcher<T>(endpointToUrl(ef));
}

// =============================================================================
// SWR Key Helpers
// =============================================================================

function getEndpointKey(ef: EndpointFetch | null | undefined): EndpointFetch | null {
  return isValidEndpoint(ef) ? ef : null;
}

function getStringKey(endpoint: string | null | undefined): string | null {
  return isValidStringEndpoint(endpoint) ? endpoint : null;
}

// =============================================================================
// Main API Hook
// =============================================================================

export function useApi() {
  return {
    /**
     * Build URL without trailing slash
     */
    noSlashUrl(endpoint: string): string {
      return `/api/proxy/ccp4i2/${endpoint}`;
    },

    /**
     * Generic GET request with SWR caching
     * Pass null to skip fetching (conditional fetch pattern)
     */
    get<T>(endpoint: string | null, refreshInterval: number = 0) {
      return useSWR<T>(getStringKey(endpoint), jsonFetcher, {
        refreshInterval,
        dedupingInterval: 5000,  // Dedupe identical requests within 5 seconds
        keepPreviousData: true,  // Keep showing old data while revalidating
      });
    },

    /**
     * Fetch app config
     */
    config<T>() {
      return useSWR<T>("config", () => apiJson("/api/config"));
    },

    /**
     * Fetch from typed endpoint (returns raw response)
     * Pass null to skip fetching (conditional fetch pattern)
     */
    get_endpoint<T>(ef: EndpointFetch | null, refreshInterval: number = 0) {
      return useSWR<T>(getEndpointKey(ef), endpointFetcher as any, { refreshInterval });
    },

    /**
     * Fetch XML endpoint, parse to XMLDocument
     */
    get_endpoint_xml(ef: EndpointFetch, refreshInterval: number = 0) {
      return useSWR<XMLDocument | null>(getEndpointKey(ef), xmlFetcher, { refreshInterval });
    },

    /**
     * Fetch XML endpoint, parse and prettify
     * Pass null to skip fetching (conditional fetch pattern)
     */
    get_pretty_endpoint_xml(ef: EndpointFetch | null) {
      return useSWR<string | null>(getEndpointKey(ef), prettyXmlFetcher, {
        shouldRetryOnError: false,
        // Silently handle errors - caller can check error state
        onError: () => {},
      });
    },

    /**
     * Fetch JSON endpoint with legacy format handling (container, etc.)
     */
    get_wrapped_endpoint_json<T>(ef: EndpointFetch) {
      return useSWR<T>(getEndpointKey(ef), wrappedJsonFetcher);
    },

    /**
     * Fetch validation endpoint, transform to error map
     */
    get_validation(ef: EndpointFetch) {
      return useSWR<ValidationErrors>(getEndpointKey(ef), validationFetcher);
    },

    /**
     * Fetch digest with error handling
     */
    digest<T>(endpoint: string) {
      const swrConfig: SWRConfiguration = {
        onError: (error) => console.warn(`Digest error for "${endpoint}":`, error),
        fallbackData: null as T,
        shouldRetryOnError: false,
      };
      return useSWR<T>(getStringKey(endpoint), jsonFetcher, swrConfig);
    },

    /**
     * POST request
     */
    async post<T>(endpoint: string, body: any = {}): Promise<T> {
      return apiPost<T>(endpoint, body);
    },

    /**
     * DELETE request
     */
    async delete(endpoint: string): Promise<void> {
      await apiDelete(endpoint);
    },

    /**
     * PATCH request
     */
    async patch<T>(endpoint: string, body: any = {}): Promise<T> {
      return apiPatch<T>(endpoint, body);
    },

    /**
     * Fetch file text content by UUID
     */
    fileTextContent(djangoFile: any) {
      const swrKey = djangoFile?.dbFileId
        ? `files/${djangoFile.dbFileId}/download_by_uuid`
        : null;
      return useSWR<string>(swrKey, apiText);
    },

    // =========================================================================
    // Specialized hooks for common patterns (replaces direct useSWR calls)
    // =========================================================================

    /**
     * Fetch job report XML with optional polling.
     * Handles both wrapped {success, data: {xml}} and direct {xml} response formats.
     *
     * @param jobId - Job ID to fetch report for (null to skip)
     * @param shouldPoll - Whether to enable polling (e.g., when job is active)
     */
    jobReportXml(jobId: number | null | undefined, shouldPoll: boolean = false): SWRResponse<any> {
      const swrKey = jobId ? `jobs/${jobId}/report_xml` : null;
      return useSWR<any>(
        getStringKey(swrKey),
        jsonFetcher,
        { refreshInterval: shouldPoll ? POLL_INTERVAL.ACTIVE : POLL_INTERVAL.DISABLED }
      );
    },

    /**
     * Fetch project directory with optional polling.
     *
     * @param projectId - Project ID to fetch directory for (null to skip)
     * @param shouldPoll - Whether to enable polling (e.g., when dialog is open)
     */
    projectDirectory(projectId: number | null | undefined, shouldPoll: boolean = false): SWRResponse<any> {
      const swrKey = projectId ? `projects/${projectId}/directory` : null;
      return useSWR<any>(
        getStringKey(swrKey),
        jsonFetcher,
        { refreshInterval: shouldPoll ? POLL_INTERVAL.ACTIVE : POLL_INTERVAL.DISABLED }
      );
    },

    /**
     * Fetch a file by ID.
     *
     * @param fileId - File ID to fetch (null to skip)
     */
    file<T>(fileId: number | string | null | undefined): SWRResponse<T> {
      const swrKey = fileId ? `files/${fileId}` : null;
      return useSWR<T>(getStringKey(swrKey), jsonFetcher);
    },

    /**
     * Call an object method on a job with SWR caching.
     * Uses array keys for proper cache invalidation based on dependencies.
     *
     * @param jobId - Job ID
     * @param objectPath - Object path for the method call
     * @param methodName - Method name to call
     * @param kwargs - Optional keyword arguments for the method
     * @param deps - Optional dependencies that invalidate the cache when changed
     * @param enabled - Whether to enable fetching (for conditional fetches)
     */
    objectMethod<T>(
      jobId: number | null | undefined,
      objectPath: string,
      methodName: string,
      kwargs: Record<string, any> = {},
      deps: any[] = [],
      enabled: boolean = true
    ): SWRResponse<T> {
      // Build a cache key that includes all dependencies
      const swrKey = enabled && jobId
        ? [`jobs/${jobId}/object_method`, methodName, ...deps]
        : null;

      const fetcher = async ([url]: [string]) => {
        return apiPost<T>(url, {
          object_path: objectPath,
          method_name: methodName,
          kwargs,
        });
      };

      return useSWR<T>(swrKey, fetcher, { keepPreviousData: true });
    },
  };
}

// =============================================================================
// File Download Utilities
// =============================================================================

/**
 * Download a file using the browser's native download capability.
 * This streams directly to disk without loading into memory - safe for large files.
 *
 * For authenticated endpoints, the access token is appended as a query parameter
 * since anchor clicks don't include Authorization headers.
 *
 * Note: This doesn't provide progress feedback, but avoids memory issues with large files.
 */
export const doDownload = async (
  theURL: string,
  targetName: string,
  _optionsIn?: any,
  _onProgress?: (bytesRead: number) => void
) => {
  // Get access token to append to URL (anchor clicks don't send Authorization header)
  const token = await getAccessToken();

  // Append token as query parameter if available
  let downloadUrl = theURL;
  if (token) {
    const separator = theURL.includes("?") ? "&" : "?";
    downloadUrl = `${theURL}${separator}access_token=${encodeURIComponent(token)}`;
  }

  // Use direct anchor approach - browser handles streaming to disk
  // This avoids loading the entire file into memory
  const link = document.createElement("a");
  link.href = downloadUrl;
  link.download = targetName;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
};

/**
 * Download a file with progress tracking.
 * WARNING: This loads the entire file into browser memory before downloading.
 * Only use for small files where progress feedback is important.
 */
export const doDownloadWithProgress = (
  theURL: string,
  targetName: string,
  optionsIn?: any,
  onProgress: (bytesRead: number) => void = (bytesRead) => console.log(bytesRead)
) => {
  const options = optionsIn ?? {};

  if (!onProgress) return;

  return apiFetch(theURL, options).then(async (response) => {
    const reader = response.body?.getReader();
    if (!reader) return;

    const chunks: Uint8Array[] = [];
    let receivedLength = 0;

    while (true) {
      const { done, value } = await reader.read();
      if (value) {
        chunks.push(value);
        receivedLength += value.length;
        onProgress(receivedLength);
      }
      if (done) break;
    }

    // Combine chunks
    const combined = new Uint8Array(receivedLength);
    let position = 0;
    for (const chunk of chunks) {
      combined.set(chunk, position);
      position += chunk.length;
    }

    // Trigger download
    const blob = new Blob([combined]);
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.setAttribute("download", targetName);
    document.body.appendChild(link);
    link.click();
    link.parentNode?.removeChild(link);
    URL.revokeObjectURL(url);
  });
};

export const doRetrieve = async (
  theURL: string,
  _targetName: string,
  optionsIn?: any
): Promise<ArrayBuffer> => {
  const options = optionsIn ?? {};
  const response = await apiFetch(theURL, options);
  return response.arrayBuffer();
};
