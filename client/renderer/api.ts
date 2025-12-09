import useSWR, { SWRConfiguration } from "swr";
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
 * @deprecated No longer needed - apiGet/apiPost/etc. now auto-prefix with /api/proxy/
 * Kept for backward compatibility. Just use the endpoint directly with api functions.
 */
export function makeApiUrl(endpoint: string): string {
  let api_path = `/api/proxy/${endpoint}`;
  if (api_path.charAt(api_path.length - 1) !== "/") api_path += "/";
  return api_path;
}

function endpointToUrl(ef: EndpointFetch): string {
  return `${ef.type}/${ef.id}/${ef.endpoint}`;
}

function isValidEndpoint(ef: EndpointFetch | null | undefined): ef is EndpointFetch {
  return !!(ef?.id && ef?.type);
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
      return `/api/proxy/${endpoint}`;
    },

    /**
     * Generic GET request with SWR caching
     * Pass null to skip fetching (conditional fetch pattern)
     */
    get<T>(endpoint: string | null, refreshInterval: number = 0) {
      return useSWR<T>(getStringKey(endpoint), jsonFetcher, { refreshInterval });
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
  };
}

// =============================================================================
// File Download Utilities
// =============================================================================

export const doDownload = (
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
