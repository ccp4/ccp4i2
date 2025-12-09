/**
 * Centralized API fetching utilities
 * All fetch() calls should go through these functions to enable:
 * - Consistent error handling
 * - Authentication token injection
 * - Request/response logging
 * - Retry logic
 * - Automatic /api/proxy/ prefix for relative API endpoints
 */

/**
 * Configuration for API requests
 */
interface ApiFetchConfig {
  headers?: Record<string, string>;
  timeout?: number;
  retries?: number;
}

/**
 * Default configuration
 */
const DEFAULT_CONFIG: ApiFetchConfig = {
  timeout: 30000,
  retries: 1,
};

/**
 * Normalize a URL to include the /api/proxy/ prefix for API endpoints.
 *
 * - If URL starts with "http://" or "https://", return as-is (absolute URL)
 * - If URL starts with "/api/", return as-is (already prefixed)
 * - Otherwise, prepend "/api/proxy/"
 *
 * NOTE: Trailing slashes are NOT added here. The proxy route (route.ts) handles
 * adding trailing slashes when forwarding to Django. This avoids 308 redirects
 * from Next.js when trailingSlash config doesn't match the URL.
 *
 * @example
 * normalizeApiUrl("jobs/123/digest") => "/api/proxy/jobs/123/digest"
 * normalizeApiUrl("jobs/123/digest?foo=bar") => "/api/proxy/jobs/123/digest?foo=bar"
 * normalizeApiUrl("/api/proxy/jobs/123") => "/api/proxy/jobs/123"
 * normalizeApiUrl("/api/config") => "/api/config"
 * normalizeApiUrl("https://example.com/data") => "https://example.com/data"
 */
function normalizeApiUrl(url: string): string {
  // Absolute URLs pass through unchanged
  if (url.startsWith("http://") || url.startsWith("https://")) {
    return url;
  }

  // Already has /api/ prefix (covers /api/proxy/, /api/config, etc.)
  if (url.startsWith("/api/")) {
    return url;
  }

  // Strip leading slash if present for consistency
  const endpoint = url.startsWith("/") ? url.slice(1) : url;

  return `/api/proxy/${endpoint}`;
}

/**
 * Core fetch wrapper with error handling and timeout
 */
async function coreFetch(
  url: string,
  options: RequestInit = {},
  config: ApiFetchConfig = {}
): Promise<Response> {
  const finalConfig = { ...DEFAULT_CONFIG, ...config };

  // Normalize URL to include /api/proxy/ prefix for API endpoints
  const normalizedUrl = normalizeApiUrl(url);

  // Prepare headers - don't assume JSON content type
  const headers: Record<string, string> = {
    ...finalConfig.headers,
  };

  // Convert options.headers to Record<string, string> format
  if (options.headers) {
    if (options.headers instanceof Headers) {
      // Convert Headers object to Record<string, string>
      options.headers.forEach((value, key) => {
        headers[key] = value;
      });
    } else if (Array.isArray(options.headers)) {
      // Convert string[][] to Record<string, string>
      options.headers.forEach(([key, value]) => {
        headers[key] = value;
      });
    } else {
      // Already Record<string, string>, safe to spread
      Object.assign(headers, options.headers);
    }
  }

  // Only set JSON content-type if:
  // 1. No Content-Type is already set
  // 2. Body is not FormData (which needs multipart/form-data with boundary)
  // 3. Body exists and appears to be JSON
  const hasContentType = Object.keys(headers).some(
    (key) => key.toLowerCase() === "content-type"
  );

  if (!hasContentType && options.body) {
    // Don't set Content-Type for FormData - browser will set multipart/form-data with boundary
    if (!(options.body instanceof FormData)) {
      // Only set JSON content-type for non-FormData bodies
      if (
        typeof options.body === "string" ||
        options.body instanceof URLSearchParams
      ) {
        // Assume JSON for string bodies, form-encoded for URLSearchParams
        headers["Content-Type"] =
          typeof options.body === "string"
            ? "application/json"
            : "application/x-www-form-urlencoded";
      }
    }
  }

  // Set up request options
  const requestOptions: RequestInit = {
    ...options,
    headers,
  };

  // Add timeout if specified
  if (finalConfig.timeout) {
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), finalConfig.timeout);
    requestOptions.signal = controller.signal;

    try {
      const response = await fetch(normalizedUrl, requestOptions);
      clearTimeout(timeoutId);
      return response;
    } catch (error) {
      clearTimeout(timeoutId);
      throw error;
    }
  }

  return fetch(normalizedUrl, requestOptions);
}

/**
 * Generic fetch that returns the Response object
 * Use this when you need access to headers, status codes, etc.
 */
export async function apiFetch(
  url: string,
  options: RequestInit = {},
  config: ApiFetchConfig = {}
): Promise<Response> {
  try {
    const response = await coreFetch(url, options, config);

    if (!response.ok) {
      throw new Error(`HTTP ${response.status}: ${response.statusText}`);
    }

    return response;
  } catch (error) {
    console.error(`API fetch failed for ${url}:`, error);
    throw error;
  }
}

/**
 * Fetch and return JSON data
 * Most common use case for API calls
 */
export async function apiJson<T = any>(
  url: string,
  options: RequestInit = {},
  config: ApiFetchConfig = {}
): Promise<T> {
  const response = await apiFetch(url, options, config);
  return response.json();
}

/**
 * Fetch and return text data
 * Useful for XML, CSV, or other text-based responses
 */
export async function apiText(
  url: string,
  options: RequestInit = {},
  config: ApiFetchConfig = {}
): Promise<string> {
  const response = await apiFetch(url, options, config);
  return response.text();
}

/**
 * Fetch and return blob data
 * Useful for file downloads
 */
export async function apiBlob(
  url: string,
  options: RequestInit = {},
  config: ApiFetchConfig = {}
): Promise<Blob> {
  const response = await apiFetch(url, options, config);
  return response.blob();
}

/**
 * Fetch and return array buffer data
 * Useful for binary file processing
 */
export async function apiArrayBuffer(
  url: string,
  options: RequestInit = {},
  config: ApiFetchConfig = {}
): Promise<ArrayBuffer> {
  const response = await apiFetch(url, options, config);
  return response.arrayBuffer();
}

/**
 * POST request with automatic body type handling
 * Supports JSON objects, FormData, strings, etc.
 * Content-Type is automatically set based on body type
 */
export async function apiPost<T = any>(
  url: string,
  data: any,
  config: ApiFetchConfig = {}
): Promise<T> {
  let body: any;

  // Let coreFetch handle FormData and other types automatically
  if (
    data instanceof FormData ||
    typeof data === "string" ||
    data instanceof URLSearchParams
  ) {
    body = data;
  } else {
    // For plain objects, JSON stringify them
    body = JSON.stringify(data);
  }

  return apiJson<T>(
    url,
    {
      method: "POST",
      body,
      // Don't set Content-Type here - let coreFetch decide based on body type
    },
    config
  );
}

/**
 * PUT request with automatic body type handling
 * Supports JSON objects, FormData, strings, etc.
 * Content-Type is automatically set based on body type
 */
export async function apiPut<T = any>(
  url: string,
  data: any,
  config: ApiFetchConfig = {}
): Promise<T> {
  let body: any;

  // Let coreFetch handle FormData and other types automatically
  if (
    data instanceof FormData ||
    typeof data === "string" ||
    data instanceof URLSearchParams
  ) {
    body = data;
  } else {
    // For plain objects, JSON stringify them
    body = JSON.stringify(data);
  }

  return apiJson<T>(
    url,
    {
      method: "PUT",
      body,
      // Don't set Content-Type here - let coreFetch decide based on body type
    },
    config
  );
}

/**
 * PATCH request with automatic body type handling
 * Supports JSON objects, FormData, strings, etc.
 * Content-Type is automatically set based on body type
 */
export async function apiPatch<T = any>(
  url: string,
  data: any,
  config: ApiFetchConfig = {}
): Promise<T> {
  let body: any;

  // Let coreFetch handle FormData and other types automatically
  if (
    data instanceof FormData ||
    typeof data === "string" ||
    data instanceof URLSearchParams
  ) {
    body = data;
  } else {
    // For plain objects, JSON stringify them
    body = JSON.stringify(data);
  }

  return apiJson<T>(
    url,
    {
      method: "PATCH",
      body,
      // Don't set Content-Type here - let coreFetch decide based on body type
    },
    config
  );
}

/**
 * DELETE request
 */
export async function apiDelete<T = any>(
  url: string,
  config: ApiFetchConfig = {}
): Promise<T> {
  return apiJson<T>(
    url,
    {
      method: "DELETE",
    },
    config
  );
}

/**
 * GET request for JSON data (convenience method)
 */
export async function apiGet<T = any>(
  url: string,
  config: ApiFetchConfig = {}
): Promise<T> {
  return apiJson<T>(url, {}, config);
}

/**
 * Upload file with multipart/form-data
 */
export async function apiUpload<T = any>(
  url: string,
  file: File | FormData,
  config: ApiFetchConfig = {}
): Promise<T> {
  const formData = file instanceof FormData ? file : new FormData();
  if (file instanceof File) {
    formData.append("file", file);
  }

  // Use apiFetch directly to avoid any JSON content-type assumptions
  const response = await apiFetch(
    url,
    {
      method: "POST",
      body: formData,
      // No Content-Type header - browser will set multipart/form-data with boundary
    },
    config
  );

  // Try to parse as JSON, fallback to text
  const contentType = response.headers.get("content-type");
  if (contentType && contentType.includes("application/json")) {
    return response.json();
  } else {
    return response.text() as any;
  }
}

/**
 * Helper for SWR fetcher functions
 * Usage: const { data, error } = useSWR('/api/endpoint', swrFetcher);
 */
export const swrFetcher = (url: string) => apiJson(url);

/**
 * Helper for creating SWR fetcher with POST data
 */
export const swrPostFetcher = (data: any) => (url: string) =>
  apiPost(url, data);
