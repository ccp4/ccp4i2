/**
 * Centralized API fetching utilities — shared library factory.
 *
 * Consumers call ``createApiFetch({ baseUrl })`` to bind a set of fetch
 * helpers (apiFetch, apiJson, apiPost, apiPut, apiPatch, apiDelete,
 * apiGet, apiUpload, apiText, apiBlob, apiArrayBuffer, swrFetcher,
 * swrPostFetcher) to a specific URL prefix.
 *
 * The generic implementation handles:
 * - Authentication token injection via ``getAccessToken()`` from the
 *   companion auth-token module (set up via ``setTokenGetter``).
 * - URL normalisation: relative paths get prefixed with ``baseUrl``;
 *   ``/api/...`` paths and absolute URLs pass through unchanged.
 * - 401/403 detection that emits ``AUTH_ERROR_EVENT`` on ``window`` so
 *   UI components can drive re-auth UX.
 * - Body-type-aware Content-Type negotiation (JSON vs form-encoded vs
 *   FormData).
 * - Optional request timeout.
 *
 * Each consumer typically has a thin local file like
 * ``client/renderer/api-fetch.ts`` that calls ``createApiFetch`` once
 * and re-exports the bound helpers, so call sites can keep their
 * existing relative imports.
 */

import { getAccessToken } from "./auth-token.js";

// =============================================================================
// Universal exports — auth-error event surface
// =============================================================================

/**
 * Event name dispatched on ``window`` when an API call receives a 401
 * or 403. UI components can listen and show re-auth prompts.
 */
export const AUTH_ERROR_EVENT = "ccp4i2:auth-error";

export interface AuthErrorDetail {
  status: 401 | 403;
  url: string;
  message: string;
}

function emitAuthError(detail: AuthErrorDetail): void {
  if (typeof window !== "undefined") {
    window.dispatchEvent(new CustomEvent(AUTH_ERROR_EVENT, { detail }));
  }
}

// =============================================================================
// Per-request configuration + factory options
// =============================================================================

export interface ApiFetchConfig {
  headers?: Record<string, string>;
  timeout?: number;
  retries?: number;
}

const DEFAULT_CONFIG: ApiFetchConfig = {
  timeout: 30000,
  retries: 1,
};

export interface CreateApiFetchOptions {
  /**
   * URL prefix for relative API paths, e.g. ``"/api/proxy/ccp4i2/"`` or
   * ``"/api/proxy/compounds/"``. Must include trailing slash. URLs that
   * already start with ``/api/`` or ``http(s)://`` pass through
   * unchanged. The proxy is responsible for adding any trailing slash
   * Django requires; this layer never adds one.
   */
  baseUrl: string;
}

/**
 * The set of fetch helpers a consumer gets from ``createApiFetch``. All
 * helpers share the same baseUrl and authentication chain.
 */
export interface ApiFetcher {
  apiFetch(
    url: string,
    options?: RequestInit,
    config?: ApiFetchConfig,
  ): Promise<Response>;
  apiJson<T = any>(
    url: string,
    options?: RequestInit,
    config?: ApiFetchConfig,
  ): Promise<T>;
  apiText(
    url: string,
    options?: RequestInit,
    config?: ApiFetchConfig,
  ): Promise<string>;
  apiBlob(
    url: string,
    options?: RequestInit,
    config?: ApiFetchConfig,
  ): Promise<Blob>;
  apiArrayBuffer(
    url: string,
    options?: RequestInit,
    config?: ApiFetchConfig,
  ): Promise<ArrayBuffer>;
  apiPost<T = any>(
    url: string,
    data: any,
    config?: ApiFetchConfig,
  ): Promise<T>;
  apiPut<T = any>(
    url: string,
    data: any,
    config?: ApiFetchConfig,
  ): Promise<T>;
  apiPatch<T = any>(
    url: string,
    data: any,
    config?: ApiFetchConfig,
  ): Promise<T>;
  apiDelete<T = any>(url: string, config?: ApiFetchConfig): Promise<T>;
  apiGet<T = any>(url: string, config?: ApiFetchConfig): Promise<T>;
  apiUpload<T = any>(
    url: string,
    file: File | FormData,
    config?: ApiFetchConfig,
  ): Promise<T>;
  swrFetcher(url: string): Promise<any>;
  swrPostFetcher(data: any): (url: string) => Promise<any>;
}

/**
 * Build a set of fetch helpers bound to a given URL prefix and the
 * shared auth-token chain. Consumers call this once at module load.
 */
export function createApiFetch(options: CreateApiFetchOptions): ApiFetcher {
  const { baseUrl } = options;

  function normalizeApiUrl(url: string): string {
    // Absolute URLs pass through unchanged.
    if (url.startsWith("http://") || url.startsWith("https://")) {
      return url;
    }
    // Already has /api/ prefix (covers /api/proxy/, /api/config, etc.).
    if (url.startsWith("/api/")) {
      return url;
    }
    // Strip leading slash if present for consistency.
    const endpoint = url.startsWith("/") ? url.slice(1) : url;
    return `${baseUrl}${endpoint}`;
  }

  /**
   * Core fetch wrapper with error handling and optional timeout.
   * Token injection happens here; URL normalisation happens here;
   * everything else is HTTP plumbing.
   */
  async function coreFetch(
    url: string,
    requestOptions: RequestInit = {},
    config: ApiFetchConfig = {},
  ): Promise<Response> {
    const finalConfig = { ...DEFAULT_CONFIG, ...config };
    const normalizedUrl = normalizeApiUrl(url);

    const headers: Record<string, string> = {
      ...finalConfig.headers,
    };

    const token = await getAccessToken();
    if (token) {
      headers["Authorization"] = `Bearer ${token}`;
      if (
        process.env.NODE_ENV === "development" &&
        process.env.DEBUG_AUTH
      ) {
        console.log(
          "[FETCH] Token injected for:",
          normalizedUrl.substring(0, 60),
        );
      }
    }

    if (requestOptions.headers) {
      if (requestOptions.headers instanceof Headers) {
        requestOptions.headers.forEach((value, key) => {
          headers[key] = value;
        });
      } else if (Array.isArray(requestOptions.headers)) {
        requestOptions.headers.forEach(([key, value]) => {
          headers[key] = value;
        });
      } else {
        Object.assign(headers, requestOptions.headers);
      }
    }

    // Set JSON content-type only when:
    //   1. None already set
    //   2. Body isn't FormData (browser sets multipart with boundary)
    //   3. Body exists and is a string-ish payload
    const hasContentType = Object.keys(headers).some(
      (key) => key.toLowerCase() === "content-type",
    );
    if (!hasContentType && requestOptions.body) {
      if (!(requestOptions.body instanceof FormData)) {
        if (
          typeof requestOptions.body === "string" ||
          requestOptions.body instanceof URLSearchParams
        ) {
          headers["Content-Type"] =
            typeof requestOptions.body === "string"
              ? "application/json"
              : "application/x-www-form-urlencoded";
        }
      }
    }

    const finalRequestOptions: RequestInit = {
      ...requestOptions,
      headers,
    };

    if (finalConfig.timeout) {
      const controller = new AbortController();
      const timeoutId = setTimeout(
        () => controller.abort(),
        finalConfig.timeout,
      );
      finalRequestOptions.signal = controller.signal;
      try {
        const response = await fetch(normalizedUrl, finalRequestOptions);
        clearTimeout(timeoutId);
        return response;
      } catch (error) {
        clearTimeout(timeoutId);
        throw error;
      }
    }

    return fetch(normalizedUrl, finalRequestOptions);
  }

  /** Returns the Response object. Use when you need headers / status. */
  async function apiFetch(
    url: string,
    options: RequestInit = {},
    config: ApiFetchConfig = {},
  ): Promise<Response> {
    try {
      const response = await coreFetch(url, options, config);

      if (!response.ok) {
        let errorMessage = `HTTP ${response.status}: ${response.statusText}`;
        try {
          const errorData = await response.json();
          if (errorData?.error) {
            errorMessage = errorData.error;
          }
        } catch {
          // Body wasn't JSON; default message stands.
        }

        if (response.status === 401 || response.status === 403) {
          emitAuthError({
            status: response.status as 401 | 403,
            url,
            message:
              response.status === 401
                ? "Your session has expired. Please sign in again."
                : "You don't have permission to access this resource.",
          });
        }

        throw new Error(errorMessage);
      }

      return response;
    } catch (error) {
      console.error(`API fetch failed for ${url}:`, error);
      throw error;
    }
  }

  async function apiJson<T = any>(
    url: string,
    options: RequestInit = {},
    config: ApiFetchConfig = {},
  ): Promise<T> {
    const response = await apiFetch(url, options, config);
    return response.json();
  }

  async function apiText(
    url: string,
    options: RequestInit = {},
    config: ApiFetchConfig = {},
  ): Promise<string> {
    const response = await apiFetch(url, options, config);
    return response.text();
  }

  async function apiBlob(
    url: string,
    options: RequestInit = {},
    config: ApiFetchConfig = {},
  ): Promise<Blob> {
    const response = await apiFetch(url, options, config);
    return response.blob();
  }

  async function apiArrayBuffer(
    url: string,
    options: RequestInit = {},
    config: ApiFetchConfig = {},
  ): Promise<ArrayBuffer> {
    const response = await apiFetch(url, options, config);
    return response.arrayBuffer();
  }

  function bodyForMethod(data: any): BodyInit {
    if (
      data instanceof FormData ||
      typeof data === "string" ||
      data instanceof URLSearchParams
    ) {
      return data as BodyInit;
    }
    return JSON.stringify(data);
  }

  async function apiPost<T = any>(
    url: string,
    data: any,
    config: ApiFetchConfig = {},
  ): Promise<T> {
    return apiJson<T>(
      url,
      { method: "POST", body: bodyForMethod(data) },
      config,
    );
  }

  async function apiPut<T = any>(
    url: string,
    data: any,
    config: ApiFetchConfig = {},
  ): Promise<T> {
    return apiJson<T>(
      url,
      { method: "PUT", body: bodyForMethod(data) },
      config,
    );
  }

  async function apiPatch<T = any>(
    url: string,
    data: any,
    config: ApiFetchConfig = {},
  ): Promise<T> {
    return apiJson<T>(
      url,
      { method: "PATCH", body: bodyForMethod(data) },
      config,
    );
  }

  async function apiDelete<T = any>(
    url: string,
    config: ApiFetchConfig = {},
  ): Promise<T> {
    return apiJson<T>(url, { method: "DELETE" }, config);
  }

  async function apiGet<T = any>(
    url: string,
    config: ApiFetchConfig = {},
  ): Promise<T> {
    return apiJson<T>(url, {}, config);
  }

  async function apiUpload<T = any>(
    url: string,
    file: File | FormData,
    config: ApiFetchConfig = {},
  ): Promise<T> {
    const formData = file instanceof FormData ? file : new FormData();
    if (file instanceof File) {
      formData.append("file", file);
    }

    const response = await apiFetch(
      url,
      { method: "POST", body: formData },
      config,
    );

    const contentType = response.headers.get("content-type");
    if (contentType && contentType.includes("application/json")) {
      return response.json();
    }
    return (await response.text()) as unknown as T;
  }

  const swrFetcher = (url: string) => apiJson(url);
  const swrPostFetcher = (data: any) => (url: string) =>
    apiPost(url, data);

  return {
    apiFetch,
    apiJson,
    apiText,
    apiBlob,
    apiArrayBuffer,
    apiPost,
    apiPut,
    apiPatch,
    apiDelete,
    apiGet,
    apiUpload,
    swrFetcher,
    swrPostFetcher,
  };
}
