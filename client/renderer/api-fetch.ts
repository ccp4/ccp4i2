/**
 * CCP4i2-specific binding of the shared api-fetch wrapper.
 *
 * The generic implementation lives in @ccp4/ccp4i2-api's
 * ``createApiFetch`` factory; this file binds it to the
 * ``/api/proxy/ccp4i2/`` proxy prefix used by the desktop renderer and
 * the cloud web app, and re-exports the bound helpers under their
 * legacy names so existing import sites
 * (``import { apiFetch, apiJson, ... } from "../api-fetch"``) keep
 * working unchanged.
 *
 * AUTH_ERROR_EVENT and AuthErrorDetail are universal (the event name is
 * the same across consumers) and so come straight from the package.
 */

import {
  AUTH_ERROR_EVENT,
  type ApiFetchConfig,
  type AuthErrorDetail,
  createApiFetch,
  getAccessToken,
} from "@ccp4/ccp4i2-api";

export { AUTH_ERROR_EVENT };
export type { AuthErrorDetail, ApiFetchConfig };

/**
 * Proxy prefix for every relative CCP4i2 API path. Shared by the
 * factory below and by ``apiUploadWithProgress``, which resolves URLs
 * itself because it does not go through the factory's fetch path.
 */
const PROXY_BASE = "/api/proxy/ccp4i2/";

const fetcher = createApiFetch({ baseUrl: PROXY_BASE });

export const apiFetch = fetcher.apiFetch;
export const apiJson = fetcher.apiJson;
export const apiText = fetcher.apiText;
export const apiBlob = fetcher.apiBlob;
export const apiArrayBuffer = fetcher.apiArrayBuffer;
export const apiPost = fetcher.apiPost;
export const apiPut = fetcher.apiPut;
export const apiPatch = fetcher.apiPatch;
export const apiDelete = fetcher.apiDelete;
export const apiGet = fetcher.apiGet;
export const swrFetcher = fetcher.swrFetcher;
export const swrPostFetcher = fetcher.swrPostFetcher;

// =============================================================================
// Uploads
// =============================================================================

/**
 * Upload a file, with no wall-clock ceiling by default.
 *
 * ``createApiFetch`` puts a 30 000 ms AbortController around every
 * request. That is a sane ceiling for a JSON round trip and a fatal one
 * for a file upload: the clock starts when the request is issued and
 * keeps running while the bytes are still going up the wire, so any file
 * that takes longer than 30 s to transfer is aborted by the browser
 * itself, part-way through the body. Nothing reaches Django — no request
 * is logged there, no status comes back — the socket simply resets and
 * the caller gets an "aborted" error with no indication that size was
 * the problem. An unmerged Diamond MTZ is routinely 50-200 MB, well past
 * 30 s on any normal uplink, so this made batch import of real fragment
 * data impossible.
 *
 * The meaningful bound on an upload is *progress*, not elapsed time, so
 * the default here is no timeout; callers wanting a ceiling pass one in
 * ``config``. Prefer ``apiUploadWithProgress`` anywhere a person is
 * watching: it reports bytes sent and aborts on a genuine stall.
 */
export const apiUpload = <T = any>(
  url: string,
  file: File | FormData,
  config: ApiFetchConfig = {},
): Promise<T> => fetcher.apiUpload<T>(url, file, { timeout: 0, ...config });

export interface UploadProgress {
  /** Bytes sent so far. */
  loaded: number;
  /** Total bytes to send, or 0 when the browser can't say. */
  total: number;
  /** 0-1, or null when ``total`` is unknown. */
  fraction: number | null;
}

export interface ApiUploadWithProgressOptions {
  /** Called as bytes go up the wire. */
  onProgress?: (progress: UploadProgress) => void;
  /** Caller-controlled cancellation. */
  signal?: AbortSignal;
  /**
   * Abort if no bytes move for this many ms. This is the honest form of
   * an upload timeout — it fires on a dead connection but never on a
   * merely large file. 0 disables it. Default 120 s.
   */
  stallTimeout?: number;
}

/** Mirrors the factory's URL rule: absolute and /api/... pass through. */
function resolveProxyUrl(url: string): string {
  if (url.startsWith("http://") || url.startsWith("https://")) return url;
  if (url.startsWith("/api/")) return url;
  return `${PROXY_BASE}${url.startsWith("/") ? url.slice(1) : url}`;
}

function formatBytes(bytes: number): string {
  if (!Number.isFinite(bytes) || bytes < 0) return "?";
  if (bytes < 1024) return `${bytes} B`;
  const units = ["kB", "MB", "GB"];
  let value = bytes / 1024;
  let unit = 0;
  while (value >= 1024 && unit < units.length - 1) {
    value /= 1024;
    unit += 1;
  }
  return `${value.toFixed(value < 10 ? 1 : 0)} ${units[unit]}`;
}

function emitAuthError(detail: AuthErrorDetail): void {
  if (typeof window !== "undefined") {
    window.dispatchEvent(new CustomEvent(AUTH_ERROR_EVENT, { detail }));
  }
}

/**
 * Upload a file and report progress as it goes.
 *
 * Uses XMLHttpRequest rather than fetch because fetch still gives no
 * upload-progress events in any shipping browser, and a multi-hundred-
 * megabyte upload with no visible progress is indistinguishable from a
 * hang. Auth, URL resolution and the 401/403 auth-error event match what
 * ``createApiFetch`` does, so this is a drop-in for ``apiUpload``.
 *
 * Errors carry what actually happened, including how far the upload got,
 * so a UI can say something better than "aborted".
 */
export function apiUploadWithProgress<T = any>(
  url: string,
  file: File | FormData,
  options: ApiUploadWithProgressOptions = {},
): Promise<T> {
  const { onProgress, signal, stallTimeout = 120000 } = options;

  const formData =
    file instanceof FormData
      ? file
      : (() => {
          const fd = new FormData();
          fd.append("file", file);
          return fd;
        })();

  return new Promise<T>((resolve, reject) => {
    if (signal?.aborted) {
      reject(new DOMException("Upload cancelled", "AbortError"));
      return;
    }

    getAccessToken()
      .then((token) => {
        const xhr = new XMLHttpRequest();
        const resolvedUrl = resolveProxyUrl(url);

        let loaded = 0;
        let total = 0;
        let lastMovement = Date.now();
        let stallTimer: ReturnType<typeof setInterval> | undefined;
        let settled = false;

        const cleanup = () => {
          if (stallTimer !== undefined) clearInterval(stallTimer);
          signal?.removeEventListener("abort", onAbortRequested);
        };
        const fail = (error: Error) => {
          if (settled) return;
          settled = true;
          cleanup();
          console.error(`API upload failed for ${url}:`, error);
          reject(error);
        };
        const succeed = (value: T) => {
          if (settled) return;
          settled = true;
          cleanup();
          resolve(value);
        };

        function onAbortRequested() {
          xhr.abort();
          fail(new DOMException("Upload cancelled", "AbortError"));
        }
        signal?.addEventListener("abort", onAbortRequested);

        xhr.upload.addEventListener("progress", (event) => {
          loaded = event.loaded;
          total = event.lengthComputable ? event.total : 0;
          lastMovement = Date.now();
          onProgress?.({
            loaded,
            total,
            fraction: total > 0 ? loaded / total : null,
          });
        });

        if (stallTimeout > 0) {
          stallTimer = setInterval(() => {
            if (Date.now() - lastMovement < stallTimeout) return;
            xhr.abort();
            fail(
              new Error(
                `Upload stalled: no data sent for ${Math.round(
                  stallTimeout / 1000,
                )}s after ${formatBytes(loaded)}${
                  total ? ` of ${formatBytes(total)}` : ""
                }. The connection to the server was lost.`,
              ),
            );
          }, 5000);
        }

        xhr.addEventListener("error", () => {
          fail(
            new Error(
              `Network error during upload after ${formatBytes(loaded)}${
                total ? ` of ${formatBytes(total)}` : ""
              }. The connection closed before the server replied.`,
            ),
          );
        });

        xhr.addEventListener("timeout", () => {
          fail(new Error("Upload timed out."));
        });

        xhr.addEventListener("load", () => {
          const contentType = xhr.getResponseHeader("content-type") || "";
          const isJson = contentType.includes("application/json");

          if (xhr.status >= 200 && xhr.status < 300) {
            if (!isJson) {
              succeed(xhr.responseText as unknown as T);
              return;
            }
            try {
              succeed(JSON.parse(xhr.responseText) as T);
            } catch {
              fail(
                new Error(
                  "Server returned a malformed JSON response to the upload.",
                ),
              );
            }
            return;
          }

          // Prefer the backend's own message; it is far more useful than
          // the status line when a task rejects a file.
          let message = `HTTP ${xhr.status}: ${xhr.statusText || "upload failed"}`;
          if (isJson) {
            try {
              const body = JSON.parse(xhr.responseText);
              if (body?.error) message = body.error;
              else if (body?.detail) message = body.detail;
            } catch {
              /* keep the status line */
            }
          }

          if (xhr.status === 401 || xhr.status === 403) {
            emitAuthError({
              status: xhr.status as 401 | 403,
              url,
              message:
                xhr.status === 401
                  ? "Your session has expired. Please sign in again."
                  : "You don't have permission to access this resource.",
            });
          }

          fail(new Error(message));
        });

        xhr.open("POST", resolvedUrl, true);
        if (token) xhr.setRequestHeader("Authorization", `Bearer ${token}`);
        // Content-Type is deliberately left unset: the browser must add
        // the multipart boundary itself.
        xhr.send(formData);
      })
      .catch((error) => reject(error));
  });
}
