/**
 * API utilities for the compounds frontend.
 * Uses SWR for data fetching with caching.
 *
 * This module mirrors the pattern in client/renderer/api.ts but is specialized
 * for compounds endpoints.
 *
 * Authentication Integration:
 * - Standalone (development): No auth required, works without tokens
 * - Integrated (Docker): Dockerfile runs sed to uncomment the auth import
 */

import useSWR, { SWRConfiguration, SWRResponse } from 'swr';

// =============================================================================
// Authentication Integration
// =============================================================================

// authFetch was a compounds-side fetch wrapper that injected the
// auth token + X-User-Email header. The implementation moved to
// ./api-fetch (a thin binding around @ccp4/ccp4i2-auth's
// createApiFetch factory). Imported + re-exported here so existing
// call sites (`import { authFetch } from "../lib/compounds/api"`)
// keep working.
//
// IMPORTANT BEHAVIOUR CHANGE: the new authFetch throws on non-OK
// responses (the canonical apiFetch semantic). Callers that wrap
// authFetch in `if (!res.ok) throw` blocks are now defensively dead
// — the throw fires inside authFetch, which suits the auth-error-
// event flow but makes those manual checks redundant. Cleanup of
// such dead checks is a separate concern from this migration.
import { authFetch } from "./api-fetch";
export { authFetch };

// =============================================================================
// Configuration
// =============================================================================

// API base path - uses /api/proxy/compounds to match CCP4i2 client pattern
// In standalone mode, next.config.js rewrites this to Django's /api/compounds/
// When integrated with CCP4i2 client, it goes through the same proxy as other APIs
const API_BASE = '/api/proxy/compounds';

// =============================================================================
// API Request Helpers
// =============================================================================

/**
 * Standard fetcher for SWR - returns parsed JSON
 */
async function fetcher<T>(url: string): Promise<T> {
  const res = await authFetch(url);
  if (!res.ok) {
    const errorText = await res.text().catch(() => 'Unknown error');
    const error = new Error(`API request failed: ${res.status} ${errorText}`);
    (error as any).status = res.status;
    throw error;
  }
  return res.json();
}

/**
 * GET request helper (for one-off requests outside of SWR)
 */
export async function apiGet<T>(endpoint: string): Promise<T> {
  const res = await authFetch(`${API_BASE}/${endpoint}`);
  if (!res.ok) {
    const errorText = await res.text().catch(() => 'Unknown error');
    const error = new Error(`API request failed: ${res.status} ${errorText}`);
    (error as any).status = res.status;
    throw error;
  }
  return res.json();
}

/**
 * POST request helper
 */
export async function apiPost<T>(endpoint: string, body: any): Promise<T> {
  const res = await authFetch(`${API_BASE}/${endpoint}`, {
    method: 'POST',
    body: body instanceof FormData ? body : JSON.stringify(body),
  });
  if (!res.ok) {
    const errorText = await res.text().catch(() => 'Unknown error');
    const error = new Error(`API request failed: ${res.status} ${errorText}`);
    (error as any).status = res.status;
    throw error;
  }
  return res.json();
}

/**
 * DELETE request helper
 */
export async function apiDelete(endpoint: string): Promise<void> {
  const res = await authFetch(`${API_BASE}/${endpoint}`, {
    method: 'DELETE',
  });
  if (!res.ok) {
    const errorText = await res.text().catch(() => 'Unknown error');
    const error = new Error(`API request failed: ${res.status} ${errorText}`);
    (error as any).status = res.status;
    throw error;
  }
}

/**
 * PATCH request helper
 */
export async function apiPatch<T>(endpoint: string, body: any): Promise<T> {
  const res = await authFetch(`${API_BASE}/${endpoint}`, {
    method: 'PATCH',
    body: JSON.stringify(body),
  });
  if (!res.ok) {
    const errorText = await res.text().catch(() => 'Unknown error');
    const error = new Error(`API request failed: ${res.status} ${errorText}`);
    (error as any).status = res.status;
    throw error;
  }
  return res.json();
}

/**
 * File upload helper (multipart/form-data)
 */
export async function apiUpload<T>(endpoint: string, formData: FormData): Promise<T> {
  const res = await authFetch(`${API_BASE}/${endpoint}`, {
    method: 'POST',
    body: formData,
    // Don't set Content-Type - let browser set it with boundary
  });
  if (!res.ok) {
    const errorText = await res.text().catch(() => 'Unknown error');
    const error = new Error(`API request failed: ${res.status} ${errorText}`);
    (error as any).status = res.status;
    throw error;
  }
  return res.json();
}

/**
 * Fetch an authenticated file and open it in a new tab (or trigger a save
 * dialog if `filename` is provided). The token travels in the Authorization
 * header, never in the URL — this avoids 431 Request Header Too Large when
 * the JWT's groups claim is bulky, and keeps the token out of browser history
 * and proxy logs.
 *
 * Must be invoked from a user-gesture handler (onClick) so popup blockers
 * allow `window.open`.
 */
export async function openAuthenticatedDownload(url: string, filename?: string | null): Promise<void> {
  const res = await authFetch(url);
  if (!res.ok) {
    const errorText = await res.text().catch(() => 'Unknown error');
    const error = new Error(`Download failed: ${res.status} ${errorText}`);
    (error as any).status = res.status;
    throw error;
  }
  const blob = await res.blob();
  const objectUrl = URL.createObjectURL(blob);

  if (filename) {
    const a = document.createElement('a');
    a.href = objectUrl;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
  } else {
    window.open(objectUrl, '_blank');
  }

  setTimeout(() => URL.revokeObjectURL(objectUrl), 60_000);
}

// =============================================================================
// Main API Hook
// =============================================================================

/**
 * Hook for fetching data from the compounds API.
 * Uses SWR for caching and revalidation.
 */
export function useCompoundsApi() {
  return {
    /**
     * GET request with SWR caching
     * Pass null endpoint to skip fetching (conditional fetch pattern)
     */
    get<T>(endpoint: string | null, config?: SWRConfiguration): SWRResponse<T> {
      const url = endpoint ? `${API_BASE}/${endpoint}` : null;
      return useSWR<T>(url, fetcher, config);
    },

    /**
     * POST request (mutation)
     */
    post: apiPost,

    /**
     * DELETE request (mutation)
     */
    delete: apiDelete,

    /**
     * PATCH request (mutation)
     */
    patch: apiPatch,

    /**
     * File upload (mutation)
     */
    upload: apiUpload,
  };
}

// =============================================================================
// Integration Notes
// =============================================================================
//
// When integrating with the main ccp4i2 client:
//
// 1. Copy this file to: client/renderer/lib/compounds-api.ts
//
// 2. Replace getAccessToken() with:
//    import { getAccessToken } from '../utils/auth-token';
//
// 3. The API_BASE path ('/api/proxy/compounds') will work as-is because
//    the main ccp4i2 proxy at /api/proxy/[...proxy]/route.ts can be
//    extended to handle 'compounds' paths, OR a separate proxy route
//    can be added at /api/proxy/compounds/[...path]/route.ts
//
// 4. All authentication tokens will flow through automatically because
//    authFetch injects the Authorization header.
