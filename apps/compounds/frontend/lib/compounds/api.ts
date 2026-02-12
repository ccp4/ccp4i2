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

// For Docker integration: Try to import auth helpers from ccp4i2 client's auth-token
// Falls back to no-op for standalone development
let getAccessToken: () => Promise<string | null>;
let getUserEmail: () => string | null;

try {
  // This path works when overlaid into ccp4i2 client (renderer/lib/compounds/)
  // The dynamic import is resolved at build time by webpack
  const authModule = require('../../utils/auth-token');
  getAccessToken = authModule.getAccessToken;
  getUserEmail = authModule.getUserEmail || (() => null);
} catch {
  // Fallback for standalone development - no auth needed
  getAccessToken = async () => null;
  getUserEmail = () => null;
}

// =============================================================================
// Configuration
// =============================================================================

// API base path - uses /api/proxy/compounds to match CCP4i2 client pattern
// In standalone mode, next.config.js rewrites this to Django's /api/compounds/
// When integrated with CCP4i2 client, it goes through the same proxy as other APIs
const API_BASE = '/api/proxy/compounds';

/**
 * Core fetch wrapper with authentication support.
 * Mirrors the pattern in client/renderer/api-fetch.ts
 *
 * Use this for fetching any protected resource (images, files, API endpoints)
 * that requires authentication headers.
 */
export async function authFetch(
  url: string,
  options: RequestInit = {}
): Promise<Response> {
  const headers: Record<string, string> = {};

  // Inject authentication token if available
  const token = await getAccessToken();
  if (token) {
    headers['Authorization'] = `Bearer ${token}`;
  }

  // Include user email as fallback for backends where access tokens don't have email claims
  const email = getUserEmail();
  if (email) {
    headers['X-User-Email'] = email;
  }

  // Merge with provided headers
  if (options.headers) {
    if (options.headers instanceof Headers) {
      options.headers.forEach((value, key) => {
        headers[key] = value;
      });
    } else if (Array.isArray(options.headers)) {
      options.headers.forEach(([key, value]) => {
        headers[key] = value;
      });
    } else {
      Object.assign(headers, options.headers);
    }
  }

  // Set Content-Type for JSON bodies (but not FormData)
  if (options.body && !(options.body instanceof FormData)) {
    if (!Object.keys(headers).some(k => k.toLowerCase() === 'content-type')) {
      headers['Content-Type'] = 'application/json';
    }
  }

  return fetch(url, {
    ...options,
    headers,
  });
}

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
 * Get an authenticated download URL.
 * Appends the access token as a query parameter for anchor-based downloads
 * where headers cannot be set.
 */
export async function getAuthenticatedDownloadUrl(url: string): Promise<string> {
  const token = await getAccessToken();
  if (!token) {
    return url; // No auth needed in standalone mode
  }
  const separator = url.includes('?') ? '&' : '?';
  return `${url}${separator}access_token=${encodeURIComponent(token)}`;
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
