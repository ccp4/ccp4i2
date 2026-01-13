/**
 * API utilities for the users frontend.
 * Uses the same authentication pattern as compounds/api.ts
 */

// =============================================================================
// Authentication Integration
// =============================================================================

// For Docker integration: Try to import getAccessToken from ccp4i2 client's auth-token
// Falls back to no-op for standalone development
let getAccessToken: () => Promise<string | null>;

try {
  // This path works when overlaid into ccp4i2 client (renderer/lib/users/)
  const authModule = require('../../utils/auth-token');
  getAccessToken = authModule.getAccessToken;
} catch {
  // Fallback for standalone development - no auth needed
  getAccessToken = async () => null;
}

// =============================================================================
// Configuration
// =============================================================================

// API base path - uses /api/proxy/users to match CCP4i2 client pattern
const API_BASE = '/api/proxy/users';

/**
 * Core fetch wrapper with authentication support.
 */
async function coreFetch(
  url: string,
  options: RequestInit = {}
): Promise<Response> {
  const headers: Record<string, string> = {};

  // Inject authentication token if available
  const token = await getAccessToken();
  if (token) {
    headers['Authorization'] = `Bearer ${token}`;
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
 * GET request helper
 */
export async function apiGet<T>(endpoint: string): Promise<T> {
  const res = await coreFetch(`${API_BASE}/${endpoint}`);
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
export async function apiPost<T>(endpoint: string, body?: any): Promise<T> {
  const res = await coreFetch(`${API_BASE}/${endpoint}`, {
    method: 'POST',
    body: body ? JSON.stringify(body) : undefined,
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
 * PATCH request helper
 */
export async function apiPatch<T>(endpoint: string, body: any): Promise<T> {
  const res = await coreFetch(`${API_BASE}/${endpoint}`, {
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
 * DELETE request helper
 */
export async function apiDelete(endpoint: string): Promise<void> {
  const res = await coreFetch(`${API_BASE}/${endpoint}`, {
    method: 'DELETE',
  });
  if (!res.ok) {
    const errorText = await res.text().catch(() => 'Unknown error');
    const error = new Error(`API request failed: ${res.status} ${errorText}`);
    (error as any).status = res.status;
    throw error;
  }
}
