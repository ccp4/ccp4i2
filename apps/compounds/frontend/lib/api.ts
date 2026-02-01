/**
 * API utilities for the compounds frontend.
 * Uses SWR for data fetching with caching.
 */

import useSWR, { SWRConfiguration, SWRResponse } from 'swr';

// Use /api/proxy/compounds to match CCP4i2 client pattern
// In standalone mode, next.config.js rewrites this to Django's /compounds/
// When integrated with CCP4i2 client, it goes through the same proxy as other APIs
const API_BASE = '/api/proxy/compounds';

/**
 * Standard fetcher for SWR
 */
async function fetcher<T>(url: string): Promise<T> {
  const res = await fetch(url);
  if (!res.ok) {
    const error = new Error('API request failed');
    throw error;
  }
  return res.json();
}

/**
 * POST request helper
 */
export async function apiPost<T>(endpoint: string, body: any): Promise<T> {
  const res = await fetch(`${API_BASE}/${endpoint}`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    body: JSON.stringify(body),
  });
  if (!res.ok) {
    const error = new Error('API request failed');
    throw error;
  }
  return res.json();
}

/**
 * DELETE request helper
 */
export async function apiDelete(endpoint: string): Promise<void> {
  const res = await fetch(`${API_BASE}/${endpoint}`, {
    method: 'DELETE',
  });
  if (!res.ok) {
    const error = new Error('API request failed');
    throw error;
  }
}

/**
 * PATCH request helper
 */
export async function apiPatch<T>(endpoint: string, body: any): Promise<T> {
  const res = await fetch(`${API_BASE}/${endpoint}`, {
    method: 'PATCH',
    headers: {
      'Content-Type': 'application/json',
    },
    body: JSON.stringify(body),
  });
  if (!res.ok) {
    const error = new Error('API request failed');
    throw error;
  }
  return res.json();
}

/**
 * Hook for fetching data from the compounds API
 */
export function useCompoundsApi() {
  return {
    /**
     * GET request with SWR caching
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
  };
}
