/**
 * Aggregation API utilities.
 * Provides functions and hooks for the aggregation endpoint.
 *
 * Authentication Integration:
 * - Standalone (development): No auth required, works without tokens
 * - Integrated (Docker): Uses auth-token module from ccp4i2 client
 */

import useSWR, { SWRConfiguration, SWRResponse } from 'swr';
import {
  AggregationRequest,
  AggregationResponse,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import { Target } from '@/types/compounds/models';

// =============================================================================
// Authentication Integration
// =============================================================================

// For Docker integration: Try to import auth helpers from ccp4i2 client's auth-token
// Falls back to no-op for standalone development
let getAccessToken: () => Promise<string | null>;
let getUserEmail: () => string | null;

try {
  // This path works when overlaid into ccp4i2 client (renderer/lib/compounds/)
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

const API_BASE = '/api/proxy/compounds';

/**
 * Core fetch wrapper with authentication support.
 * Mirrors the pattern in api.ts
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

/**
 * Fetch aggregation data.
 * Uses POST because the request body can be complex.
 */
export async function fetchAggregation(
  request: AggregationRequest
): Promise<AggregationResponse> {
  const res = await coreFetch(`${API_BASE}/aggregations/aggregate/`, {
    method: 'POST',
    body: JSON.stringify(request),
  });

  if (!res.ok) {
    const error = await res.json().catch(() => ({ error: 'Request failed' }));
    throw new Error(error.error || 'Aggregation request failed');
  }

  return res.json();
}

/**
 * Stable key generator for SWR caching.
 * Converts request object to a stable string key.
 */
function getAggregationKey(request: AggregationRequest | null): string | null {
  if (!request) return null;
  return `aggregation:${JSON.stringify(request)}`;
}

/**
 * SWR hook for aggregation data.
 * Caches results and handles loading/error states.
 */
export function useAggregation(
  request: AggregationRequest | null,
  config?: SWRConfiguration
): SWRResponse<AggregationResponse | null> {
  return useSWR<AggregationResponse | null>(
    getAggregationKey(request),
    () => (request ? fetchAggregation(request) : null),
    {
      revalidateOnFocus: false,
      ...config,
    }
  );
}

/**
 * Fetch available protocols for predicate builder.
 */
export async function fetchProtocols(params?: {
  target?: string;
  search?: string;
}): Promise<ProtocolInfo[]> {
  const searchParams = new URLSearchParams();
  if (params?.target) searchParams.set('target', params.target);
  if (params?.search) searchParams.set('search', params.search);

  const url = `${API_BASE}/aggregations/protocols/?${searchParams}`;
  const res = await coreFetch(url);

  if (!res.ok) {
    throw new Error('Failed to fetch protocols');
  }

  return res.json();
}

/**
 * Fetch available targets for predicate builder.
 */
export async function fetchTargets(params?: {
  search?: string;
}): Promise<Target[]> {
  const searchParams = new URLSearchParams();
  if (params?.search) searchParams.set('search', params.search);

  const url = `${API_BASE}/aggregations/targets/?${searchParams}`;
  const res = await coreFetch(url);

  if (!res.ok) {
    throw new Error('Failed to fetch targets');
  }

  return res.json();
}

/**
 * Format a KPI value for display.
 * Uses exponential notation for very small or very large numbers.
 */
export function formatKpiValue(value: number | null | undefined): string {
  if (value === null || value === undefined) {
    return '-';
  }

  // Use exponential for very small or very large numbers
  if (value !== 0 && (Math.abs(value) < 0.01 || Math.abs(value) > 100000)) {
    return value.toExponential(2);
  }

  // Regular formatting with 2 decimal places
  return value.toFixed(2);
}

/**
 * Generate CSV content from compact aggregation data.
 */
export function generateCompactCsv(
  data: AggregationResponse & { protocols: ProtocolInfo[] },
  aggregations: string[]
): string {
  const { protocols, data: rows } = data;

  // Build header row
  const headers = ['Compound ID', 'SMILES', 'Target'];

  // Add columns for each protocol and aggregation
  for (const protocol of protocols) {
    for (const agg of aggregations) {
      headers.push(`${protocol.name} (${agg})`);
    }
  }

  const csvRows = [headers.join(',')];

  // Build data rows
  for (const row of rows as any[]) {
    const values = [
      `"${row.formatted_id}"`,
      `"${row.smiles || ''}"`,
      `"${row.target_name || ''}"`,
    ];

    for (const protocol of protocols) {
      const protocolData = row.protocols[protocol.id] || {};
      for (const agg of aggregations) {
        const value = protocolData[agg];
        if (agg === 'list') {
          values.push(`"${value || ''}"`);
        } else {
          values.push(value !== null && value !== undefined ? String(value) : '');
        }
      }
    }

    csvRows.push(values.join(','));
  }

  return csvRows.join('\n');
}

/**
 * Generate CSV content from long aggregation data.
 */
export function generateLongCsv(data: AggregationResponse): string {
  const { data: rows } = data;

  const headers = [
    'Compound ID',
    'Compound Name',
    'SMILES',
    'Target',
    'Protocol',
    'Assay Date',
    'KPI Value',
    'Status',
  ];

  const csvRows = [headers.join(',')];

  for (const row of rows as any[]) {
    const values = [
      `"${row.formatted_id || ''}"`,
      `"${row.compound_name || ''}"`,
      `"${row.smiles || ''}"`,
      `"${row.target_name || ''}"`,
      `"${row.protocol_name || ''}"`,
      `"${row.assay_date || ''}"`,
      row.kpi_value !== null && row.kpi_value !== undefined ? String(row.kpi_value) : '',
      `"${row.status || ''}"`,
    ];

    csvRows.push(values.join(','));
  }

  return csvRows.join('\n');
}

/**
 * Download CSV file.
 */
export function downloadCsv(content: string, filename: string): void {
  const blob = new Blob([content], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}
