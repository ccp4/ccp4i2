/**
 * Aggregation API utilities.
 * Provides functions and hooks for the aggregation endpoint.
 */

import useSWR, { SWRConfiguration, SWRResponse } from 'swr';
import {
  AggregationRequest,
  AggregationResponse,
  ProtocolInfo,
} from '@/types/aggregation';
import { Target } from '@/types/models';

const API_BASE = '/api/proxy/compounds';

/**
 * Fetch aggregation data.
 * Uses POST because the request body can be complex.
 */
export async function fetchAggregation(
  request: AggregationRequest
): Promise<AggregationResponse> {
  const res = await fetch(`${API_BASE}/aggregations/aggregate/`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
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
  const res = await fetch(url);

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
  const res = await fetch(url);

  if (!res.ok) {
    throw new Error('Failed to fetch targets');
  }

  return res.json();
}

/**
 * Format a KPI value for display.
 * Uses exponential notation for very small or very large numbers.
 */
export function formatKpiValue(value: number | string | null | undefined): string {
  if (value === null || value === undefined || value === '') {
    return '-';
  }

  // Convert string values to numbers (handles legacy data stored as strings)
  const numValue = typeof value === 'string' ? parseFloat(value) : value;
  if (isNaN(numValue)) {
    return String(value);
  }

  // Use exponential for very small or very large numbers
  if (numValue !== 0 && (Math.abs(numValue) < 0.01 || Math.abs(numValue) > 100000)) {
    return numValue.toExponential(2);
  }

  // Regular formatting with 2 decimal places
  return numValue.toFixed(2);
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
