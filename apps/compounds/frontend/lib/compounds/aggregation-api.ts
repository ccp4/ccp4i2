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
import { Target, SavedAggregationView } from '@/types/compounds/models';

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
 * Format KPI unit for display.
 * Converts ASCII unit codes to proper display format (e.g., 'uM' -> 'μM').
 */
export function formatKpiUnit(unit: string | null | undefined): string {
  if (!unit) return '';

  // Convert common unit formats to display format with proper symbols
  const unitMap: Record<string, string> = {
    // Molar concentrations
    'nM': 'nM',
    'uM': 'μM',  // ASCII 'u' to Greek mu
    'mM': 'mM',
    'pM': 'pM',
    'M': 'M',
    // Time units
    'min': 'min',
    's': 's',
    'h': 'h',
    // Percentage
    '%': '%',
    // Rate units
    'uL/min/mg': 'μL/min/mg',
    'mL/min/kg': 'mL/min/kg',
    // Permeability
    '1e-6 cm/s': '10⁻⁶ cm/s',
    'cm/s': 'cm/s',
  };

  return unitMap[unit] || unit;
}

// =============================================================================
// Concentration Display Conversion
// =============================================================================

import { ConcentrationDisplayMode } from '@/types/compounds/aggregation';

/** Molar concentration units that can be converted */
const CONCENTRATION_UNITS = ['pM', 'nM', 'uM', 'mM', 'M'] as const;

/** Conversion factors to molar (M) */
const UNIT_TO_MOLAR: Record<string, number> = {
  'pM': 1e-12,
  'nM': 1e-9,
  'uM': 1e-6,
  'mM': 1e-3,
  'M': 1,
};

/**
 * Check if a unit is a molar concentration unit.
 */
export function isConcentrationUnit(unit: string | null | undefined): boolean {
  if (!unit) return false;
  return CONCENTRATION_UNITS.includes(unit as typeof CONCENTRATION_UNITS[number]);
}

/**
 * Convert a concentration value from one unit to another.
 * Returns the original value if units are unknown or the same.
 */
export function convertConcentration(
  value: number,
  fromUnit: string,
  toUnit: string
): number {
  if (fromUnit === toUnit) return value;

  const fromFactor = UNIT_TO_MOLAR[fromUnit];
  const toFactor = UNIT_TO_MOLAR[toUnit];

  if (!fromFactor || !toFactor) return value;

  // Convert to molar, then to target unit
  const molar = value * fromFactor;
  return molar / toFactor;
}

/**
 * Convert a concentration value to pConc (-log10 of molar concentration).
 * Returns null if the unit is unknown or value is non-positive.
 */
export function toPConc(value: number, unit: string): number | null {
  const factor = UNIT_TO_MOLAR[unit];
  if (!factor || value <= 0) return null;

  const molar = value * factor;
  return -Math.log10(molar);
}

/**
 * Format a concentration value according to the display mode.
 * Returns the formatted value string and the display unit.
 */
export function formatConcentrationValue(
  value: number | null | undefined,
  unit: string | null | undefined,
  displayMode: ConcentrationDisplayMode
): { displayValue: string; displayUnit: string } {
  // Handle null/undefined values
  if (value === null || value === undefined) {
    return { displayValue: '-', displayUnit: '' };
  }

  // If not a concentration unit or mode is natural, use standard formatting
  if (!unit || !isConcentrationUnit(unit) || displayMode === 'natural') {
    return {
      displayValue: formatKpiValue(value),
      displayUnit: unit ? formatKpiUnit(unit) : '',
    };
  }

  // Handle pConc mode
  if (displayMode === 'pConc') {
    const pValue = toPConc(value, unit);
    if (pValue === null) {
      return { displayValue: '-', displayUnit: '' };
    }
    return {
      displayValue: pValue.toFixed(2),
      displayUnit: '', // pConc doesn't have a unit suffix
    };
  }

  // Handle standardized unit modes (nM, uM, mM)
  const convertedValue = convertConcentration(value, unit, displayMode);
  return {
    displayValue: formatKpiValue(convertedValue),
    displayUnit: formatKpiUnit(displayMode),
  };
}

/**
 * Get the display unit for a column header based on the concentration display mode.
 * Returns the formatted unit string for the header.
 */
export function getConcentrationHeaderUnit(
  unit: string | null | undefined,
  displayMode: ConcentrationDisplayMode
): string {
  if (!unit || !isConcentrationUnit(unit) || displayMode === 'natural') {
    return unit ? formatKpiUnit(unit) : '';
  }

  if (displayMode === 'pConc') {
    return ''; // pConc uses "p" prefix on value name, not a unit suffix
  }

  return formatKpiUnit(displayMode);
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

  // Add columns for each protocol and aggregation (include unit for value-based aggs)
  for (const protocol of protocols) {
    for (const agg of aggregations) {
      const unitSuffix = agg !== 'count' && agg !== 'list' && protocol.kpi_unit
        ? ` [${formatKpiUnit(protocol.kpi_unit)}]`
        : '';
      headers.push(`${protocol.name} (${agg})${unitSuffix}`);
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
    'KPI Unit',
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
      `"${row.kpi_unit ? formatKpiUnit(row.kpi_unit) : ''}"`,
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

// =============================================================================
// Saved Aggregation View API
// =============================================================================

/**
 * Save an aggregation view configuration to a target.
 * Requires admin operating level.
 */
export async function saveAggregationView(
  targetId: string,
  config: SavedAggregationView
): Promise<{ success: boolean; saved_view: SavedAggregationView }> {
  const res = await coreFetch(`${API_BASE}/targets/${targetId}/saved_view/`, {
    method: 'POST',
    body: JSON.stringify(config),
  });

  if (!res.ok) {
    const error = await res.json().catch(() => ({ error: 'Request failed' }));
    throw new Error(error.error || 'Failed to save aggregation view');
  }

  return res.json();
}

/**
 * Delete a saved aggregation view from a target.
 * Requires admin operating level.
 */
export async function deleteAggregationView(targetId: string): Promise<void> {
  const res = await coreFetch(`${API_BASE}/targets/${targetId}/saved_view/`, {
    method: 'DELETE',
  });

  if (!res.ok) {
    const error = await res.json().catch(() => ({ error: 'Request failed' }));
    throw new Error(error.error || 'Failed to delete aggregation view');
  }
}
