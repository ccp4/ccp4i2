/**
 * Aggregation API type definitions.
 * Types for data aggregation queries and responses.
 */

import { AnalysisStatus } from './models';

/** Aggregation function types */
export type AggregationType = 'geomean' | 'count' | 'stdev' | 'list';

/** Output format options */
export type OutputFormat = 'compact' | 'medium' | 'long';

/** Filter predicates for querying data series */
export interface Predicates {
  /** Target UUIDs to filter by */
  targets?: string[];
  /** Compound UUIDs to filter by */
  compounds?: string[];
  /** Text search for compound formatted_id (e.g., NCL-00026) */
  compound_search?: string;
  /** Protocol UUIDs to filter by */
  protocols?: string[];
  /** Analysis status filter (default: 'valid') */
  status?: AnalysisStatus | '';
}

/** Request body for aggregation API */
export interface AggregationRequest {
  predicates?: Predicates;
  output_format: OutputFormat;
  aggregations: AggregationType[];
}

/** Protocol info in response */
export interface ProtocolInfo {
  id: string;
  name: string;
}

/** Aggregated values for a single protocol */
export interface ProtocolAggregation {
  geomean?: number | null;
  count?: number;
  stdev?: number | null;
  list?: string;
}

/** Summary metadata for response */
export interface AggregationMeta {
  compound_count: number;
  protocol_count: number;
  total_measurements: number;
}

/** A single compound row in compact format */
export interface CompactRow {
  compound_id: string;
  formatted_id: string;
  smiles: string | null;
  target_name: string | null;
  /** Protocol aggregations keyed by protocol UUID */
  protocols: Record<string, ProtocolAggregation>;
}

/** Compact format response (one row per compound) */
export interface CompactAggregationResponse {
  meta: AggregationMeta;
  protocols: ProtocolInfo[];
  data: CompactRow[];
}

/** A single row in medium format (one per compound-protocol pair) */
export interface MediumRow {
  compound_id: string;
  formatted_id: string;
  smiles: string | null;
  target_name: string | null;
  protocol_id: string;
  protocol_name: string;
  /** Aggregated values inline */
  geomean?: number | null;
  count?: number;
  stdev?: number | null;
  list?: string;
}

/** Medium format response (one row per compound-protocol pair) */
export interface MediumAggregationResponse {
  meta: AggregationMeta;
  data: MediumRow[];
}

/** A single measurement row in long format */
export interface LongRow {
  data_series_id: string;
  compound_id: string | null;
  formatted_id: string | null;
  compound_name: string | null;
  smiles: string | null;
  target_name: string | null;
  protocol_id: string;
  protocol_name: string;
  assay_id: string;
  assay_date: string | null;
  kpi_value: number | null;
  status: AnalysisStatus | null;
}

/** Long format response (one row per measurement) */
export interface LongAggregationResponse {
  meta: AggregationMeta;
  data: LongRow[];
}

/** Union type for aggregation responses */
export type AggregationResponse = CompactAggregationResponse | MediumAggregationResponse | LongAggregationResponse;

/** Type guard for compact response */
export function isCompactResponse(
  response: AggregationResponse
): response is CompactAggregationResponse {
  return 'protocols' in response;
}

/** Type guard for medium response */
export function isMediumResponse(
  response: AggregationResponse
): response is MediumAggregationResponse {
  return !('protocols' in response) && response.data.length > 0 && 'protocol_name' in response.data[0];
}

/** Type guard for long response */
export function isLongResponse(
  response: AggregationResponse
): response is LongAggregationResponse {
  return !('protocols' in response) && (response.data.length === 0 || 'data_series_id' in response.data[0]);
}
