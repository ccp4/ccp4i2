/**
 * Aggregation API type definitions.
 * Types for data aggregation queries and responses.
 */

import { AnalysisStatus, AnalysisStatusFilter } from './models';

/** Aggregation function types */
export type AggregationType = 'geomean' | 'count' | 'stdev' | 'list';

/** Molecular property names that can be included in aggregation tables */
export type MolecularPropertyName =
  | 'molecular_weight'
  | 'heavy_atom_count'
  | 'hbd'
  | 'hba'
  | 'clogp'
  | 'tpsa'
  | 'rotatable_bonds'
  | 'fraction_sp3';

/** Molecular property values for a compound */
export interface MolecularPropertyValues {
  molecular_weight?: number | null;
  heavy_atom_count?: number | null;
  hbd?: number | null;
  hba?: number | null;
  clogp?: number | null;
  tpsa?: number | null;
  rotatable_bonds?: number | null;
  fraction_sp3?: number | null;
}

/** RAG threshold direction */
export type ThresholdDirection = 'above' | 'below';

/** RAG status */
export type RagStatus = 'green' | 'amber' | 'red';

/** RAG threshold configuration for a molecular property */
export interface MolecularPropertyThreshold {
  id: string;
  property_name: MolecularPropertyName;
  property_display: string;
  direction: ThresholdDirection;
  direction_display: string;
  amber_threshold: number;
  red_threshold: number;
  enabled: boolean;
}

/** Molecular property metadata */
export interface MolecularPropertyMeta {
  name: MolecularPropertyName;
  display_name: string;
  description: string;
}

/** Concentration display mode for KPI values */
export type ConcentrationDisplayMode = 'natural' | 'nM' | 'uM' | 'mM' | 'pConc';

/** Output format options */
export type OutputFormat = 'compact' | 'medium' | 'long' | 'pivot' | 'cards';

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
  /** Analysis status filter (default: 'valid'). Use 'no_analysis' for DataSeries where analysis failed. */
  status?: AnalysisStatusFilter | '';
}

/** Request body for aggregation API */
export interface AggregationRequest {
  predicates?: Predicates;
  output_format: OutputFormat;
  aggregations: AggregationType[];
  /** When true, split results by batch (creates separate rows for each batch) */
  group_by_batch?: boolean;
  /** When true, include compounds tested but with no valid KPI values (shown with count=0) */
  include_tested_no_data?: boolean;
  /** Molecular properties to include as columns in the aggregation table */
  include_properties?: MolecularPropertyName[];
}

/** Protocol info in response */
export interface ProtocolInfo {
  id: string;
  name: string;
  /** KPI unit for this protocol (e.g., 'nM', 'uM', 'mM') */
  kpi_unit?: string | null;
}

/** Aggregated values for a single protocol */
export interface ProtocolAggregation {
  geomean?: number | null;
  count?: number;
  stdev?: number | null;
  /** Standard deviation in log10 space, for pConc display */
  stdev_log?: number | null;
  list?: string;
  /** Total number of DataSeries tested (regardless of analysis status) */
  tested?: number;
  /** Number of DataSeries where analysis failed (no AnalysisResult) */
  no_analysis?: number;
  /** Number of DataSeries with status='invalid' */
  invalid?: number;
  /** Number of DataSeries with status='unassigned' */
  unassigned?: number;
}

/** Summary metadata for response */
export interface AggregationMeta {
  compound_count: number;
  /** Number of rows (may differ from compound_count when group_by_batch=true) */
  row_count?: number;
  protocol_count: number;
  total_measurements: number;
  /** Whether results are grouped by batch */
  group_by_batch?: boolean;
  /** Whether results include compounds tested but with no valid KPI values */
  include_tested_no_data?: boolean;
  /** Molecular properties included in this response */
  include_properties?: MolecularPropertyName[];
}

/** A single compound row in compact format */
export interface CompactRow {
  compound_id: string;
  formatted_id: string;
  smiles: string | null;
  target_name: string | null;
  /** Batch UUID (only present when group_by_batch=true) */
  batch_id?: string | null;
  /** Batch number (only present when group_by_batch=true) */
  batch_number?: number | null;
  /** Protocol aggregations keyed by protocol UUID */
  protocols: Record<string, ProtocolAggregation>;
  /** Molecular properties (only present when include_properties is specified) */
  properties?: MolecularPropertyValues;
}

/** Compact format response (one row per compound) */
export interface CompactAggregationResponse {
  meta: AggregationMeta;
  protocols: ProtocolInfo[];
  data: CompactRow[];
  /** RAG thresholds for included properties (only present when include_properties is specified) */
  property_thresholds?: MolecularPropertyThreshold[];
}

/** A single row in medium format (one per compound-protocol pair) */
export interface MediumRow {
  compound_id: string;
  formatted_id: string;
  smiles: string | null;
  target_name: string | null;
  /** Batch UUID (only present when group_by_batch=true) */
  batch_id?: string | null;
  /** Batch number (only present when group_by_batch=true) */
  batch_number?: number | null;
  protocol_id: string;
  protocol_name: string;
  /** KPI unit for this row (e.g., 'nM', 'uM', 'mM') */
  kpi_unit?: string | null;
  /** Aggregated values inline */
  geomean?: number | null;
  count?: number;
  stdev?: number | null;
  /** Standard deviation in log10 space, for pConc display */
  stdev_log?: number | null;
  list?: string;
  /** Total number of DataSeries tested (regardless of analysis status) */
  tested?: number;
  /** Number of DataSeries where analysis failed (no AnalysisResult) */
  no_analysis?: number;
  /** Number of DataSeries with status='invalid' */
  invalid?: number;
  /** Number of DataSeries with status='unassigned' */
  unassigned?: number;
  /** Molecular properties (only present when include_properties is specified) */
  properties?: MolecularPropertyValues;
}

/** Medium format response (one row per compound-protocol pair) */
export interface MediumAggregationResponse {
  meta: AggregationMeta;
  data: MediumRow[];
  /** RAG thresholds for included properties (only present when include_properties is specified) */
  property_thresholds?: MolecularPropertyThreshold[];
}

/** A single measurement row in long format */
export interface LongRow {
  data_series_id: string;
  compound_id: string | null;
  formatted_id: string | null;
  compound_name: string | null;
  smiles: string | null;
  target_name: string | null;
  /** Batch UUID (always present in long format) */
  batch_id?: string | null;
  /** Batch number (always present in long format) */
  batch_number?: number | null;
  protocol_id: string;
  protocol_name: string;
  assay_id: string;
  assay_date: string | null;
  kpi_value: number | null;
  /** KPI unit for this measurement (e.g., 'nM', 'uM', 'mM') */
  kpi_unit?: string | null;
  status: AnalysisStatus | null;
  /** Molecular properties (only present when include_properties is specified) */
  properties?: MolecularPropertyValues;
}

/** Long format response (one row per measurement) */
export interface LongAggregationResponse {
  meta: AggregationMeta;
  data: LongRow[];
  /** RAG thresholds for included properties (only present when include_properties is specified) */
  property_thresholds?: MolecularPropertyThreshold[];
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

/**
 * Determine RAG status for a property value based on thresholds.
 * @param value The property value to evaluate
 * @param threshold The threshold configuration for this property
 * @returns RAG status ('green', 'amber', or 'red')
 */
export function getRagStatus(
  value: number | null | undefined,
  threshold: MolecularPropertyThreshold | undefined
): RagStatus {
  if (value == null || !threshold || !threshold.enabled) {
    return 'green';
  }

  if (threshold.direction === 'above') {
    if (value >= threshold.red_threshold) return 'red';
    if (value >= threshold.amber_threshold) return 'amber';
    return 'green';
  } else {
    // direction === 'below'
    if (value <= threshold.red_threshold) return 'red';
    if (value <= threshold.amber_threshold) return 'amber';
    return 'green';
  }
}

/** Response from /aggregations/molecular_properties/ endpoint */
export interface MolecularPropertiesResponse {
  properties: MolecularPropertyMeta[];
  thresholds: MolecularPropertyThreshold[];
}
