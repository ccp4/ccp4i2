/**
 * Compounds Registry and Assays type definitions.
 *
 * These map to the Django models in apps/compounds/registry and apps/compounds/assays.
 */

// =============================================================================
// Registry Models
// =============================================================================

export interface Supplier {
  id: string;  // UUID
  name: string;
  initials: string | null;
}

export interface Target {
  id: string;  // UUID
  name: string;
  parent: string | null;  // UUID of parent target
  created_at: string;  // ISO datetime
  // Computed/related fields from serializer
  compound_count?: number;
  children?: Target[];
}

export interface Compound {
  id: string;  // UUID
  reg_number: number;
  formatted_id: string;  // NCL-XXXXXXXX
  target: string;  // UUID
  target_name?: string;  // From serializer

  // Chemistry
  smiles: string;
  rdkit_smiles: string | null;
  inchi: string | null;
  molecular_weight: number | null;
  stereo_comment: StereoComment;

  // Provenance
  supplier: string | null;  // UUID
  supplier_name?: string;  // From serializer
  supplier_ref: string | null;
  labbook_number: number | null;
  page_number: number | null;
  compound_number: number;
  barcode: string | null;

  // Audit
  registered_by: number | null;  // User ID
  registered_by_email?: string;  // From serializer
  legacy_registered_by: string | null;
  registered_at: string;  // ISO datetime
  modified_at: string | null;
  comments: string | null;

  // Files
  svg_file: string | null;  // URL

  // Related counts
  batch_count?: number;
}

export type StereoComment =
  | 'unset'
  | 'achiral'
  | 'racemic'
  | 'single_unknown'
  | 'r_enantiomer'
  | 's_enantiomer'
  | 'non_racemic_mixture'
  | 'four_diastereomers'
  | 'two_diastereomers'
  | 'single_diastereomer_unknown'
  | 'rr_diastereomer'
  | 'rs_diastereomer'
  | 'sr_diastereomer'
  | 'ss_diastereomer'
  | 'epimer_mixture'
  | 'ez_mixture'
  | 'e_isomer'
  | 'z_isomer';

export interface Batch {
  id: string;  // UUID
  compound: string;  // UUID
  compound_formatted_id?: string;  // From serializer
  batch_number: number;

  // Provenance
  supplier: string | null;  // UUID
  supplier_name?: string;  // From serializer
  supplier_ref: string | null;
  labbook_number: number | null;
  page_number: number | null;

  // Properties
  amount: string | null;  // Decimal as string
  salt_code: string | null;
  molecular_weight: number | null;

  // Metadata
  registered_at: string;  // ISO datetime
  comments: string | null;

  // Related counts
  qc_file_count?: number;
}

export interface BatchQCFile {
  id: string;  // UUID
  batch: string;  // UUID
  batch_display?: string;  // From serializer
  file: string;  // URL
  filename: string | null;
  comments: string | null;
  uploaded_at: string | null;  // ISO datetime
}

export interface CompoundTemplate {
  id: string;  // UUID
  target: string | null;  // UUID
  target_name?: string;  // From serializer
  name: string;
  mol2d: string;  // 2D MOL block
  svg_file: string | null;  // URL
}

// =============================================================================
// Assay Models
// =============================================================================

export interface DilutionSeries {
  id: string;  // UUID
  concentrations: number[];
  unit: 'nM' | 'uM' | 'mM';
}

export type AnalysisMethod =
  | 'hill_langmuir'
  | 'hill_langmuir_fix_hill'
  | 'hill_langmuir_fix_hill_minmax'
  | 'hill_langmuir_fix_minmax'
  | 'ms_intact'
  | 'table_of_values';

export interface Protocol {
  id: string;  // UUID
  name: string;
  analysis_method: AnalysisMethod;
  preferred_dilutions: string | null;  // UUID
  created_by: number | null;
  created_at: string;  // ISO datetime
  comments: string | null;
}

export interface Assay {
  id: string;  // UUID
  protocol: string;  // UUID
  protocol_name?: string;  // From serializer
  target: string | null;  // UUID
  target_name?: string;  // From serializer
  data_file: string;  // URL
  data_filename: string | null;
  labbook_number: number | null;
  page_number: number | null;
  created_by: number | null;
  created_at: string;  // ISO datetime
  comments: string | null;

  // Related counts
  data_series_count?: number;
}

export interface DataSeries {
  id: string;  // UUID
  assay: string;  // UUID
  compound: string | null;  // UUID
  compound_formatted_id?: string;  // From serializer
  compound_name: string | null;

  // Position in source
  row: number;
  start_column: number;
  end_column: number;

  // Data
  dilution_series: string | null;  // UUID
  extracted_data: Record<string, any>;
  skip_points: number[];

  // Results
  analysis: string | null;  // UUID
  analysis_result?: AnalysisResult;  // Nested from serializer
  plot_image: string | null;  // URL to fitted curve plot image
}

export type AnalysisStatus = 'valid' | 'invalid' | 'unassigned';

export interface AnalysisResult {
  id: string;  // UUID
  status: AnalysisStatus;
  results: Record<string, any>;  // {EC50, Hill, minVal, maxVal, KPI, ...}
  kpi_value?: any;  // Computed property
}

export type HypothesisStatus = 'pending' | 'rejected' | 'chemistry' | 'shelved' | 'made';

export interface Hypothesis {
  id: string;  // UUID
  target: string;  // UUID
  target_name?: string;  // From serializer
  parent: string | null;  // UUID of parent hypothesis
  smiles: string;
  rationale: string | null;
  model_url: string | null;
  status: HypothesisStatus;
  product_compound: string | null;  // UUID
  product_compound_formatted_id?: string;  // From serializer
  completion_notes: string | null;
  svg_file: string | null;  // URL
  created_at: string;  // ISO datetime
  updated_at: string;  // ISO datetime
}

// =============================================================================
// API Response Types
// =============================================================================

export interface PaginatedResponse<T> {
  count: number;
  next: string | null;
  previous: string | null;
  results: T[];
}
