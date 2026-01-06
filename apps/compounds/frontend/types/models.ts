/**
 * Compounds Registry type definitions.
 * Maps to Django models in apps/compounds/registry.
 */

export interface Supplier {
  id: string;
  name: string;
  initials: string | null;
}

export interface Target {
  id: string;
  name: string;
  parent: string | null;
  created_at: string;
  compound_count?: number;
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

export interface Compound {
  id: string;
  reg_number: number;
  formatted_id: string;
  target: string;
  target_name?: string;
  smiles: string;
  rdkit_smiles: string | null;
  inchi: string | null;
  molecular_weight: number | null;
  stereo_comment: StereoComment;
  supplier: string | null;
  supplier_name?: string;
  supplier_ref: string | null;
  labbook_number: number | null;
  page_number: number | null;
  compound_number: number;
  barcode: string | null;
  registered_by: number | null;
  registered_by_email?: string;
  legacy_registered_by: string | null;
  registered_at: string;
  modified_at: string | null;
  comments: string | null;
  svg_file: string | null;
  batch_count?: number;
}

export interface Batch {
  id: string;
  compound: string;
  compound_formatted_id?: string;
  batch_number: number;
  supplier: string | null;
  supplier_name?: string;
  supplier_ref: string | null;
  labbook_number: number | null;
  page_number: number | null;
  amount: string | null;
  salt_code: string | null;
  molecular_weight: number | null;
  registered_at: string;
  comments: string | null;
  qc_file_count?: number;
}

export interface BatchQCFile {
  id: string;
  batch: string;
  batch_display?: string;
  file: string;
  filename: string | null;
  comments: string | null;
  uploaded_at: string | null;
}

export interface PaginatedResponse<T> {
  count: number;
  next: string | null;
  previous: string | null;
  results: T[];
}

/**
 * Compounds Assays type definitions.
 * Maps to Django models in apps/compounds/assays.
 */

export interface DilutionSeries {
  id: string;
  concentrations: number[];
  unit: 'nM' | 'uM' | 'mM';
  display_name?: string;
}

export type AnalysisMethod =
  | 'hill_langmuir'
  | 'hill_langmuir_fix_hill'
  | 'hill_langmuir_fix_hill_minmax'
  | 'hill_langmuir_fix_minmax'
  | 'ms_intact'
  | 'table_of_values';

export interface Protocol {
  id: string;
  name: string;
  analysis_method: AnalysisMethod;
  pherastar_table: string | null;
  preferred_dilutions: string | null;
  preferred_dilutions_display?: string;
  created_by: number | null;
  created_by_email?: string;
  created_at: string;
  comments: string | null;
  assays_count?: number;
}

export interface Assay {
  id: string;
  protocol: string;
  protocol_name?: string;
  target: string | null;
  target_name?: string;
  data_file: string;
  data_filename?: string;
  labbook_number: number | null;
  page_number: number | null;
  created_by: number | null;
  created_by_email?: string;
  created_at: string;
  comments: string | null;
  data_series_count?: number;
  data_series?: DataSeries[];
}

export type AnalysisStatus = 'valid' | 'invalid' | 'unassigned';

export interface AnalysisResult {
  id: string;
  status: AnalysisStatus;
  results: Record<string, any>;
  kpi_value?: number | null;
}

export interface DataSeries {
  id: string;
  assay: string;
  compound: string | null;
  compound_formatted_id?: string;
  compound_name: string | null;
  row: number;
  start_column: number;
  end_column: number;
  dilution_series?: DilutionSeries;
  extracted_data: number[];  // Response values aligned with dilution_series.concentrations
  skip_points: number[];     // Indices of points to exclude from curve fitting
  analysis?: AnalysisResult;
  analysis_status?: AnalysisStatus;
  analysis_kpi?: number | null;
  svg_file: string | null;
}

export type HypothesisStatus = 'pending' | 'rejected' | 'chemistry' | 'shelved' | 'made';

export interface Hypothesis {
  id: string;
  target: string;
  target_name?: string;
  parent: string | null;
  parent_smiles?: string;
  smiles: string;
  rationale: string | null;
  model_url: string | null;
  status: HypothesisStatus;
  product_compound: string | null;
  product_formatted_id?: string;
  completion_notes: string | null;
  svg_file: string | null;
  created_at: string;
  updated_at: string;
}
