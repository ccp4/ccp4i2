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
  fitting_method?: string | null;
  fitting_method_name?: string;
  plate_layout?: Partial<PlateLayout> | null;
  fitting_parameters?: Record<string, any> | null;
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

/**
 * Plate Layout Configuration types.
 * Used to define the physical arrangement of assay plates.
 */

export type PlateFormat = 24 | 96 | 384 | 1536;

export type ReplicatePattern =
  | 'adjacent_rows'      // A1,B1 are replicates of compound 1
  | 'adjacent_columns'   // A1,A2 are replicates of compound 1
  | 'grouped_rows'       // A1-A12, B1-B12 are replicates of compounds 1-12
  | 'interleaved_rows';  // Alternating rows for same compound

export type ControlPlacement =
  | 'edge_columns'       // Controls in dedicated columns at plate edges
  | 'edge_rows'          // Controls in dedicated rows at top/bottom
  | 'per_compound';      // Controls interspersed with each compound's data strip

export type CompoundSourceType =
  | 'row_order'          // Compounds assigned by row order in sample region
  | 'column_header'      // Compound IDs from column headers in import file
  | 'row_header'         // Compound IDs from row labels
  | 'adjacent_column'    // Compound IDs in column immediately after data region
  | 'plate_map_file'     // Separate plate map defines compound positions
  | 'explicit_wells';    // Per-well compound assignment in layout

export type DilutionDirection = 'horizontal' | 'vertical';

export interface WellRegion {
  columns: number[];      // 1-indexed column numbers
  rows: string[];         // Row letters: A, B, C, etc.
}

export interface ControlsConfig {
  placement: ControlPlacement;
  max: WellRegion;
  min: WellRegion;
}

/**
 * Strip layout configuration for "per_compound" control placement.
 *
 * Example: 24-column plate with pattern [min×2][data×8][max×2] repeated twice per row
 * - strip_width: 12 (total columns per compound including controls)
 * - min_wells: 2 (control wells at start of strip)
 * - data_wells: 8 (concentration points)
 * - max_wells: 2 (control wells at end of strip)
 * - strips_per_row: 2 (A1-A12 and A13-A24 are replicates of same compound)
 */
export interface StripLayoutConfig {
  strip_width: number;        // Total columns per strip (min + data + max)
  min_wells: number;          // Number of min control wells at start
  data_wells: number;         // Number of data/concentration wells
  max_wells: number;          // Number of max control wells at end
  strips_per_row: number;     // Number of replicate strips per row
}

export interface SampleRegion {
  start_column: number;
  end_column: number;
  start_row: string;
  end_row: string;
}

export interface DilutionConfig {
  direction: DilutionDirection;
  num_concentrations: number;
}

export interface ReplicateConfig {
  count: number;
  pattern: ReplicatePattern;
}

export interface CompoundSourceConfig {
  type: CompoundSourceType;
  id_column?: string;     // Column name for column_header type
}

/**
 * Spreadsheet origin configuration.
 * Defines where the plate data starts in imported Excel files.
 */
export interface SpreadsheetOrigin {
  column: string;  // Excel column letter (e.g., "A", "B", "AA")
  row: number;     // 1-indexed row number
}

export interface PlateLayout {
  plate_format: PlateFormat;
  controls: ControlsConfig;
  sample_region: SampleRegion;
  dilution: DilutionConfig;
  replicate: ReplicateConfig;
  compound_source: CompoundSourceConfig;
  // For strip-based layouts with embedded controls
  strip_layout?: StripLayoutConfig;
  // Where plate data starts in imported spreadsheets
  spreadsheet_origin?: SpreadsheetOrigin;
}
