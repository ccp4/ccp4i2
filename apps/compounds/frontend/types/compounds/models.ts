/**
 * Compounds Registry type definitions.
 * Maps to Django models in apps/compounds/registry.
 */

// =============================================================================
// User & Authorization Types
// =============================================================================

/**
 * User role levels (matches backend UserProfile.ROLE_CHOICES)
 */
export type UserRole = 'user' | 'contributor' | 'admin';

/**
 * User profile from the API
 */
export interface UserProfile {
  role: UserRole;
  is_platform_admin: boolean;
  legacy_username: string;
  legacy_display_name: string;
  imported_at: string | null;
  first_login_at: string | null;
  last_seen_at: string | null;
}

/**
 * Current user info from /api/users/me/
 */
export interface CurrentUser {
  id: number;
  username: string;
  email: string;
  first_name: string;
  last_name: string;
  display_name: string;
  is_admin: boolean;
  role: UserRole;
  operating_level: UserRole;
  can_contribute: boolean;
  can_administer: boolean;
  profile: UserProfile;
}

/**
 * Operating level info from /api/users/me/operating-level/
 */
export interface OperatingLevelInfo {
  operating_level: UserRole;
  role: UserRole;
  available_levels: UserRole[];
}

// =============================================================================
// Registry Types
// =============================================================================

export interface Supplier {
  id: string;
  name: string;
  initials: string | null;
  user: number | null;
  is_current_user?: boolean;
  compound_count?: number;
  batch_count?: number;
}

/** Concentration display mode for KPI values */
export type ConcentrationDisplayMode = 'natural' | 'nM' | 'uM' | 'mM' | 'pConc';

/** Saved aggregation view configuration stored on a Target */
export interface SavedAggregationView {
  protocol_names: string[];
  compound_search: string;
  output_format: 'compact' | 'medium' | 'long';
  aggregations: ('geomean' | 'count' | 'stdev' | 'list')[];
  status: 'valid' | 'invalid' | 'unassigned' | '';
  /** Concentration display mode (default: 'natural') */
  concentration_display?: ConcentrationDisplayMode;
}

export interface Target {
  id: string;
  name: string;
  parent: string | null;
  created_at: string;
  compound_count?: number;
  assay_count?: number;
  has_recent_compounds?: boolean;
  has_recent_assays?: boolean;
  latest_activity?: string | null;
  image?: string | null;
  saved_aggregation_view?: SavedAggregationView | null;
}

/**
 * Dashboard types for target landing page.
 */

export interface DashboardCompound {
  id: string;
  formatted_id: string;
  smiles: string;
  registered_at: string;
  molecular_weight: number | null;
}

export interface DashboardAssay {
  id: string;
  protocol_name: string;
  created_at: string;
  data_series_count: number;
}

export interface DashboardProject {
  id: number;
  name: string;
  last_access: string;
  job_count: number;
  matching_compound_ids: string[];
}

export interface TargetDashboard extends Target {
  recent_compounds: DashboardCompound[];
  recent_assays: DashboardAssay[];
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
  /** Alternative names/identifiers for this compound (supplier codes, abbreviations, etc.) */
  aliases?: string[];
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

export interface ProtocolDocument {
  id: string;
  protocol: string;
  file: string;
  filename: string | null;
  created_by: number | null;
  created_by_email?: string;
  created_at: string;
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

export interface FittingMethod {
  id: string;
  name: string;
  slug: string;
  version: string;
  description: string;
  is_active: boolean;
  is_builtin: boolean;
}

export interface FittingMethodDetail extends FittingMethod {
  script: string;
  input_schema: Record<string, unknown>;
  output_schema: Record<string, unknown>;
  created_at: string;
  updated_at: string;
}

/**
 * Reusable plate layout configuration.
 * Maps to Django PlateLayout model in apps/compounds/assays.
 */
export interface PlateLayoutRecord {
  id: string;
  name: string;
  description?: string | null;
  config: PlateLayout;  // The actual plate layout configuration
  plate_format?: PlateFormat;  // Convenience field derived from config
  is_builtin: boolean;
  created_by: number | null;
  created_by_email?: string;
  created_at: string;
  updated_at: string;
  protocols_count?: number;  // Number of protocols using this layout
}

/**
 * Import type - what kind of data is being imported.
 * This replaces the legacy analysis_method for categorizing data.
 */
export type ImportType =
  | 'raw_data'        // Raw dose-response data requiring curve fitting
  | 'ms_intact'       // MS-Intact pre-analyzed import
  | 'table_of_values' // Pre-analyzed table import
  | 'pharmaron_adme'; // Pharmaron ADME import

/**
 * @deprecated Use ImportType instead. Kept for backward compatibility.
 */
export type AnalysisMethod =
  | 'hill_langmuir'
  | 'hill_langmuir_fix_hill'
  | 'hill_langmuir_fix_hill_minmax'
  | 'hill_langmuir_fix_minmax'
  | 'ms_intact'
  | 'table_of_values'
  | 'pharmaron_adme';

/**
 * Tight-binding analysis parameters for Wang equation fitting.
 * Required when using the 'tight-binding-wang' fitting method.
 */
export interface TightBindingParameters {
  protein_conc: number;   // Total protein [P]t in nM
  ligand_conc: number;    // Total labeled ligand [L]t in nM
  ligand_kd: number;      // Kd of labeled ligand in nM
}

/**
 * Validation rules for curve fitting quality flags.
 * Determines which flags should cause automatic invalidation.
 */
export interface ValidationRules {
  invalidating_flags?: string[];
}

/**
 * Fitting parameters that can be stored on a Protocol.
 * Different fitting methods use different subsets of these parameters.
 */
export interface FittingParameters {
  // Standard 4PL constraint parameters
  fix_hill?: number | null;    // Specific value to fix Hill coefficient to (e.g., 1.0)
  fix_top?: boolean | null;    // If true, use control max as fixed top asymptote
  fix_bottom?: boolean | null; // If true, use control min as fixed bottom asymptote
  // Tight-binding parameters (for Wang equation)
  protein_conc?: number;
  ligand_conc?: number;
  ligand_kd?: number;
  // Validation rules for curve fit quality flags
  validation_rules?: ValidationRules;
}

export interface Protocol {
  id: string;
  name: string;
  import_type: ImportType;
  /** @deprecated Use import_type instead */
  analysis_method: AnalysisMethod;
  fitting_method?: string | null;
  fitting_method_name?: string;
  target?: string | null;
  target_name?: string;
  // Plate layout is now a FK to PlateLayoutRecord
  plate_layout?: string | null;          // FK ID to PlateLayoutRecord
  plate_layout_name?: string;            // Denormalized name for display
  plate_layout_config?: PlateLayout | null;  // Denormalized config for convenience
  fitting_parameters?: FittingParameters | null;
  preferred_dilutions: string | null;
  preferred_dilutions_display?: string;
  created_by: number | null;
  created_by_email?: string;
  created_at: string;
  comments: string | null;
  assays_count?: number;
  has_recent_assays?: boolean;
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
  extracted_data: number[];  // Legacy: response values only. New: [min_control, resp1, ..., respN, max_control]
  skip_points: number[];     // Indices of points to exclude from curve fitting
  analysis?: AnalysisResult;
  analysis_status?: AnalysisStatus;
  analysis_kpi?: number | null;
  plot_image: string | null;
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
  | 'interleaved_rows'   // Alternating rows for same compound
  | 'explicit';          // Each strip named explicitly; replicates inferred from matching names

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
  /**
   * Absolute Excel row (1-indexed) where compound names start.
   * For stripe layouts, compound names are read from the left-most control column
   * of each stripe at this row. One compound name per data row.
   */
  compound_name_row?: number;
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
