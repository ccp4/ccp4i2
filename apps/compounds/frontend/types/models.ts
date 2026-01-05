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
