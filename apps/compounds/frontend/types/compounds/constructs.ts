/**
 * Constructs Database type definitions.
 * Maps to Django models in apps/compounds/constructs.
 */

// =============================================================================
// Reference Data Types
// =============================================================================

export interface ExpressionTagType {
  id: string;
  name: string;
  created_at: string;
}

export interface Protease {
  id: string;
  name: string;
  created_at: string;
}

export interface ExpressionTagLocation {
  id: string;
  name: string;
  created_at: string;
}

// =============================================================================
// Construct Project Types
// =============================================================================

export interface ConstructProject {
  id: string;
  name: string;
  parent: string | null;
  parent_name?: string;
  plasmid_count?: number;
  created_at: string;
  updated_at: string;
  created_by: number | null;
  created_by_email?: string;
}

export interface ConstructProjectDetail extends ConstructProject {
  children: ConstructProject[];
}

// =============================================================================
// Protein Types
// =============================================================================

export interface ProteinSynonym {
  id: string;
  name: string;
  protein: string;
  protein_uniprot_id?: string;
  created_at: string;
}

export interface Protein {
  id: string;
  uniprot_id: string;
  synonym_count?: number;
  cassette_count?: number;
  created_at: string;
}

export interface ProteinDetail extends Protein {
  synonyms: ProteinSynonym[];
  updated_at: string;
  created_by: number | null;
  created_by_email?: string;
}

export interface ProteinUse {
  id: string;
  protein: string;
  protein_uniprot_id?: string;
  project: string;
  project_name?: string;
  created_at: string;
}

// =============================================================================
// Cassette Types
// =============================================================================

export interface Cassette {
  id: string;
  protein: string;
  protein_uniprot_id?: string;
  start: number;
  end: number;
  display_name?: string;
  plasmid_count?: number;
  created_at: string;
}

export interface CassetteDetail extends Cassette {
  updated_at: string;
  created_by: number | null;
  created_by_email?: string;
}

// =============================================================================
// Expression Tag Types
// =============================================================================

export interface ExpressionTag {
  id: string;
  expression_tag_type: string;
  expression_tag_type_name?: string;
  protease: string | null;
  protease_name?: string;
  cassette_use: string;
  location: string;
  location_name?: string;
  created_at: string;
}

// =============================================================================
// Cassette Use Types
// =============================================================================

export interface CassetteUse {
  id: string;
  cassette: string;
  cassette_display?: string;
  plasmid: string;
  plasmid_formatted_id?: string;
  expression_tag_count?: number;
  created_at: string;
}

export interface CassetteUseDetail extends CassetteUse {
  alignment_file: string | null;
  expression_tags: ExpressionTag[];
  updated_at: string;
  created_by: number | null;
  created_by_email?: string;
}

// =============================================================================
// Sequencing Result Types
// =============================================================================

export interface SequencingResult {
  id: string;
  cassette_use: string;
  cassette_display?: string;
  plasmid: string;
  plasmid_formatted_id?: string;
  file: string;
  filename?: string;
  created_at: string;
}

// =============================================================================
// Plasmid Types
// =============================================================================

export interface Plasmid {
  id: string;
  ncn_id: number;
  formatted_id: string;
  name: string;
  project: string | null;
  project_name?: string;
  parent: string | null;
  parent_formatted_id?: string;
  cassette_count?: number;
  created_at: string;
  created_by_email?: string;
}

export interface PlasmidDetail extends Plasmid {
  genbank_file: string | null;
  genbank_file_url?: string;
  cassette_uses: CassetteUseDetail[];
  sequencing_results: SequencingResult[];
  updated_at: string;
  created_by: number | null;
}

// =============================================================================
// API Response Types
// =============================================================================

export interface GenbankContentResponse {
  content: string;
}
