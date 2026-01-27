/**
 * Campaign (Fragment Screening) type definitions.
 *
 * These map to the Django models ProjectGroup and ProjectGroupMembership
 * in server/ccp4i2/db/models.py.
 *
 * Campaigns organize multiple projects for batch crystallographic processing,
 * typically for fragment screening where a parent project provides reference
 * coordinates and FreeR flags, and member projects each represent a dataset
 * soaked with a different compound.
 */

import { Project, File as DbFile, Job } from "./models";

/**
 * Browser native File type (for uploads).
 * Aliased to avoid confusion with the DbFile model type.
 */
export type NativeFile = globalThis.File;

// =============================================================================
// Core Models
// =============================================================================

export type ProjectGroupType = "general_set" | "fragment_set";

export type MembershipType = "parent" | "member";

export interface ProjectGroup {
  id: number;
  name: string;
  type: ProjectGroupType;
}

export interface ProjectGroupMembership {
  id: number;
  group: number;
  project: number;
  type: MembershipType;
}

// =============================================================================
// Extended Types for API Responses
// =============================================================================

/**
 * ProjectGroup with nested memberships (from retrieve endpoint)
 */
export interface ProjectGroupDetail extends ProjectGroup {
  memberships: ProjectGroupMembership[];
}

/**
 * Job summary for a project in a campaign
 */
export interface JobSummary {
  total: number;
  finished: number;
  failed: number;
  running: number;
  pending: number;
}

/**
 * Key Performance Indicators from job results
 */
export interface ProjectKPIs {
  RFactor?: number;
  RFree?: number;
  highResLimit?: number;
  [key: string]: string | number | undefined;
}

/**
 * Minimal job info for display in campaign view
 */
export interface CampaignJobInfo {
  id: number;
  uuid: string;
  number: string;
  task_name: string;
  status: number;
  title: string;
}

/**
 * Member project with job summary, jobs list, and KPIs (from member_projects endpoint)
 * Note: We use Omit to exclude the base Project.jobs (number[]) as the API returns
 * detailed job info (CampaignJobInfo[]) for campaign views.
 */
export interface MemberProjectWithSummary extends Omit<Project, "jobs"> {
  job_summary: JobSummary;
  /** Detailed job info for campaign display (replaces base Project.jobs: number[]) */
  jobs: CampaignJobInfo[];
  kpis: ProjectKPIs;
}

/**
 * Parent files response (from parent_files endpoint)
 */
export interface ParentFilesResponse {
  coordinates: DbFile[];
  freer: DbFile[];
}

// =============================================================================
// Batch Import Types
// =============================================================================

export type BatchFileStatus =
  | "idle"
  | "creating_project"
  | "creating_job"
  | "uploading_coords"
  | "uploading_reflections"
  | "fetching_smiles"
  | "setting_params"
  | "queuing"
  | "done"
  | "error";

/**
 * Represents a file in the batch import queue with extracted metadata
 */
export interface BatchFileItem {
  /** The browser File object for upload */
  file: NativeFile;
  /** Visit code extracted from filename */
  visit?: string;
  /** Crystal name extracted from filename */
  crystal?: string;
  /** Registry ID (NCL ID) extracted from filename */
  nclId?: string;
  /** Processing method extracted from filename */
  processing?: string;
  /** Current processing status */
  status: BatchFileStatus;
  /** Error message if status is 'error' */
  error?: string;
  /** SMILES string if compound was found */
  smiles?: string;
  /** Created project ID after successful import */
  projectId?: number;
  /** Created job ID after successful import */
  jobId?: number;
}

/**
 * Regex pattern for extracting NCL ID from dataset filenames.
 * Case-insensitive, matches NCL-XXXXX or NCL_XXXXX anywhere in the string.
 * Accepts variable-length numbers (1-8 digits), with or without zero-padding.
 *
 * Examples that match:
 * - dls12092024_x0203_NCL-00028921_xia_3dii
 * - dls041225_CDK4D1_12_NCL-00030882_AutoProcST
 * - mx12345-1_xtal1_ncl-28921_autoproc.mtz
 * - test_NCL_12345_data
 */
export const NCL_ID_PATTERN = /ncl[-_](\d{1,8})/i;

/**
 * Legacy strict pattern for backwards compatibility.
 * Only matches the original expected format.
 * @deprecated Use NCL_ID_PATTERN for more flexible matching
 */
export const DATASET_FILENAME_PATTERN =
  /^(?<visit>[^_]+)_(?<crystal>[^_]+)_NCL[-_](?<nclId>\d{8})_(?<processing>.+)$/;

/**
 * Parse a filename to extract campaign metadata.
 * The NCL ID is the critical piece - extracted case-insensitively from anywhere in the filename.
 * Other fields (visit, crystal, processing) are extracted on a best-effort basis.
 */
export function parseDatasetFilename(filename: string): {
  visit?: string;
  crystal?: string;
  nclId?: string;
  processing?: string;
} {
  // Remove extension
  const baseName = filename.replace(/\.[^.]+$/, "");

  // Extract NCL ID - this is the critical piece (case-insensitive, anywhere in string)
  const nclMatch = NCL_ID_PATTERN.exec(baseName);
  const nclId = nclMatch ? nclMatch[1] : undefined;

  if (!nclId) {
    return {};
  }

  // Try to extract other metadata by splitting around the NCL pattern
  const nclStart = baseName.toLowerCase().indexOf("ncl");
  const beforeNcl = baseName.substring(0, nclStart).replace(/_$/, ""); // Remove trailing underscore
  // Remove the NCL-XXXXX part and any following underscore
  const afterNclMatch = baseName.substring(nclStart).match(/^ncl[-_]\d+_?(.*)/i);
  const afterNcl = afterNclMatch ? afterNclMatch[1] : "";

  // Before NCL: split by underscore, first is visit, rest is crystal
  const beforeParts = beforeNcl.split("_").filter(Boolean);
  const visit = beforeParts[0] || undefined;
  const crystal = beforeParts.slice(1).join("_") || undefined;

  // After NCL: everything remaining is processing
  const processing = afterNcl || undefined;

  return {
    visit,
    crystal,
    nclId,
    processing,
  };
}

// =============================================================================
// PANDDA Export Types
// =============================================================================

export interface PanddaExportProject {
  projectName: string;
  dimpleJobId: number | null;
  acedrgJobId: number | null;
}

export interface PanddaExportRequest {
  projects: string[];
  dimples: string[];
  acedrgs: string[];
}

export interface PanddaExportStatus {
  Done: boolean;
  Result?: {
    size: number;
  };
}

// =============================================================================
// Form/Dialog Types
// =============================================================================

export interface CreateCampaignFormData {
  name: string;
  parentProjectId?: number;
}

export interface ImportCoordsFormData {
  file: NativeFile;
}

export interface ImportFreeRFormData {
  fSigFFile: NativeFile;
  freeRFile?: NativeFile;
  maxResolution?: number;
}
