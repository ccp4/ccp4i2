/**
 * CCP4i2 service contract — v0 TypeScript types.
 *
 * These types describe the response shapes of the documented stable
 * subset of CCP4i2's REST API (under ``/api/ccp4i2/``). External
 * consumers — Materia, third-party integrators — depend on these
 * shapes; CCP4i2 commits to not removing documented fields without a
 * major version bump.
 *
 * See ``docs/CCP4I2_SERVICE_CONTRACT.md`` (in the ccp4/ccp4i2 repo)
 * for the narrative documentation, stability promises, and
 * out-of-scope surface notes.
 *
 * Endpoints NOT covered by these types are CCP4i2-internal and may
 * change without notice.
 */

// ---------------------------------------------------------------------------
// Common building blocks
// ---------------------------------------------------------------------------

/** ISO 8601 timestamp string (e.g. "2026-04-29T13:21:24.854154Z"). */
export type Iso8601 = string;

/** UUID v4 string. */
export type UuidV4 = string;

/**
 * Canonical error response shape returned by CCP4i2's auth middleware
 * for 401 / 403 responses. The renderer's auth-error event handler and
 * Materia's equivalent both pattern-match against this.
 */
export interface AuthErrorBody {
  success: false;
  error: string;
}

// ---------------------------------------------------------------------------
// /version/ — server version info (public; no auth required)
// ---------------------------------------------------------------------------

/**
 * Server version info, returned by ``GET /api/ccp4i2/version/``.
 * Consumers use this for compatibility-banner rendering at app load.
 */
export interface VersionInfo {
  /** CCP4i2 software version string. */
  version: string;
  /** CCP4 suite version baked into the deployment. */
  ccp4_version?: string;
  /** Build timestamp (ISO 8601) of the running image, if available. */
  build_timestamp?: Iso8601;
  /** Short git commit hash of the running image, if available. */
  git_commit?: string;
}

// ---------------------------------------------------------------------------
// /health/ — server health (public; no auth required)
// ---------------------------------------------------------------------------

/**
 * Health check response, returned by ``GET /api/ccp4i2/health/``.
 * Used by container orchestration (Azure Container Apps, k8s) and by
 * the renderer's connectivity diagnostics.
 */
export interface HealthStatus {
  status: "ok" | "degraded" | "down";
  /** Free-text detail; not parsed by consumers. */
  detail?: string;
}

// ---------------------------------------------------------------------------
// Projects
// ---------------------------------------------------------------------------

/**
 * Lightweight project shape returned by ``GET /api/ccp4i2/projects/``.
 * Excludes ``directory`` to keep list responses small. Tags are
 * inlined with id+text only (no reverse projects relation).
 */
export interface ProjectListItem {
  id: number;
  uuid: UuidV4;
  name: string;
  creation_time: Iso8601;
  last_access: Iso8601;
  tags: ProjectTagSummary[];
}

/**
 * Full project shape returned by ``GET /api/ccp4i2/projects/{id}/``
 * and the create/update endpoints. Fields beyond those listed here are
 * subject to change.
 */
export interface Project {
  id: number;
  uuid: UuidV4;
  name: string;
  description: string;
  directory: string;
  creation_time: Iso8601;
  creation_user: string;
  creation_host: string;
  last_access: Iso8601;
  last_job_number: number;
  follow_from_job: number | null;
  i1_project_name: string;
  i1_project_directory: string;
  tags: ProjectTagDetail[];
}

export interface ProjectTagSummary {
  id: number;
  text: string;
  parent: number | null;
}

export interface ProjectTagDetail extends ProjectTagSummary {
  /** Project IDs this tag is applied to. */
  projects: number[];
}

// ---------------------------------------------------------------------------
// Jobs
// ---------------------------------------------------------------------------

/**
 * Job-status enum. Values match the integer codes on the Django Job
 * model. Numeric stability is part of the contract — consumers may
 * pattern-match on these values.
 */
export const JobStatus = {
  Unknown: 0,
  Pending: 1,
  Queued: 2,
  Running: 3,
  Interrupted: 4,
  Failed: 5,
  Finished: 6,
  RunningRemotely: 7,
  FileHolder: 8,
  ToDelete: 9,
  Unsatisfactory: 10,
} as const;
export type JobStatusValue = (typeof JobStatus)[keyof typeof JobStatus];

/** Job-evaluation enum. Numeric stability is part of the contract. */
export const JobEvaluation = {
  Unknown: 0,
  Best: 1,
  Good: 2,
  Rejected: 3,
} as const;
export type JobEvaluationValue =
  (typeof JobEvaluation)[keyof typeof JobEvaluation];

/**
 * Job shape returned by ``GET /api/ccp4i2/jobs/`` (list) and
 * ``GET /api/ccp4i2/jobs/{id}/`` (detail). The shape is identical for
 * both endpoints — all fields below are guaranteed.
 */
export interface Job {
  id: number;
  uuid: UuidV4;
  project: number;
  parent: number | null;
  number: string;
  title: string;
  status: JobStatusValue;
  evaluation: JobEvaluationValue;
  comments: string;
  creation_time: Iso8601;
  finish_time: Iso8601 | null;
  task_name: string;
  process_id: number | null;
}

/**
 * Active jobs summary, returned by ``GET /api/ccp4i2/active_jobs/``.
 * Lightweight — used for per-second polling of the running set.
 */
export interface ActiveJobsResponse {
  jobs: Pick<Job, "id" | "uuid" | "title" | "task_name" | "status">[];
}

// ---------------------------------------------------------------------------
// Files
// ---------------------------------------------------------------------------

/** File-directory enum. Numeric stability is part of the contract. */
export const FileDirectory = {
  JobDir: 1,
  ImportDir: 2,
} as const;
export type FileDirectoryValue =
  (typeof FileDirectory)[keyof typeof FileDirectory];

/**
 * File shape returned by ``GET /api/ccp4i2/files/{id}/`` and embedded
 * in project / job detail responses.
 */
export interface File {
  id: number;
  uuid: UuidV4;
  name: string;
  directory: FileDirectoryValue;
  /** FK to FileType (the type's ``name`` is the primary key). */
  type: string;
  sub_type: number | null;
  content: number | null;
  annotation: string;
  job: number | null;
  job_param_name: string;
  /** Computed: full filesystem path, derivable from directory+name+job. */
  path: string | null;
}

// ---------------------------------------------------------------------------
// Project groups (campaigns)
// ---------------------------------------------------------------------------

/**
 * Project-group kind. Membership of this union is contracted: adding a
 * new value is a minor change; removing or renaming an existing value
 * is a major change.
 */
export type ProjectGroupType = "general_set" | "fragment_set";

/** Membership type within a project group. */
export type MembershipType = "parent" | "member";

/**
 * Saved view-state landmark within a campaign. Used by the Moorhen
 * viewer to navigate between binding-site origins of interest.
 */
export interface CampaignSite {
  /** Display name for the site. */
  name: string;
  /** View origin coordinates [x, y, z]. */
  origin: [number, number, number];
  /** View orientation quaternion [x, y, z, w] (optional). */
  quat?: [number, number, number, number];
  /** Zoom level (optional). */
  zoom?: number;
}

/**
 * Project-group shape returned by ``GET /api/ccp4i2/projectgroups/``
 * (list) and embedded in detail responses.
 */
export interface ProjectGroup {
  id: number;
  name: string;
  type: ProjectGroupType;
  /** Saved binding sites for Moorhen viewer navigation. */
  sites?: CampaignSite[];
}

/**
 * A membership tying a Project to a ProjectGroup.
 */
export interface ProjectGroupMembership {
  id: number;
  group: number;
  project: number;
  type: MembershipType;
}

/**
 * Full project-group shape returned by
 * ``GET /api/ccp4i2/projectgroups/{id}/``, including its memberships
 * inline. Memberships are bounded per-group (not an unbounded child
 * collection — they're inlined per the *Nested children* convention).
 */
export interface ProjectGroupDetail extends ProjectGroup {
  memberships: ProjectGroupMembership[];
}

/**
 * Counts of jobs in a member project, broken down by status. Inlined
 * into MemberProjectWithSummary.
 */
export interface JobSummary {
  total: number;
  finished: number;
  failed: number;
  running: number;
  pending: number;
}

/**
 * Key Performance Indicators extracted from job results. Documented
 * keys are stable; the dictionary may carry additional task-specific
 * keys, which consumers should treat as opaque.
 */
export interface ProjectKPIs {
  RFactor?: number;
  RFree?: number;
  highResLimit?: number;
  [key: string]: string | number | undefined;
}

/**
 * Minimal per-job info inlined into MemberProjectWithSummary for
 * campaign-dashboard rendering.
 */
export interface CampaignJobInfo {
  id: number;
  uuid: UuidV4;
  number: string;
  task_name: string;
  status: JobStatusValue;
  title: string;
}

/**
 * Member project enriched with job summary, per-job info, and KPIs.
 * Returned by ``GET /api/ccp4i2/projectgroups/{id}/member_projects/``.
 * This is an enrichment-bearing endpoint, not a filtered list of
 * Project — the embedded summaries are computed server-side.
 *
 * Uses Omit<Project, "jobs"> defensively: if Project ever gains a
 * "jobs" field, the enriched shape stays self-consistent.
 */
export interface MemberProjectWithSummary extends Omit<Project, "jobs"> {
  job_summary: JobSummary;
  jobs: CampaignJobInfo[];
  kpis: ProjectKPIs;
}

/**
 * Coordinate and FreeR file references from a campaign's parent
 * project. Returned by ``GET /api/ccp4i2/projectgroups/{id}/parent_files/``.
 * Used by member projects to inherit reference inputs.
 */
export interface ParentFilesResponse {
  coordinates: File[];
  freer: File[];
}
