// Wire-contract types for legacy CCP4i2 -> Django migration.
// Mirrors docs/LEGACY_MIGRATION_WIRE_CONTRACT.md and the shapes returned by
// server/ccp4i2/db/import_sqlite.py + api/admin_views.py.

export type IssueSeverity = "blocking" | "warning" | "info";

export interface StructureIssue {
  type: string; // nested_project | dest_collision | case_clash | ...
  severity: IssueSeverity;
  mode_scope: "copy" | "schlurp" | "both";
  projects: string[]; // legacy project ids (32-hex)
  detail: string;
  resolution: string | null; // null => needs user action; else what migration will do
}

export interface Structure {
  ok: boolean;
  issues: StructureIssue[];
}

export interface PlanEntry {
  legacy_project_id: string;
  name: string;
  source_dir: string;
  mode: "in_place" | "copy";
  dest: string;
  reason: "copy_files" | "nested" | "in_place";
  source_exists: boolean;
  renamed_to: string | null;
}

export interface ValidateSummary {
  ok: boolean;
  projects_on_disk: string;
  jobs_on_disk: string;
  files_on_disk: string;
  // Informational: original external import sources still resolvable. Not a
  // health signal (files_on_disk already covers the imported copies).
  import_sources_present: string;
  // Missing files split by the preservation contract: top-level job files are
  // guaranteed (a miss is a real violation), sub-job files are not.
  top_level_files_missing: number;
  sub_job_files_missing: number;
  integrity_issues: number;
  data_quality_issues: number;
  structure_issues: number;
  blocking_issues: number;
  plan_summary: {
    in_place: number;
    copy: number;
    copied_due_to_nesting: number;
  };
}

export interface ValidateReport {
  source: string;
  counts: Record<string, number | null>;
  projects: {
    total: number;
    dir_exists: number;
    dir_missing: { project: string; directory: string }[];
    dir_empty: string[];
  };
  jobs: {
    total: number;
    dir_exists: number;
    dir_missing_count: number;
    dir_missing: { job_number: string; task_name: string; expected_dir: string }[];
  };
  files: {
    total: number;
    exists: number;
    missing_count: number;
    missing: {
      filename: string;
      expected_path: string;
      job_number: string;
      top_level: boolean;
    }[];
    top_missing_count: number;
    top_missing: { filename: string; expected_path: string; job_number: string }[];
    sub_missing_count: number;
  };
  imported_files: {
    total: number;
    // Original external import sources still present (informational only).
    source_exists: number;
    source_missing_count: number;
  };
  integrity: { ok: boolean; issues: string[] };
  data_quality: { ok: boolean; issues: string[] };
  structure: Structure;
  plan: PlanEntry[];
  summary: ValidateSummary;
}

export interface ImportResult {
  dry_run: boolean;
  source: string;
  stats: Record<string, number>;
  errors: string[];
  plan: PlanEntry[];
  structure: Structure;
}

export interface ImportStatus {
  projects: { total: number; groups: number; memberships: number; tags: number };
  jobs: {
    total: number;
    value_keys: number;
    float_values: number;
    char_values: number;
    xdata: number;
  };
  files: { total: number; types: number; uses: number; imports: number; exports: number };
}
