# Legacy Migration — Wire Contract (UI ↔ REST)

The handshake between the migration wizard
([migrate-legacy-content.tsx]) and the two Django admin endpoints. Everything
downstream (§2.5 structural checks, the importer, the wizard steps) is built to
satisfy *this*. See `docs/LEGACY_MIGRATION_UI_PLAN.md` for rationale.

**Contract stance:** the endpoints already exist and return most of this today.
Fields are tagged:
- **(current)** — already returned by
  `SQLiteValidator.run()` / `SQLiteImporter.summary()` verbatim. Do not change.
- **(new)** — must be added (the structural checks, the per-project plan, and the
  richer per-issue shape). Additive only — no current field is renamed or removed,
  so the UI can be written against the full contract before the backend catches
  up, and old CLI callers keep working.

All requests are `multipart/form-data` (so a `.sqlite` upload and scalar fields
share one body), admin-gated (`IsAdminUser`). Base path `/api/proxy/ccp4i2/`.

---

## Shared vocabulary

### `Issue` (new) — one structural/quality problem the UI renders

```jsonc
{
  "type": "nested_project",        // enum, see §Issue types
  "severity": "blocking" | "warning" | "info",
  "mode_scope": "copy" | "schlurp" | "both",  // when this issue actually bites
  "projects": ["<legacy_project_id>", ...],   // legacy ProjectIDs involved (32-hex)
  "detail": "Project 'B' lies under 'A/CCP4_JOBS/job_3'",  // human sentence
  "resolution": "B will be copied to CCP4X_PROJECTS/B and repointed"
                                    // null if user action required; else what migration will do
}
```

- `severity: "blocking"` + `resolution: null` ⇒ the wizard must not let the user
  proceed until they tick the acknowledgement (which sets
  `allow_structural_issues=true`) or fix the DB. `blocking` **with** a non-null
  `resolution` (e.g. nesting) is auto-handled — surface it, require a confirm, but
  it does not hard-stop.
- `mode_scope` lets the UI hide copy-only issues when the user picked schlurp.

#### Issue `type` enum (new)
| type | §2.5 | severity (default) | mode_scope |
|------|------|--------------------|------------|
| `nested_project` | 1 | blocking (auto-resolved) | both |
| `dest_collision` | 2 | blocking | copy |
| `shared_directory` | 3 | warning | both |
| `not_project_root` | 4 | warning | both |
| `dest_unwritable` | 5 | blocking | copy |
| `insufficient_space` | 5 | blocking | copy |
| `i1_dir_overlap` | 6 | warning | both |
| `case_clash` | 7 | blocking | copy |
| `unicode_clash` | 8 | blocking | copy |
| `shared_subtree` | 9 | blocking-if-under-migrating-dest / else warning | copy |
| `unreadable_source` | 10 | blocking | copy |

### `PlanEntry` (new) — what migration will do to one project

```jsonc
{
  "legacy_project_id": "<32-hex>",
  "name": "beta_lactamase",         // legacy ProjectName
  "source_dir": "/home/u/.CCP4I2/CCP4X_PROJECTS/beta_lactamase",  // post-remap
  "mode": "in_place" | "copy",
  "dest": "/home/u/CCP4X_PROJECTS/beta_lactamase",  // where directory ends up; == source_dir when in_place
  "reason": "copy_files" | "nested" | "in_place",   // why this mode
  "source_exists": true,            // legacy dir present on disk?
  "renamed_to": "beta_lactamase_legacy"  // null unless a name/dir collision forced a rename
}
```

The `mode`/`dest` are decided **by the validator** and consumed unchanged by the
importer, so validate and import never disagree.

---

## 1. `POST admin/validate-sqlite/` — read-only pre-flight

Runs `SQLiteValidator`. Never writes to Django, the sqlite, or disk.

### Request (multipart)
| field | req? | notes |
|-------|------|-------|
| `sqlite_db` | one-of | uploaded `.sqlite` file **OR**… |
| `db_path` | one-of | server-side path, e.g. `~/.CCP4I2/db/database.sqlite` (desktop default) |
| `copy_files` | no | `"true"`/`"false"` — validator needs it to compute the plan + which issues are in scope. Default `"false"`. |
| `dest_root` | no | override for copy destination; default `CCP4I2_PROJECTS_DIR`. |
| `remap_from`, `remap_to` | no | prefix remap applied before every disk check. |

### Response `200` — the full report

```jsonc
{
  "source": "/home/u/.CCP4I2/db/database.sqlite",   // (current)

  "counts": {                                        // (current) raw table row counts
    "Projects": 12, "Jobs": 340, "Files": 1201, "FileUses": 1800,
    "ImportFiles": 44, "ExportFiles": 3, "XData": 12, "JobKeyValues": 5000,
    "JobKeyCharValues": 800, "Tags": 6, "ProjectTags": 9, "Comments": 0,
    "ProjectComments": 2, "ServerJobs": 0, "FileTypes": 60, "KeyTypes": 40,
    "FileRoles": 2, "ProjectExports": 3, "ProjectImports": 1
    // a table absent in the legacy DB is null
  },

  "projects": {                                      // (current)
    "total": 12, "dir_exists": 11,
    "dir_missing": [ { "project": "old1", "directory": "/gone/old1" } ],
    "dir_empty":   [ "unnamed_proj" ]                // projects with blank directory
  },
  "jobs": {                                          // (current)
    "total": 340, "dir_exists": 338,
    "dir_missing_count": 2,
    "dir_missing": [ { "job_number": "3.1", "task_name": "refmac",
                       "expected_dir": "/.../CCP4_JOBS/job_3/job_1" } ]  // capped at 50 unless verbose
  },
  "files": {                                         // (current)
    "total": 1201, "exists": 1198,
    "missing_count": 3, "orphan_job": 0,
    "missing": [ { "filename": "XYZOUT.pdb",
                   "expected_path": "/.../job_3/XYZOUT.pdb" } ] // capped at 50
  },
  "imported_files": {                                // (current)
    "total": 44, "source_exists": 40,
    "source_missing_count": 4,
    "source_missing": [ "/vanished/data.mtz" ]       // capped at 50
  },
  "integrity": {                                     // (current)
    "ok": true,
    "issues": [ "2 jobs reference non-existent parent jobs" ]  // free-text strings
  },
  "data_quality": {                                  // (current)
    "ok": true,
    "issues": [ "1 projects have NULL creation time" ]         // free-text strings
  },

  "structure": {                                     // (new) §2.5
    "ok": false,                                     // false if any blocking issue not auto-resolved
    "issues": [ /* Issue[] */ ]
  },
  "plan": [ /* PlanEntry[] — one per legacy project */ ],   // (new)

  "summary": {                                       // (current) + (new) fields
    "ok": false,                                     // (current) — now ALSO false on blocking structural issues
    "projects_on_disk": "11/12",                     // (current) "n/total" strings
    "jobs_on_disk": "338/340",                        // (current)
    "files_on_disk": "1198/1201",                     // (current)
    "import_sources_on_disk": "40/44",                // (current)
    "integrity_issues": 1,                            // (current) counts
    "data_quality_issues": 1,                         // (current)
    "structure_issues": 3,                            // (new)
    "blocking_issues": 1,                             // (new) count of severity=blocking && resolution==null
    "plan_summary": { "in_place": 10, "copy": 2, "copied_due_to_nesting": 1 }  // (new)
  }
}
```

### Errors (standardised on `{error: <machine_string>, reason: <sentence>}`)
- `400` `{ "error": "bad_input", "reason": "Must provide sqlite_db file upload or db_path" }`
- `404` `{ "error": "db_not_found", "reason": "Database not found: <path>" }` (from `FileNotFoundError`)
- `500` `{ "error": "internal_error", "reason": "…" }`

(Today these are bare `error` sentence strings — see the tightening note under
§2 Errors; same change applies here.)

**UI contract:** Step 2 renders `summary` chips + the four `*_on_disk` ratios, the
`integrity`/`data_quality` strings, then the `structure.issues` grouped by `type`,
then the `plan` table. **Next→Dry-run is disabled while `summary.blocking_issues >
0` and the acknowledgement box is unticked.**

---

## 2. `POST admin/import-sqlite/` — dry-run and real import

Runs `SQLiteImporter`. Same endpoint for dry-run and commit; `dry_run` decides.

### Request (multipart)
| field | req? | notes |
|-------|------|-------|
| `sqlite_db` / `db_path` | one-of | as above |
| `dry_run` | no | `"true"` runs inside a rolled-back transaction, copies nothing (default `"false"`). |
| `copy_files` | no | default `"false"`. |
| `dest_root` | no | default `CCP4I2_PROJECTS_DIR`. |
| `remap_from`, `remap_to` | no | as above |
| `allow_structural_issues` | no | **(new)** `"true"` only after the user acknowledged blocking issues in Step 2. If false/absent and the run has blocking-unresolved issues ⇒ `400` with `error: "structural_issues_unacknowledged"` (below). |

### Response `200` — import summary

```jsonc
{
  "dry_run": true,                                   // (current)
  "source": "/home/u/.CCP4I2/db/database.sqlite",    // (current)
  "stats": {                                         // (current keys) + (new keys)
    "projects": 10, "projects_skipped": 2, "projects_renamed": 1,   // (current)
    "jobs": 340, "jobs_skipped": 0,
    "files": 1198, "files_skipped": 3,
    "filetypes": 60, "keytypes": 40, "tags": 6,
    "groups": 1, "memberships": 3, "projecttags": 9,
    "fileuses": 1800, "fileimports": 44, "fileexports": 3,
    "jobfloatvalues": 5000, "jobcharvalues": 800, "xdata": 12,
    "projectexports": 3, "projectimports": 1,

    "projects_copied": 2,               // (new) copy-mode / nested
    "projects_in_place": 10,            // (new)
    "projects_nested": 1,              // (new) copied out because nested
    "copy_excluded_subtrees": 1,       // (new) inner dirs skipped during outer copytree
    "projects_copy_skipped": 0,        // (new) dest already existed
    "copy_errors": 0                   // (new)
  },
  "errors": [                                        // (current) — list of strings
    "File 0a1b… : <exception text>"
  ],
  "plan": [ /* PlanEntry[] — what was (or would be) done */ ],  // (new) echoes validate's plan, with post-run renamed_to filled in
  "structure": { "ok": true, "issues": [] }          // (new) recomputed; empty on a clean run
}
```

- HTTP status: `200` when `errors` is empty, else `400` (current behaviour —
  `resp_status = 400 if result['errors'] else 200`).

### Errors
All errors use `400` with a **machine-readable string in `error`** plus context —
the house convention across the suite (e.g. JobViewSet's
`{"error": "local_execution_unavailable", "reason": …}`). The status code stays
conventional; the *discriminator lives in the body*. The wizard branches on
`error`, not on the status code.
- `400` `{ "error": "bad_input", "reason": "Must provide sqlite_db file upload or db_path" }`
- `400` `{ "error": "structural_issues_unacknowledged", "reason": "…",`
  `"structure": {…}, "plan": [...] }` **(new)** — blocking-unresolved issues exist
  and `allow_structural_issues` is not `"true"`. The wizard branches on this
  `error` value and returns to Step 2 with the issues shown. (No structural gate
  exists today, so this body is new; the `400` status is not.)
- `400` `{ "error": "import_failed", "stats": {…}, "errors": [...] }` — the run
  completed but `errors` is non-empty (current behaviour: 400 when `errors`).
- `404`, `500` as validate.

**UI contract:** Step 3 (dry-run) and Step 4 (commit) both read `stats` + `plan` +
`errors`. On a `400` whose `error == "structural_issues_unacknowledged"`, the
wizard returns to Step 2 with the issues shown. On Step-4 success it then
`GET admin/import-status/` for confirmation counts.

> **Note — a small tightening of the current envelope.** Today the admin views
> return bare-sentence `error` strings (`{"error": "Must provide …"}`) and, on a
> run with `errors`, the raw importer summary at status 400 (no `error` key). This
> contract standardises *all* migration errors on `{error: <machine_string>,
> reason: <sentence>, …context}` to match JobViewSet. That's a **behaviour change
> to the two existing bodies** (`bad_input`, `import_failed`) — minor, but call it
> out in the PR since a CLI/script author could be parsing the old strings.

---

## 3. `GET admin/import-status/` — post-import confirmation (current, unchanged)

```jsonc
{
  "projects": { "total": 10, "groups": 1, "memberships": 3, "tags": 6 },
  "jobs":     { "total": 340, "value_keys": 40, "float_values": 5000, "char_values": 800, "xdata": 12 },
  "files":    { "total": 1198, "types": 60, "uses": 1800, "imports": 44, "exports": 3 }
}
```

Wizard shows these as a "migration landed" panel after Step 4.

---

## Delta summary — what the backend must add

1. `SQLiteValidator`: new `_validate_structure(cur)` producing `structure` +
   the per-project `plan`; extend `summary` with `structure_issues`,
   `blocking_issues`, `plan_summary`; make `summary.ok` false on blocking issues.
2. `SQLiteImporter`: consume `plan`; `copy_files`/`dest_root`; new `stats` keys;
   echo `plan` + recomputed `structure`; honour `allow_structural_issues` → `400`
   with `error: "structural_issues_unacknowledged"`.
3. `admin_views`: pass `copy_files`, `dest_root`, `allow_structural_issues`
   through; add the structural-gate branch; standardise error bodies on
   `{error: <machine_string>, reason: …}`. Success envelopes unchanged.

No existing field is renamed or removed, so the frontend can be built against this
contract in full before the backend delta lands (it just sees `structure`/`plan`
absent until then — treat missing as "no structural data yet").
