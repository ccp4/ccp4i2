# Service contract — surface audit (v1 gap analysis)

*Companion to [`CCP4I2_SERVICE_CONTRACT.md`](CCP4I2_SERVICE_CONTRACT.md). The contract document is the *promise*. This document is the *evidence*: every CCP4i2 REST endpoint hit by the two known consumer codebases — the React renderer (`client/renderer/`) and the [`i2remote` CLI](https://github.com/newcastleuniversity/materia/blob/main/Docker/cli/i2remote.py) — tabulated against v0 to expose the gap.*

## Why this exists

The v0 contract documents 9 endpoints. The renderer alone consumes ~67 distinct endpoints; `i2remote` consumes ~46 distinct endpoints (substantial overlap with the renderer + a `/uploads/*` family of its own). v0 is therefore far too thin to be useful as a contract for any operational consumer — at present it's read-only-summary-only, and even that is incomplete (no field-level promise on the file detail shape, no contracted query parameters, etc.).

The purpose of this audit is to give the CCP4i2 dev team a working artefact for deciding what gets promoted into the contract. Three buckets:

- **v0 today** — the 9 endpoints already documented as stable.
- **v1 candidates** — heavily used by both consumers, semantics mature, low risk of churn. Lifting these into the contract is mostly a documentation exercise.
- **Out of scope** — admin / internal / plumbing surfaces that should *not* be contract.

## Methodology

- **Renderer**: grepped `client/renderer/` for `apiPost`, `apiPatch`, `apiPut`, `apiDelete`, `apiFetch`, `apiJson`, `useApi().get<...>(...)`, and direct `fetch('/api/proxy/ccp4i2/...')`. Endpoint paths are relative — `createApiFetch` prefixes `/api/proxy/ccp4i2/` at runtime.
- **i2remote**: grepped `Docker/cli/i2remote.py` (lives in `newcastleuniversity/materia` post-cut) for `self.get(...)`, `self.post(...)`, `self.delete(...)`, and direct `requests.post(...)` calls.
- Endpoints are normalised to their path templates (`projectgroups/{id}/sites/`, etc.) regardless of whether the call uses an interpolated string or a literal one.

Counts deliberately err toward over-reporting — listing an endpoint twice (once per consumer) when it's the same wire is informative, not noise.

## Endpoint inventory

Legend: **R** = renderer, **i2** = i2remote, **v0** = in current v0 contract, **rec** = recommended action (`v0+` = promote to v0, `v1` = candidate for v1, `internal` = leave undocumented).

### `/projects/`

| Method | Path | R | i2 | v0 | rec | Notes |
|---|---|---|---|---|---|---|
| GET | `/projects/` | ✓ | ✓ | ✓ | — | List projects |
| GET | `/projects/{id}/` | ✓ | ✓ | ✓ | — | Detail |
| POST | `/projects/` | ✓ | ✓ | ✓ | — | Create |
| PATCH | `/projects/{id}/` | ✓ | | | **v0+** | Rename / update metadata. Materia would need this. |
| DELETE | `/projects/{id}/` | ✓ | | | **v0+** | Project deletion. Operational must-have. |
| GET | `/projects/{id}/jobs/` | ✓ | ✓ | | **v0+** | List jobs in project. Clear semantics; both consumers use it. |
| GET | `/projects/{id}/job_tree/` | | ✓ | | v1 | Hierarchical job view. CLI-shaped output; possibly renderer too. |
| GET | `/projects/{id}/files/` | ✓ | | | **v0+** | List files in project. |
| GET | `/projects/{id}/directory/` | ✓ | | | v1 | Filesystem directory tree; OS-leaky shape. Worth deliberation. |
| GET | `/projects/{id}/tags/` | ✓ | | | v1 | Project-tag listing. |
| POST | `/projects/{id}/tags/` | ✓ | | | v1 | Add tag. |
| DELETE | `/projects/{id}/tags/{tagId}/` | ✓ | | | v1 | Remove tag. |
| POST | `/projects/{id}/add_tag_by_text/` | ✓ | | | v1 | Get-or-create tag by text. Used by Campaigns batch-import. |
| POST | `/projects/{id}/create_task/` | ✓ | | | **v0+** | **Create a job in a project.** Without this, the contract is read-only. |
| POST | `/projects/{id}/export/` | ✓ | ✓ | | **v0+** | Start project export. Both consumers use it. |
| GET | `/projects/{id}/resolve_fileuse/` | | ✓ | | v1 | Resolve a `fileUse:projN/jobM/paramX` reference to a concrete file ID. |

### `/jobs/`

| Method | Path | R | i2 | v0 | rec | Notes |
|---|---|---|---|---|---|---|
| GET | `/jobs/` | ✓ | ✓ | ✓ | — | List |
| GET | `/jobs/{id}/` | ✓ | ✓ | ✓ | — | Detail |
| GET | `/active_jobs/` | ✓ | | ✓ | — | Currently running |
| PATCH | `/jobs/{id}/` | ✓ | | | **v0+** | Update job metadata (title, comments, evaluation). |
| DELETE | `/jobs/{id}/` | ✓ | | | **v0+** | Delete a job. |
| POST | `/jobs/{id}/run/` | ✓ | ✓ | | **v0+** | **Submit job to queue.** Operational must-have. |
| POST | `/jobs/{id}/run_local/` | ✓ | | | **v0+** | Synchronous (non-queued) run. Used by the Campaigns import dialogs for short-running imports that must complete before the next step. Promoted to v0+ on the "Campaigns scope" framing below. |
| POST | `/jobs/{id}/clone/` | ✓ | ✓ | | **v0+** | Clone a job for re-run with same/edited parameters. |
| POST | `/jobs/{id}/cancel/` | ✓ | | | **v0+** | Cancel a running job. |
| POST | `/jobs/{id}/preview/` | ✓ | | | v1 | Launch a viewer (coot / viewhkl / ccp4mg / terminal). Desktop-Electron-flavoured; less obvious for cloud consumers. |
| POST | `/jobs/{id}/regenerate_report/` | ✓ | | | v1 | Re-render the job's HTML report. |
| POST | `/jobs/{id}/set_parameter/` | ✓ | ✓ | | **v0+** | Set a single parameter on a job (object-path-keyed). |
| GET | `/jobs/{id}/get_parameter/` | | ✓ | | v1 | Read a single parameter. CLI uses for scripted edit-then-run. |
| PUT | `/jobs/{id}/params_xml/` | ✓ | | | v1 | Replace parameters as a wholesale XML document. Power-user surface. |
| POST | `/jobs/{id}/upload_file_param/` | ✓ | ✓ | | **v0+** | Multipart upload of a file parameter. Different shape from JSON-body endpoints. |
| POST | `/jobs/{id}/set_context_job/` | ✓ | | | v1 | Set "context" parent job for inherited parameters. |
| GET | `/jobs/{id}/container/` | ✓ | ✓ | | **v0+** | Fetch the parameter container as JSON. Heavily used. |
| GET | `/jobs/{id}/files/` | ✓ | ✓ | | **v0+** | List files associated with the job (inputs + outputs). |
| GET | `/jobs/{id}/dependent_jobs/` | ✓ | | | v1 | Jobs that depend on this job's outputs. |
| GET | `/jobs/{id}/validation/` | ✓ | | | v1 | Pre-run validation results (errors + warnings, severity-classified). |
| GET | `/jobs/{id}/run_time_validation/` | ✓ | | | v1 | Pre-flight validation including expensive checks (file content). |
| GET | `/jobs/{id}/what_next/` | ✓ | | | v1 | Suggested follow-on tasks given this job's outputs. |
| GET | `/jobs/{id}/object_method/` | ✓ | | | v1 | Generic RPC into a Python object on the server (e.g., weight calc). Powerful but smelly; review before contracting. |
| GET | `/jobs/{id}/report_xml/` | ✓ | ✓ | | v1 | Job-report XML document. Both consumers; shape is XML-flavoured. |
| GET | `/jobs/{id}/export_job/` | | ✓ | | v1 | Single-job export bundle (CLI-flavoured). |
| POST | `/jobs/bulk_dependent_jobs/` | ✓ | | | v1 | Batch dependent-jobs lookup; used by jobs-list rendering. |
| POST | `/jobs/bulk_delete/` | ✓ | | | v1 | Batch delete. |

### `/files/`

| Method | Path | R | i2 | v0 | rec | Notes |
|---|---|---|---|---|---|---|
| GET | `/files/{id}/` | ✓ | ✓ | ✓ | — | File metadata |
| POST | `/files/` | ✓ | | | **v0+** | Upload/import file to project. |
| PATCH | `/files/{id}/` | ✓ | | | **v0+** | Update file metadata (annotation, name). |
| DELETE | `/files/{id}/` | ✓ | | | **v0+** | Delete file. |
| GET | `/files/{id}/download/` | ✓ | | | **v0+** | Download binary content. |
| GET | `/files/{id}/download_by_uuid/` | ✓ | | | v1 | Same content, addressed by UUID rather than ID. |
| GET | `/files/{id}/digest/` | | ✓ | | v1 | File digest (CLI uses for checksum-based dedup). |
| GET | `/files/{id}/usage/` | ✓ | | | v1 | Where this file is used (which jobs, as which params). |
| POST | `/files/{id}/preview/` | ✓ | | | v1 | Launch viewer (desktop-flavoured). |

### `/projectgroups/` (Campaigns)

This whole resource family is missing from v0. Both consumers use it. Materia inherits the renderer's usage straight through. **All entries below are v1 candidates** unless flagged otherwise.

| Method | Path | R | i2 | rec | Notes |
|---|---|---|---|---|---|
| GET | `/projectgroups/?type=fragment_set` | ✓ | ✓ | v1 | List campaigns (filter by type) |
| GET | `/projectgroups/{id}/` | ✓ | ✓ | v1 | Campaign detail (includes nested memberships) |
| POST | `/projectgroups/` | ✓ | | v1 | Create |
| POST | `/projectgroups/create_with_parent/` | ✓ | ✓ | v1 | Create + auto-parent project |
| PATCH | `/projectgroups/{id}/` | ✓ | | v1 | Update |
| DELETE | `/projectgroups/{id}/` | ✓ | | v1 | Delete |
| GET | `/projectgroups/{id}/parent_project/` | ✓ | ✓ | v1 | The parent project |
| GET | `/projectgroups/{id}/member_projects/` | ✓ | ✓ | v1 | Members + KPIs + per-job summaries |
| GET | `/projectgroups/{id}/parent_files/` | ✓ | ✓ | v1 | Parent's coordinate + FreeR file refs |
| GET | `/projectgroups/{id}/sites/` | ✓ | | v1 | Binding sites (Moorhen-camera-flavoured) |
| PUT | `/projectgroups/{id}/sites/` | ✓ | | v1 | Update binding sites |
| POST | `/projectgroups/{id}/set_parent/` | ✓ | | v1 | Set/replace parent |
| POST | `/projectgroups/{id}/add_member/` | ✓ | ✓ | v1 | Add member project |
| DELETE | `/projectgroups/{id}/members/{pid}/` | ✓ | ✓ | v1 | Remove member |
| GET | `/projectgroups/{id}/pandda_data/` | | ✓ | v1 | PanDDA-bundle data (CLI-flavoured) |
| POST | `/projectgroups/project_campaigns/` | ✓ | | v1 | Reverse lookup: which campaigns reference these projects (POST to dodge URL-length limits on long ID lists) |

Supporting types (in [client/renderer/types/campaigns.ts](../client/renderer/types/campaigns.ts)): `ProjectGroup`, `ProjectGroupType`, `MembershipType`, `ProjectGroupMembership`, `ProjectGroupDetail`, `JobSummary`, `ProjectKPIs`, `CampaignJobInfo`, `MemberProjectWithSummary`, `ParentFilesResponse`, `CampaignSite`. Around 11 distinct types, ~150 lines of contract.

### `/uploads/` (staged uploads)

i2remote-only family for resumable / chunked / large-file uploads. Server-side likely gateway-managed. Probably belongs in v1 if external CLI consumers will exist, otherwise leave undocumented.

| Method | Path | i2 | rec | Notes |
|---|---|---|---|---|
| POST | `/uploads/request` | ✓ | v1 | Request a staged-upload session (returns upload_id + presigned URLs) |
| GET | `/uploads/` | ✓ | v1 | List staged uploads |
| GET | `/uploads/{id}/` | ✓ | v1 | Status of a staged upload |
| POST | `/uploads/{id}/complete` | ✓ | v1 | Mark upload complete; trigger server-side ingest |
| DELETE | `/uploads/{id}/cancel` | ✓ | v1 | Cancel |
| POST | `/uploads/{id}/reset` | ✓ | v1 | Reset back to pending |
| POST | `/uploads/{id}/force-complete` | ✓ | internal | Operator escape hatch; not for external consumers |

### Smaller resource families

| Method | Path | R | i2 | rec | Notes |
|---|---|---|---|---|---|
| GET | `/projecttags/` | ✓ | | v1 | List all project tags (across projects) |
| POST | `/projecttags/` | ✓ | | v1 | Create a tag |
| GET | `/projectexports/` | ✓ | | v1 | List export tasks |
| DELETE | `/projectexports/{id}/` | ✓ | | v1 | Delete completed export |
| GET | `/task_lookup/` | ✓ | | **v0+** | Task name → metadata mapping. Stable, small, used widely; trivially v0. |
| GET | `/version/` | | | ✓ | Already in v0; not directly hit by either consumer in the audit (renderer reads from a Next.js BFF route) |
| GET | `/health/` | | | ✓ | Already in v0; framework-managed (orchestrator probes) |

## Summary by recommendation

### Promote to v0 (17 endpoints)

These are the gaps that genuinely block "operational consumer". A consumer with only the current v0 surface can list and inspect; with the additions below, they can *use* the system.

- `PATCH /projects/{id}/`, `DELETE /projects/{id}/`
- `GET /projects/{id}/jobs/`, `GET /projects/{id}/files/`
- `POST /projects/{id}/create_task/`, `POST /projects/{id}/export/`
- `PATCH /jobs/{id}/`, `DELETE /jobs/{id}/`
- `POST /jobs/{id}/run/`, `POST /jobs/{id}/run_local/`, `POST /jobs/{id}/clone/`, `POST /jobs/{id}/cancel/`
- `POST /jobs/{id}/set_parameter/`, `POST /jobs/{id}/upload_file_param/`
- `GET /jobs/{id}/container/`, `GET /jobs/{id}/files/`
- `POST /files/`, `PATCH /files/{id}/`, `DELETE /files/{id}/`, `GET /files/{id}/download/`
- `GET /task_lookup/`

### v1 candidates (~50 endpoints)

The Campaigns / `projectgroups/` family (16), plus the `uploads/` family (6), plus the longer tail of job/file/project endpoints not in the v0+ list. These are real surfaces with real consumers; the question is whether they're stable enough to commit to today, or whether some need refinement first.

Worth especially close review:
- `GET /jobs/{id}/object_method/` — generic RPC into a server-side Python object. Powerful, currently used by validation/weight-calc flows. Architecturally smelly (effectively an internal RPC tunnelled through the API). Probably should *not* go into a contract; consider replacing with named endpoints.
- `PUT /jobs/{id}/params_xml/` — wholesale XML replace. Power-user; possibly fine, but the XML format itself becomes contract.
- `GET /projects/{id}/directory/` — leaks filesystem layout. Cloud deployments hide this; desktop exposes it. Worth thinking about whether the contract should commit to it at all.

### Internal / out of scope

- `POST /uploads/{id}/force-complete` — operator escape hatch.
- Any endpoint we don't list (admin views, debug, plumbing).
- Per-task interface JSON / phil / def.xml shapes (already explicitly out of scope in the v0 doc).

## Discussion prompts for the dev team

1. **Where's the cutoff between v0 and v1?** This audit suggests v0+v1 ≈ 67 endpoints (everything except internal). That's a much bigger commitment than v0's current 9. Is the right move to promote the operational subset (~25) into v0 and treat the long tail as v1 candidates with a separate cadence? Or keep v0 narrow and commit to a single richer v1?

2. **`object_method` and `params_xml`** — these are tunnels for internal mechanisms. Promoting them to contract locks in current implementation choices. Worth deliberating whether they should be re-shaped before contracting, or contracted as-is and re-shaped behind a v2 transition.

3. **Campaigns / `projectgroups/`** — v0's framing was "projects + jobs + files". Campaigns is a *fourth* resource family of comparable substance. Either:
   - It belongs in v0 (because it's load-bearing for the only known external consumer, Materia), or
   - It belongs in v1 with explicit "consumer must opt-in to depend on this" language.

   **Position recorded** (Martin Noble, 2026-04-30): if pushed to a position, promote to v0 *whatever the Campaigns code path actions*. Scope: every endpoint hit by `client/renderer/lib/campaigns-api.ts`, `client/renderer/components/campaigns/*`, and `client/renderer/components/moorhen/campaign-*.tsx`. That captures: all `/projectgroups/*` endpoints in the table above, plus `POST /jobs/{id}/run_local/` (used by the import dialogs). Already in v0+ via this framing, no further additions needed.

4. **Bulk endpoints** (`jobs/bulk_dependent_jobs`, `jobs/bulk_delete`, `projectgroups/project_campaigns`). These exist for a reason (URL-length limits, transactional batch ops). Are they part of the contract, or implementation detail of the renderer that a different consumer might re-derive differently?

5. **Pagination** — flagged in v0 already. Particularly relevant for `/projects/`, `/jobs/`, `/projectgroups/` once datasets grow. The shape of the breaking change matters.

6. **Filter / query params** — v0 mentions `?project={id}` for `/jobs/`. The audit shows `?type=fragment_set` for `/projectgroups/`, plus an implicit `django-filter` surface (`?tags__text=foo` and similar) that consumers may have come to depend on without a written promise. Worth making explicit which query params are contracted vs available-but-volatile.

## How to use this document

- For the contract review meeting: walk the v0+ list first (small, uncontroversial). Then debate the v1 candidates one resource family at a time. Then decide what to do about the smelly bits.
- After ratification: add the agreed v0+ rows to [`CCP4I2_SERVICE_CONTRACT.md`](CCP4I2_SERVICE_CONTRACT.md) "Stable endpoints", commit a corresponding npm minor (0.2.0), and add the relevant Django shape guards to [`server/ccp4i2/tests/api/unit/test_contract.py`](../server/ccp4i2/tests/api/unit/test_contract.py).
- This audit document gets re-generated when a new significant consumer comes online or when a meaningful chunk of new endpoints lands in the codebase. Otherwise it's a one-shot snapshot.
