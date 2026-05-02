# CCP4i2 service contract — v0

_The promise: "what shape of REST API can third parties (Materia, future integrators) build against and trust will not silently change?"_

This document is the externally-facing answer to that question. The contract is **draft v0** — small surface, well-defined shapes, evolves deliberately. Endpoints not listed here are CCP4i2-internal and may change without notice.

## Status

- **Draft v0**, circulating with the CCP4i2 dev team for review.
- TypeScript types live in [`packages/ccp4i2-auth/src/contracts/ccp4i2.ts`](../packages/ccp4i2-auth/src/contracts/ccp4i2.ts) and re-export from [`@ccp4/ccp4i2-auth`](https://www.npmjs.com/package/@ccp4/ccp4i2-auth) on npm. External consumers `import type { Project, Job, File, ... } from "@ccp4/ccp4i2-auth"`.
- The first known external consumer is [Materia](https://github.com/newcastleuniversity/materia) (Newcastle's SBDD bench, which embeds CCP4i2 as a Python+JS component). The strategic context behind the contract — why split, why a stable surface, what CCP4i2 commits to — lives in Materia's [`apps/compounds/docs/`](https://github.com/newcastleuniversity/materia/tree/main/apps/compounds/docs) (specifically `CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md` and `CCP4I2_RELATIONSHIP_AND_SUSTAINABILITY.md`); read those if you want the _why_. This document is the _what_.

## Conventions

| Concern | Decision |
| --- | --- |
| Base URL | Relative to deployment, typically `/api/ccp4i2/`. The proxy / ingress prefix (`/api/proxy/ccp4i2/`) is a deployment-specific concern. |
| Authentication | `Authorization: Bearer <token>`. The bearer-token shape is the contract; how the consumer acquires it is consumer-specific. `@ccp4/ccp4i2-auth` includes provider abstractions for the two browser/Electron contexts (`MsalBearerTokenProvider` for cloud, `LocalSessionTokenProvider` for desktop); CLI clients like `i2remote` supply tokens directly (e.g. via `az account get-access-token` or a service-principal credential). The server validates per the configured middleware (LocalSession / AzureAD). |
| 401/403 response shape | `{success: false, error: string}`. Pattern-matched by `AUTH_ERROR_EVENT` listeners. Stable. |
| Date format | ISO 8601 strings, UTC (e.g., `"2026-04-29T13:21:24.854154Z"`). Stable. |
| ID types | `id` is integer (Django PK), unique within a single deployment. `uuid` is a UUID v4 string, globally unique. **Most endpoints address resources by `{id}`** (e.g. `/projects/{id}/`, `/jobs/{id}/`); a small **uuid-addressed surface exists on files** (`/files_by_uuid/{uuid}/`, `/download/`, `/digest/`) for cross-deployment references — the parameter-file format only knows uuids, not deployment-local ids. Both fields are stable. Use `id` for normal in-deployment work; use `uuid` when the reference must survive a deployment boundary or be embedded in cross-system data. |
| Unknown fields | Consumers MUST ignore fields they don't recognise. CCP4i2 reserves the right to add new fields without bumping the contract version. |
| Nested children | Detail responses do not inline **unbounded** child collections (jobs and files for a project; files for a job). Fetch those via the dedicated endpoints. Bounded / inherently-small relations (`tags` on a project, `memberships` on a project-group) **may** be inlined; the TS type definition is authoritative for which is which. |
| Live updates | v0 is **poll-based**. `/active_jobs/` is designed for per-second polling; per-job SWR refresh intervals are the per-job pattern. Push channels (WebSocket / SSE) are not currently contracted; if added in a future version they will be **additive** — polling consumers will not need to change. |
| Rate limiting | No rate-limit response headers are contracted. Deployments may impose rate limits at the ingress / gateway layer (Azure Front Door, nginx, etc.); the response for a rate-limited request is **deployment-specific** (typically a 429 with a `Retry-After` header from the gateway, but this is *not* a CCP4i2 contract). |

## Stable endpoints

The endpoints below are the v0 surface. Each entry lists method, path, request shape (where relevant), response type, and stability notes.

### `GET /version/` — server version info

Public; no auth required.

**Response:** `VersionInfo`

Used by consumers to render compatibility banners and to detect when an instance has been upgraded.

### `GET /health/` — server health

Public; no auth required.

**Response:** `HealthStatus`

Used by container-orchestration health probes and by clients' connectivity diagnostics.

### `GET /projects/` — list projects

**Response:** `ProjectListItem[]`

Lightweight list view; excludes `directory` and inlines a summary of tags. Use `/projects/{id}/` for full detail.

### `GET /projects/{id}/` — project detail

**Response:** `Project`

### `POST /projects/` — create project

**Request body:** `Partial<Project>` — minimum `{name}` required. Server fills `directory` if not supplied.

**Response:** `Project` (201).

### `GET /jobs/` — list jobs

**Response:** `Job[]`

Lists all jobs across projects the caller has access to.

**Contracted query parameters:**
- `?project={id}` — restrict to jobs of a specific project.

### `GET /jobs/{id}/` — job detail

**Response:** `Job`

The `Job` shape inlines two KPI dictionaries — `float_values` and `char_values` — keyed by KPI name (e.g. `{RfreeFinal: 0.234, ...}` and `{SpaceGroup: "P 21 21 21", ...}`). The set of names is task-specific (refmac emits `RfreeFinal`, aimless emits `ResolutionLow`, etc.); consumers should treat both dicts as open-ended and look up names they care about, not iterate exhaustively. Empty objects when the job hasn't produced KPIs yet (pending / running jobs, or tasks that don't emit any). Promoted to v0 in 0.3.0; consumers building against earlier versions of this contract should default both fields to `{}` for robustness.

### `GET /active_jobs/` — currently-running jobs

**Response:** `ActiveJobsResponse`

Lightweight summary suitable for per-second polling. Lists jobs whose `status` is in the running set (Pending, Queued, Running, RunningRemotely).

### `GET /files/{id}/` — file detail

**Response:** `File`

Includes a computed `path` field (full filesystem path, derivable from `directory` + `name` + `job`).

### `GET /files/{id}/download/` — file content

Returns the raw file content as a stream. The `Content-Type` is derived from the file's type registration; the `Content-Disposition` header carries the original filename.

### `GET /files/?job={id}` — list files of a job

**Response:** `File[]`

Lists every file (input + output) bound to a specific job.

**Contracted query parameters:**
- `?job={id}` — restrict to files of a specific job. (Required for this use case; the unfiltered `GET /files/` returns the global file list.)

### `GET /files/?project={id}` — list files of a project

**Response:** `File[]`

Lists every file across all jobs of a specific project. Equivalent to fetching all of the project's jobs and unioning their files.

**Contracted query parameters:**
- `?project={id}` — restrict to files of a specific project.

### `GET /files_by_uuid/{uuid}/` — file metadata by UUID

**Response:** `File`

Resolves a file by its UUID. CCP4i2 parameter files (`def.xml`, fileUse refs, exported job bundles) carry file references as UUIDs — they have no knowledge of deployment-local Django ids. This endpoint is the canonical resolution path for those refs.

### `GET /files_by_uuid/{uuid}/download/` — file content by UUID

Returns the raw file content as a stream. Same shape as `GET /files/{id}/download/`, addressed by UUID.

### `GET /files_by_uuid/{uuid}/digest/` — file digest by UUID

**Response:** task-specific summary JSON.

Returns a JSON digest of the file's contents — sequences for FASTA, cell + spacegroup for MTZ, ligand summary for refmac dictionary, etc. The digest *shape* is content-type-specific (each file type has its own `CDataFileContent` subclass on the server); the *availability* of the endpoint is contracted, the *contents* are not. Used by the renderer to populate parameter-summary side-panels without re-downloading the full file.

### `POST /projects/{id}/create_task/` — create a new job

**Request body:**
```json
{ "task_name": "<task-name>" }
```

Creates a new job (an instance of the named task) in the specified project. The new job has a fresh container with default parameter values; consumers typically follow up with one or more `POST /jobs/{id}/set_parameter/` calls to populate it before `POST /jobs/{id}/run/`.

**Response:**
```json
{ "success": true, "data": { "new_job": Job } }
```

The valid `task_name` values are CCP4i2-internal — see [`core/task_manager/plugin_lookup.json`](../core/task_manager/plugin_lookup.json) for the registry of installed task plugins. Names are stable but the registry list is *not* part of the contract.

### `POST /projects/{id}/add_tag_by_text/` — get-or-create a tag and apply

**Request body:**
```json
{ "text": "<tag-text>" }
```

Looks up a tag by its text; creates one if it doesn't exist; applies it to the project. Used by the renderer's batch-import flow to tag projects with site-of-origin labels without the consumer needing to know whether the tag pre-existed.

**Response:**
```json
{
  "status": "<created|already-applied|...>",
  "tag": { "id": <id>, "text": "<text>" },
  "created": <bool>,
  "message": "<human-readable>"
}
```

### `POST /jobs/{id}/run/` — submit a job to its queue

**Request body:** `{}` (or `{ "queue": "<queue-name>" }` to override the default).

Submits the job for asynchronous execution by the worker pool. Status transitions follow `Pending → Queued → Running → Finished | Failed | Interrupted`; consumers poll `GET /jobs/{id}/` or `GET /active_jobs/` for status updates.

**Response:** `{ "success": true }` on submission accepted (does NOT wait for completion).

### `POST /jobs/{id}/run_local/` — run a job synchronously

**Request body:** `{}` (or `{ "synchronous": true }`).

Runs the job in-process and returns when it completes. Used for short-running tasks (file imports, format conversions) where the caller needs to wait for the result before proceeding. Long-running tasks should always use `POST /jobs/{id}/run/` instead.

**Response:** `{ "success": true, "data": { ... } }` on completion.

### `POST /jobs/{id}/clone/` — clone a job for re-run

**Request body:** `{}`.

Creates a new job that is a parameter-by-parameter copy of this one. Useful for re-running with the same configuration, optionally after `set_parameter` tweaks. The cloned job inherits the parent project, task name, and parameter container; it does not inherit `status`, `process_id`, or output files.

**Response:**
```json
{ "success": true, "data": { "new_job": { "id": <id>, "uuid": "<uuid>" } } }
```

### `POST /jobs/{id}/set_parameter/` — set a single parameter

**Request body:**
```json
{ "object_path": "container.inputData.XYZIN", "value": <task-specific> }
```

Sets a single parameter on the job's container by dotted object-path. Endpoint and envelope shape are contracted; **`object_path` strings and `value` types are task-specific** and not contracted.

- `object_path` syntax is `container.<section>.<param>[.<sub>...]` where `<section>` is one of `inputData` / `outputData` / `controlParameters`.
- `value` is JSON whose shape depends on the parameter type: scalars for `CString`/`CFloat`/`CBoolean`; objects for file parameters (typically `{dbFileId: <id>}` referring to an existing `File` row); arrays for `CList` parameters.

Valid `object_path` values for a given job are discoverable via `GET /jobs/{id}/container/`; per-task parameter type semantics live in `tasks/<taskname>/<taskname>.def.xml` (CCP4i2-internal, not contract).

**Response:** `{ "success": true, "data": { "updated_item": ... } }` on success; error response on validation failure.

### `POST /jobs/{id}/upload_file_param/` — upload content as a file parameter

**Request body:** `multipart/form-data` with:
- `object_path`: dotted object-path string (e.g. `container.inputData.XYZIN`) — same syntax as `set_parameter`'s `object_path`.
- `file`: the file content.

Uploads new file content into a file-typed parameter slot on the job's container. The server creates a `File` record + `FileImport` record, stores the content in the project's directory, and binds the parameter to the new file. If the parameter already had a file bound, the prior file's content and `FileImport` record are deleted.

Use this endpoint when the file content is *new* (uploaded from the client). To bind a parameter to an *existing* file (output from a prior job, file already in the project), use `POST /jobs/{id}/set_parameter/` with `value: {dbFileId: <id>}`.

**Response:** `{ "success": true, "data": { "updated_item": ... } }`.

(Legacy note: older clients send the form field as `objectPath` — camelCase. The server accepts both; new clients should use snake_case `object_path` for consistency with `set_parameter`.)

### `GET /jobs/{id}/report_xml/` — job report as XML

**Response Content-Type:** `application/json`. **Response body:** `{ "success": true, "xml": "<report>...</report>" }` — the XML report is delivered as a string inside the standard JSON envelope (not as a raw `text/xml` payload). Consumers parse it client-side (DOMParser in the renderer; an XML lib in CLI consumers) into a structured tree.

The report is the canonical task-output document — what the renderer's report panel renders, what `i2run`'s `--report` flag writes to disk. The XML *schema* is task-specific (every CCP4i2 task has its own report structure) and is **not part of this contract** — only the endpoint and envelope are. The server caches the rendered report to `report_xml.xml` in the job directory and re-uses it across requests; pass `?regenerate=true` to force regeneration.

Returns `success: false` (404) if the job has not produced a report (typical for pending / failed jobs).

### `GET /projects/{id}/resolve_fileuse/?fileuse=<expr>` — resolve a fileUse DSL string

**Response:** `ResolveFileUseResponse` on success; `{ "success": false, "error": "<reason>" }` on failure (no matching job, param name not found, index out of range, etc.).

The `fileuse` query parameter is a string in CCP4i2's fileUse DSL, used to reference a previous job's output without knowing the deployment-local file id ahead of time. Four supported syntactic forms:

```
task_name[jobIndex].jobParamName[paramIndex]      # full
[jobIndex].jobParamName[paramIndex]               # no task name
task_name[jobIndex].jobParamName                  # no param index
[jobIndex].jobParamName                           # minimal
```

`jobIndex` is positive (count from the start of the project's job list) or negative (count back from the end). `paramIndex` defaults to 0. If `task_name` is given, only jobs of that task are considered. Examples:

- `[-1].XYZOUT[0]` — first XYZOUT from the most recent job in the project
- `prosmart_refmac[-1].XYZOUT` — XYZOUT from the most recent prosmart_refmac job
- `refmac[-2].HKLOUT[0]` — HKLOUT from the second-to-last refmac job

The typical caller flow is `resolve_fileuse → set_parameter({dbFileId})` — this endpoint is the wiring primitive that pairs with `POST /jobs/{id}/set_parameter/` whenever the input being bound is the output of a previous job (rather than freshly uploaded content).

Note: the response field names (`dbFileId`, `baseName`, etc.) intentionally mirror what `set_parameter` consumes; they are **not** the same as `File`'s field names (`id`, `name`). Treat the response as a resolution result, not a `File`.

### `GET /jobs/{id}/container/` — job parameter container

**Response:** `unknown` — a task-specific JSON tree.

Returns the full CCP4i2 parameter container for the job, encoded as JSON via `CCP4i2JsonEncoder`. The endpoint is contracted; the **shape of the response body is not** — every task has its own parameter schema (RefMac differs from AimlessPipe etc.) and pipelines are trees-of-trees.

Consumers wanting to read or mutate parameters generically should treat container JSON as opaque and use `POST /jobs/{id}/set_parameter/` (object-path-keyed) which abstracts the tree shape.

Consumers wanting task-specific structured access should consult the per-task `.def.xml` parameter schemas in `tasks/<taskname>/` — those are CCP4i2-internal documentation, not contract.

## Project groups (campaigns)

Project groups bundle multiple projects under a shared parent — primarily used by the fragment-screening workflow ("campaigns") where a parent project carries reference coordinates and FreeR flags, and each member project represents a single dataset soaked with one compound. The endpoints below cover the surface needed to manage that workflow end-to-end.

Types referenced (currently in `client/renderer/types/campaigns.ts`; to be lifted into `@ccp4/ccp4i2-auth/src/contracts/ccp4i2.ts` in a follow-up version): `ProjectGroup`, `ProjectGroupType` (`"general_set" | "fragment_set"`), `ProjectGroupDetail`, `ProjectGroupMembership`, `MembershipType` (`"parent" | "member"`), `MemberProjectWithSummary`, `JobSummary`, `ProjectKPIs`, `ParentFilesResponse`, `CampaignSite`.

### `GET /projectgroups/` — list project groups

**Response:** `ProjectGroup[]`

**Contracted query parameters:**
- `?type=fragment_set` (or `?type=general_set`) — restrict to a specific kind.

### `GET /projectgroups/{id}/` — project group detail

**Response:** `ProjectGroupDetail` — `ProjectGroup` augmented with a nested `memberships: ProjectGroupMembership[]` array. Memberships are bounded per-group (the inline-nesting convention); see the *Nested children* row in *Conventions*.

### `POST /projectgroups/` — create a project group

**Request body:** `Partial<ProjectGroup>` — minimum `{ "name": "<name>", "type": "fragment_set" | "general_set" }`.

**Response:** `ProjectGroup` (201).

### `POST /projectgroups/create_with_parent/` — create a group with auto-created parent project

**Request body:** `{ "name": "<name>", "type": "fragment_set" | "general_set" }`.

Creates the project group AND a fresh parent project of the same name in one call. Saves the consumer a separate `POST /projects/` round-trip when starting a new campaign from scratch.

**Response:** `ProjectGroup` (201).

### `PATCH /projectgroups/{id}/` — update a project group

**Request body:** `Partial<ProjectGroup>`.

**Response:** `ProjectGroup`.

### `DELETE /projectgroups/{id}/` — delete a project group

Deletes the group and all its memberships. Member projects are NOT deleted (they keep existing as standalone projects); only the grouping relation is removed.

**Response:** `{ "success": true }` (204 conventional).

### `GET /projectgroups/{id}/parent_project/` — the parent project

**Response:** `Project | null`. Null when the group has no designated parent (general groups; fragment-screening campaigns always have one).

### `GET /projectgroups/{id}/member_projects/` — member projects with KPIs

**Response:** `MemberProjectWithSummary[]` — each member project enriched with `JobSummary[]` and `ProjectKPIs` (counts of jobs by status, latest-success metrics, etc.).

This endpoint is *enrichment-bearing* (not just a filtered list of `Project`); the renderer's campaign dashboard depends on the joined shape.

### `GET /projectgroups/{id}/parent_files/` — parent project's coordinate / FreeR files

**Response:** `ParentFilesResponse` — references to the parent project's coordinate file(s) and FreeR flag file(s), the inputs that get propagated to each member's refinement.

### `GET /projectgroups/{id}/sites/` — binding sites

**Response:** `CampaignSite[]` — list of binding-site landmarks (name + Cartesian origin + optional camera quaternion / zoom). Used by the Moorhen viewer to navigate between sites of interest within a campaign.

### `PUT /projectgroups/{id}/sites/` — replace binding sites

**Request body:** `CampaignSite[]` — wholesale replacement of the sites list.

**Response:** `CampaignSite[]` (the persisted list).

### `POST /projectgroups/{id}/set_parent/` — set or replace the parent project

**Request body:** `{ "project_id": <id> }`.

**Response:** `ProjectGroupMembership`.

### `POST /projectgroups/{id}/add_member/` — add a member project

**Request body:** `{ "project_id": <id>, "type": "member" | "parent" }`.

**Response:** `ProjectGroupMembership`.

### `DELETE /projectgroups/{id}/members/{projectId}/` — remove a member project

Removes the membership relation; the project itself is not deleted.

**Response:** `{ "success": true }` (204 conventional).

## Type stability

The TypeScript types in [`packages/ccp4i2-auth/src/contracts/ccp4i2.ts`](../packages/ccp4i2-auth/src/contracts/ccp4i2.ts) are the contract. Specifically:

- **Documented fields are stable.** They will not be removed, renamed, or have their type changed without a major version bump of `@ccp4/ccp4i2-auth`.
- **Numeric enum values are stable.** `JobStatus.Running === 3` is part of the contract — consumers may pattern-match on values directly.
- **Undocumented fields are unstable.** They may exist in responses today (Django serializers using `fields = "__all__"`) but are not part of the contract; CCP4i2 may remove them at any time.
- **The auth-error response shape** (`{success: false, error: string}`) is stable across all endpoints.

## Out of scope (subject to change without notice)

- Plugin / wrapper-specific endpoints (anything under `/jobs/{id}/...` action paths).
- Admin / legacy-import endpoints (`/admin/import-*`). Cloud-only and operator-facing.
- Per-task interface JSON/phil/.def.xml shapes (each task has its own free-form parameter schema).
- Internal job-value endpoints (`/job_float_values/`, `/job_char_values/`).
- File serving by path (`/projects/{id}/files_by_path/...`) — convenience surface that may be refactored.
- The async execution / process-manager surface.

External consumers should NOT depend on these surfaces. If a consumer needs functionality not in v0, raise an issue against this document so the surface can be considered for inclusion in v1.

## Deliberately omitted

These surfaces *exist* server-side and are reachable, but are intentionally not promoted to the contract. The reasoning is recorded here so the question doesn't get re-litigated.

### `GET /jobs/{id}/get_parameter/?object_path=<path>`

**Status:** server endpoint exists; not contracted.

**Reason:** for read-side parameter access, use `GET /jobs/{id}/container/` (already contracted). The per-call cost of `get_parameter` is dominated by container construction (loading the plugin, instantiating its `CContainer`, overlaying `input_params.xml` and walking the object-path) — not by JSON encoding or wire size. Returning the whole tree once is, on the round-trip economics, the same cost as returning a single field, so a separate read endpoint doesn't earn its keep. `set_parameter` *is* contracted because writes inherently address a single parameter; reads do not.

Consumers wanting structured per-task parameter access should consult the task's `def.xml` and read against the container response. If a profiling case emerges where this analysis is wrong (a hot path that needs a single field many times per second and can't cache the container), reopen.

## How a consumer adopts the contract

```typescript
import type {
  Project,
  ProjectListItem,
  Job,
  JobStatus,
  File,
  VersionInfo,
} from "@ccp4/ccp4i2-auth";
import { createApiFetch } from "@ccp4/ccp4i2-auth";

// The factory binds auth + base URL; the bound apiJson<T> typesafely
// returns the contracted shape.
const api = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });

const version = await api.apiJson<VersionInfo>("version/");
const projects = await api.apiJson<ProjectListItem[]>("projects/");
const myJob = await api.apiJson<Job>(`jobs/${jobId}/`);

if (myJob.status === JobStatus.Running) {
  /* ... */
}
```

## Contract guards (Django side)

[`server/ccp4i2/tests/api/unit/test_contract.py`](../server/ccp4i2/tests/api/unit/test_contract.py) asserts the shape of contract endpoints — if a future refactor changes a documented field, the test fails. A v0 starter set of guards is included; expand as needed when surfaces are added to the contract above.

## Versioning policy

- **Patch / minor changes** to `@ccp4/ccp4i2-auth` may add new types, add optional fields, or relax types. Consumers don't need to update.
- **Major version bumps** are required for: removing a field, renaming a field, changing a field's type in an incompatible way, changing a numeric enum value, changing the 401/403 response shape.
- The CCP4i2 deployment carries a runtime `version` (visible at `/version/`). Consumers are encouraged to log or display it for debugging — but the type contract is what they build against, not the runtime version.

## Open questions

1. **Pagination.** Do `/projects/`, `/jobs/`, etc. paginate when result sets grow large? Currently they don't; future versions may add cursor pagination. If/when that lands, consumers will get an envelope shape (`{results: T[], next: string | null}`) — that's a breaking change worth flagging in advance.
2. **django-filter de facto query parameters.** Beyond the contracted query params named per-endpoint above, many list endpoints accept django-filter-style parameters (e.g. `?tags__text=foo`, `?status=3`) by virtue of how the DRF `FilterSet` is wired. These are **not contracted** — consumers may use them experimentally but should expect them to change without notice. Future contract updates may promote specific params to stable status as their semantics get scrutinised.
3. **Surface beyond v0.** v0 was selected by the **"cousin principle"**: every endpoint actioned by the Campaigns code path (the closest cousin to a future external integrator) is in scope; everything else is deferred. That captures ~37 endpoints — projects + jobs core + Campaigns family + file mutations / uuid-addressed reads. The audit at [`CCP4I2_SERVICE_CONTRACT_AUDIT.md`](CCP4I2_SERVICE_CONTRACT_AUDIT.md) tabulates the remaining ~30 endpoints (additional job-lifecycle: cancel/PATCH/DELETE/preview/regenerate_report/etc.; project-tags CRUD; project-export/import; staged uploads; bulk operations; the smelly `object_method` and `params_xml` corners). The CCP4 dev team should decide which of those to promote in a future contract revision; consumers building today against v0 can do everything Campaigns does.

These should be settled with the CCP4 dev team during contract review.
