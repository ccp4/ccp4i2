# CCP4i2 service contract — v0

*The promise: "what shape of REST API can third parties (Materia, future integrators) build against and trust will not silently change?"*

This document is the externally-facing answer to that question. The contract is **draft v0** — small surface, well-defined shapes, evolves deliberately. Endpoints not listed here are CCP4i2-internal and may change without notice.

## Status

- **Draft v0**, circulating with the CCP4i2 dev team for review.
- TypeScript types live in [`packages/ccp4i2-auth/src/contracts/ccp4i2.ts`](../packages/ccp4i2-auth/src/contracts/ccp4i2.ts) and re-export from [`@ccp4/ccp4i2-auth`](https://www.npmjs.com/package/@ccp4/ccp4i2-auth) on npm. External consumers `import type { Project, Job, File, ... } from "@ccp4/ccp4i2-auth"`.
- The first known external consumer is [Materia](https://github.com/newcastleuniversity/materia) (Newcastle's SBDD bench, which embeds CCP4i2 as a Python+JS component). The strategic context behind the contract — why split, why a stable surface, what CCP4i2 commits to — lives in Materia's [`apps/compounds/docs/`](https://github.com/newcastleuniversity/materia/tree/main/apps/compounds/docs) (specifically `CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md` and `CCP4I2_RELATIONSHIP_AND_SUSTAINABILITY.md`); read those if you want the *why*. This document is the *what*.

## Conventions

| Concern | Decision |
|---|---|
| Base URL | Relative to deployment, typically `/api/ccp4i2/`. The proxy / ingress prefix (`/api/proxy/ccp4i2/`) is a deployment-specific concern. |
| Authentication | `Authorization: Bearer <token>`. Token providers (MSAL for cloud, LocalSession for desktop) live in `@ccp4/ccp4i2-auth`. |
| 401/403 response shape | `{success: false, error: string}`. Pattern-matched by `AUTH_ERROR_EVENT` listeners. Stable. |
| Date format | ISO 8601 strings, UTC (e.g., `"2026-04-29T13:21:24.854154Z"`). Stable. |
| ID types | `id` is integer (Django PK); `uuid` is a UUID v4 string. Both are stable identifiers; `uuid` is preferred for cross-system references. |
| Unknown fields | Consumers MUST ignore fields they don't recognise. CCP4i2 reserves the right to add new fields without bumping the contract version. |

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

Lists all jobs across projects the caller has access to. Filter by `project={id}` query param.

### `GET /jobs/{id}/` — job detail

**Response:** `Job`

### `GET /active_jobs/` — currently-running jobs

**Response:** `ActiveJobsResponse`

Lightweight summary suitable for per-second polling. Lists jobs whose `status` is in the running set (Pending, Queued, Running, RunningRemotely).

### `GET /files/{id}/` — file detail

**Response:** `File`

Includes a computed `path` field (full filesystem path, derivable from `directory` + `name` + `job`).

## Type stability

The TypeScript types in [`packages/ccp4i2-auth/src/contracts/ccp4i2.ts`](../packages/ccp4i2-auth/src/contracts/ccp4i2.ts) are the contract. Specifically:

- **Documented fields are stable.** They will not be removed, renamed, or have their type changed without a major version bump of `@ccp4/ccp4i2-auth`.
- **Numeric enum values are stable.** `JobStatus.Running === 3` is part of the contract — consumers may pattern-match on values directly.
- **Undocumented fields are unstable.** They may exist in responses today (Django serializers using `fields = "__all__"`) but are not part of the contract; CCP4i2 may remove them at any time.
- **The auth-error response shape** (`{success: false, error: string}`) is stable across all endpoints.

## Out of scope (subject to change without notice)

- Plugin / wrapper-specific endpoints (anything under `/jobs/{id}/...` action paths).
- Admin / legacy-import endpoints (`/admin/import-*`). Cloud-only and operator-facing.
- Per-task interface JSON shapes (each task has its own free-form parameter schema).
- Internal job-value endpoints (`/job_float_values/`, `/job_char_values/`).
- File serving by path (`/projects/{id}/files_by_path/...`) — convenience surface that may be refactored.
- The async execution / process-manager surface.

External consumers should NOT depend on these surfaces. If a consumer needs functionality not in v0, raise an issue against this document so the surface can be considered for inclusion in v1.

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

if (myJob.status === JobStatus.Running) { /* ... */ }
```

## Contract guards (Django side)

[`server/ccp4i2/tests/api/unit/test_contract.py`](../server/ccp4i2/tests/api/unit/test_contract.py) asserts the shape of contract endpoints — if a future refactor changes a documented field, the test fails. A v0 starter set of guards is included; expand as needed when surfaces are added to the contract above.

## Versioning policy

- **Patch / minor changes** to `@ccp4/ccp4i2-auth` may add new types, add optional fields, or relax types. Consumers don't need to update.
- **Major version bumps** are required for: removing a field, renaming a field, changing a field's type in an incompatible way, changing a numeric enum value, changing the 401/403 response shape.
- The CCP4i2 deployment carries a runtime `version` (visible at `/version/`). Consumers are encouraged to log or display it for debugging — but the type contract is what they build against, not the runtime version.

## Open questions

1. **Pagination.** Do `/projects/`, `/jobs/`, etc. paginate when result sets grow large? Currently they don't; future versions may add cursor pagination. If/when that lands, consumers will get an envelope shape (`{results: T[], next: string | null}`) — that's a breaking change worth flagging in advance.
2. **Filter / search query params.** What query parameters are stable? Currently `?project={id}` is documented for `/jobs/`; others (e.g., `?tags__text=foo` from django-filter) are de facto available but not contracted.
3. **WebSocket / SSE for live job-status updates.** Currently consumers poll `/active_jobs/`. A push channel would benefit Materia's compounds-side dashboards but isn't in v0.
4. **Rate-limiting headers.** None documented; should they be?

These should be settled with the CCP4 dev team during contract review.
