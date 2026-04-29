# Repo split — implementation status & working plan

*Living document. Tracks execution of the [`CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md`](CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md) plan against the [`CCP4I2_RELATIONSHIP_AND_SUSTAINABILITY.md`](CCP4I2_RELATIONSHIP_AND_SUSTAINABILITY.md) framing. Update as work lands; treat it as the single page a fresh contributor (or a fresh Claude session) reads to know where we are.*

## Status snapshot

| Component | State |
|---|---|
| **Strategic sign-off** (CCP4i2 dev team) | obtained |
| **Materia branding decision** | locked (April 2026) |
| **Git geometry on `ccp4/ccp4i2`** | `django-sliced` branch live, transitional |
| **Deploy safety (DDU protection)** | shipped — `--env` mandatory, no `:latest`, image-lineage separation |
| **Bilingual auth library** | scaffolded + active migration in progress |
| **Lightweight desktop auth** | shipped — LocalSession end-to-end, validated on macOS Electron |
| **`materia/*` cloud lineage** | **shipped + validated** — `materia-demo-web` live and signed-in-against, AzureAD chain working through refactored middleware |
| **Compounds repo cut** | not started — gated on institutional asks |
| **Service contract for CCP4i2 REST API** | not started |
| **Standalone Materia Docker image** | not started |

## Source documents (read these first)

- [**CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md**](CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md) — the externally-facing proposal that obtained the green light. The 6-step sequencing in the "Proposed sequencing" section is the parent plan that this document executes against.
- [**CCP4I2_RELATIONSHIP_AND_SUSTAINABILITY.md**](CCP4I2_RELATIONSHIP_AND_SUSTAINABILITY.md) — the internal sustainability rationale. Recommends "option C" (two repos + shared lib) and identifies bilateral-stewardship as the load-bearing risk. Consult before any structural decision.
- This document — the implementation log.

## Locked decisions

These were settled in working sessions and shouldn't be re-litigated without explicit reason.

| Topic | Decision |
|---|---|
| Platform name (compounds-side) | **Materia** |
| Materia tagline (for README first lines, search-engine discoverability) | *"An SBDD bench for compound, assay, and construct tracking"* — proper noun + descriptive subtitle pattern. Keeps Materia's brand-ability and scope-honesty while landing the SBDD keyword anywhere a peer might Google. |
| Compounds-side repo name | `materia` |
| Shared auth library repo | `ccp4/ccp4i2-auth` (CCP4-stewarded) |
| npm scope | `@ccp4/ccp4i2-auth` |
| PyPI package name | `ccp4i2-auth` |
| CCP4i2-side branch (transitional) | `django-sliced` (not "materia") |
| Django app labels | stay flat (`compounds`, `assays`, etc.) — **never** rename to `materia_*` (would force table renames and break DDU data) |
| Azure resource naming | `materia-<instance>` for new resources only; existing `ccp4i2-bicep-*` stays |
| Image lineage in shared ACR | `ccp4i2/{web,server}` for stable; `materia/{web,server}` for development |
| Deploy script discipline | `--env <env-file>` mandatory, no implicit DDU default |
| `:latest` tag policy | never pushed; consumers reference explicit timestamp tags |
| Bilingual auth library layout | one repo, two language artifacts (TS via npm + Python via PyPI), one canonical contract |
| Auth-middleware selection (server) | fail-closed three-way switch in `settings.py` (LocalSession / AzureAD / DevAdmin / nothing) |
| `DevAdminMiddleware` activation | `DEBUG=True` gated, defence-in-depth |
| Forward-port direction | `django` → `django-sliced` regularly; rarely the other way |
| `apps/compounds/` rename | **deferred** to the moment of `git filter-repo` extraction (no in-flight rename) |
| `apps/users/` ownership | **Materia-resident.** No `ccp4i2-users` shared library. CCP4i2-side cloud admin uses Django built-in `IsAdminUser`. Cleanup tracked below. (Locked 2026-04-29; see [INVENTORY.md LOCKED #1](INVENTORY.md#locked-1-appsusers--materia).) |
| Identity primitives | **Deferred — not in shared lib for v0.** Cross-domain identity types live as Materia-internal types until empirical need arises in CCP4i2. Verified: CCP4i2 has zero construct/plasmid concept; protein sequence is `CAsuContent` / `CSeq` only. (Locked 2026-04-29; see [INVENTORY.md LOCKED #2](INVENTORY.md#locked-2-defer-identity-primitives).) |
| `mddocs/` cleanup | **Keep 6 reference docs, delete 20 process artefacts.** Materia inherits a clean docs slate. Drift-guard for the multi-inheritance pattern (`STUB_IMPLEMENTATION_INHERITANCE_PATTERN.md`) is among the keepers. (Locked 2026-04-29; see [INVENTORY.md LOCKED #3](INVENTORY.md#locked-3-mddocs-cleanup--keep-6-reference-docs-delete-20-process-artefacts).) |

## Proposal sequencing — execution status

Cross-references the 6 steps in the proposal's "Proposed sequencing" section.

### 1. Inventory ✓ (partial)

- Auth-code inventory: **done** — full audit of Python middleware (`azure_auth.py`, `azure_ad_auth.py`, `dev_auth.py` × 2) and JS-side auth (client/renderer + apps/compounds/frontend).
- Compounds-vs-CCP4i2 boundary inventory: **not done** — the broader "what moves to materia, what stays, what goes to the shared lib" exercise has only been done for auth-related code. Full inventory is a prerequisite for step 3 (the actual cut).

### 2. Shared-library proof of concept (in progress)

Living at [`packages/ccp4i2-auth/`](../../../packages/ccp4i2-auth/) on the `django-sliced` branch.

| Surface | Language | Status |
|---|---|---|
| `dev_auth.py` (was duplicated, now single-source) | Py | ✓ migrated |
| `auth-token.ts` (276 lines, was duplicated) | TS | ✓ migrated |
| `LocalSessionTokenProvider` + `hasLocalSessionToken` + email getter | TS | ✓ shipped |
| `BaseAuthMiddleware` contract | Py | ✓ shipped |
| `LocalSessionAuthMiddleware` | Py | ✓ shipped |
| `AzureADAuthMiddleware` (refactored onto base) | Py | ✓ migrated |
| `AzureADAuthentication` (DRF auth class) | Py | ✓ migrated to `ccp4i2_auth.drf` |
| `DevAdminMiddleware` (split out of AzureAD's old fallback) | Py | ✓ shipped |
| `AuthenticationFailed` / `AuthorizationFailed` exceptions | Py | ✓ shipped |
| `api-fetch` wrapper (bearer injection, 401 emit, error normalisation) | TS | ✓ migrated as `createApiFetch({ baseUrl })` factory |
| `AUTH_ERROR_EVENT` protocol | TS | ✓ migrated (universal export at module level) |
| Identity primitives (`compound-identity`, `construct-identity`) | TS | ⏳ blocked on inventory step |
| Service-contract types for CCP4i2 REST API | TS | ⏳ blocked on step 4 |
| pytest test suite | Py | ✓ scaffolded with DevAdminMiddleware smoke tests |
| vitest test suite | TS | ✓ scaffolded with createApiFetch + LocalSession smoke tests |
| CI for the package | both | ⏳ not yet added (waits on the repo cut for `ccp4/ccp4i2-auth`) |

### 3. Compounds repo cut

- `git filter-repo` plan: not yet drafted.
- New `materia` repo creation: gated on the institutional ask (CCP4 org admin to create related repo + npm/PyPI namespaces).
- Carry-along candidates: `apps/compounds/` (server + frontend), Azure deploy scripts and infrastructure specific to compounds-only instances, compounds-relevant docs.

### 4. CCP4i2 service interface documentation + contract tests

- Not started.
- Output: a versioned API contract document + TypeScript types (live in the shared lib so compounds and external integrators consume the same shape) + Django-side tests asserting the contract.
- Discussion-needed: how strictly to commit (HTTP status codes, response field stability, deprecation policy).

### 5. Lightweight auth lands in CCP4i2 desktop ✓

Done early — dovetailed with the LocalSession greenfield work in step 2.

- Per-launch token (32 random bytes hex) generated in Electron main, propagated via env to Django child + preload, exposed read-only on `window.ccp4i2LocalSession`.
- Renderer detects local session and skips MSAL initialisation; uses `setTokenGetter(createLocalSessionTokenGetter())`.
- Server-side `LocalSessionAuthMiddleware` HMAC-validates the bearer on every request, authenticates as the OS user (sanitised email).
- Cross-platform: Windows DOMAIN\user form sanitised; `.invalid` TLD avoids RFC 6761 collisions.
- **Not yet exercised**: a literal `npm run start:electron` end-to-end run. Synthetic Django RequestFactory tests + JS bundle builds prove the wiring is correct, but no live process has been run.

### 6. Standalone Materia Docker image

- Not started.
- Currently the unified web image overlays compounds onto client/renderer at Docker build time. After the cut, Materia builds its own image without that overlay.

## DDU safety contract

The `ccp4i2-bicep-*` apps in `ccp4i2-bicep-rg-uksouth` (custom domain `ddudatabase.ncl.ac.uk`) are the production CCP4i2 instance and **must not be threatened by this work**. Locked safeguards:

- Image-lineage separation: DDU pulls from `ccp4i2/{web,server}` only; `django-sliced` builds produce `materia/{web,server}` only.
- Deploy scripts default to nothing: `--env` is mandatory, so a forgotten flag errors out instead of silently aiming at DDU.
- DDU is pinned to a specific image tag (`20260429-101442` at time of writing); no moving tag follows.
- Cutover to materia is a deferred, deliberate event with rehearsed migration parity tests, an agreed maintenance window, and rollback ready. Not part of the development cadence.

## Open institutional asks

Tracked because the technical work runs ahead of these and at some point waits on them. Split into "GitHub-side asks" (require `ccp4` GitHub-org-admin = Dave Waterman / `dagewa`, the operational lead funded by STFC) and "registry-side asks" (Martin Noble can claim directly as ex-chair of the CCP4 exec committee, with co-owners added as their accounts come online).

| Ask | Right contact | State |
|---|---|---|
| Create `ccp4/ccp4i2-auth` repo on GitHub; grant maintainer access to Martin | Dave Waterman (`dagewa`) — operational org admin | initial email sent to Stuart McNicholas / Paul Bond April 2026; needs CC or follow-up to Dave (members not admins). |
| Claim `@ccp4` npm organisation; add CCP4 admins as co-owners | Martin Noble — direct (verified-free namespace; ex-chair authority) | pending; informational note to consortium when claimed |
| Claim `ccp4i2-auth` PyPI project; add CCP4 admins as co-maintainers | Martin Noble — direct | pending |
| Link `@ccp4` npm org to `ccp4` GitHub org (provenance binding) | requires both Martin (npm side) + Dave (GitHub side) | pending; sequencing — claim namespaces first, link once GitHub repo exists |
| New `materia` repo for the eventual compounds-side cut (host TBD; **not** `newcastleuniversity` per memory note) | TBD | not yet raised |

**Verified namespace state (April 2026):** `@ccp4` npm scope, `@ccp4i2` npm scope, `ccp4i2-auth` PyPI, and `ccp4i2` PyPI are all unclaimed. **`materia` (bare) is taken on both registries** — Materia-side packages will need a scope or prefix.

None of these block the in-tree workspace work. They become load-bearing only at the cut moment.

## Workstreams not in the proposal but added by this work

- **Deploy script safety** (`--env` mandatory, no `:latest`, env-controlled `IMAGE_REPO_*` in build/deploy scripts and Bicep). Shipped on `django` (not just `django-sliced`) because it protects production deploys today.
- **`materia/*` image lineage + `.env.demo-materia`**. Built and deployed to ccp4i2-demo-rg-uksouth as `materia-demo-{web,server,worker}` at tag `20260429-141934`. Validated end-to-end: refactored AzureADAuthMiddleware accepts real bearer tokens; the materia-demo web app is live at `https://materia-demo-web.agreeablebeach-f023b4ee.uksouth.azurecontainerapps.io`. DDU/kawamura `ccp4i2/*` images and apps untouched throughout.
- **`DevAdminMiddleware` split + fail-closed settings**. Closes a misconfiguration backdoor that existed pre-migration (production with `CCP4I2_REQUIRE_AUTH=false` → auto-superuser). Strictly safer.
- **`BaseAuthMiddleware` contract + `AuthenticationFailed`/`AuthorizationFailed` exceptions**. Standardises 401/403 response shape and the `_ccp4i2_auth_middleware_ran` trust flag across all auth schemes. Single point of authority for the trust contract.

## Risk register

| Risk | Mitigation | State |
|---|---|---|
| DDU production user impact | Image-lineage separation, `--env` mandatory, pinned tags, deferred cutover | implemented |
| `django-sliced` drifts too far from `django` | Forward-port discipline (weekly target) | discipline-only; no automation |
| Materia's bilingual auth contract drifts between TS and Python | Single-repo bilingual design, single contract definition | implemented |
| Cloud auth chain (refactored) untested in cloud | Demo deploy of `django-sliced` via materia/* lineage | **closed** — `materia-demo-*` deployed and signed-in-against April 2026 |
| Real Electron run with LocalSession untested | `npm run start:electron` smoke test on macOS, then Windows + Linux | macOS validated; Windows + Linux pending VM availability |
| Misconfigured cloud creates a superuser | `DevAdminMiddleware` DEBUG-gated; `settings.py` fail-closed | implemented |
| `:latest` tag re-pointed by a stray build re-introduces a moving target | `:latest` push removed from build script | implemented |
| Future renames of `apps/compounds/` make `git filter-repo` painful | Rename deferred to the moment of cut | discipline-only |

## Currently in flight

(Update this section frequently — it's the "what's on the desk right now" pointer.)

- Last commit on `django-sliced`: `594065ec3` — `Close dev_admin backdoor: split DevAdminMiddleware from AzureAD; fail-closed settings`.
- `origin/django-sliced` and `origin/django` are both up to date as of this writing.

## Recommended next actions

In rough priority order:

1. **Windows + Linux Electron smoke tests** — when VMs are available; macOS already validated.
2. **Claim `@ccp4` npm + `ccp4i2-auth` PyPI** namespaces (see Open institutional asks).
3. **Refactor `AzureADAuthMiddleware` onto `BaseAuthMiddleware`** — currently still inherits structure from the legacy module. Bringing it onto the base contract finishes the trust-flag single-sourcing across all auth middleware.
4. **Cleanup pass: remove dead `if (!res.ok) throw` blocks** at compounds-frontend call sites now that `authFetch` throws automatically. Behaviour-preserving simplification.

### Recently completed

- **Compounds-frontend `authFetch` migrated onto the createApiFetch factory** (April 2026) — 5 duplicate fetch-wrapper implementations (4 `authFetch` + 1 `coreFetch` in aggregation-api.ts, all dynamic-required `@ccp4/ccp4i2-auth`) replaced with one binding at `lib/compounds/api-fetch.ts`. `createApiFetch` grew an optional `injectUserEmail` callback so the compounds-side `X-User-Email` fallback is preserved without leaking into renderer flows. Net –100 lines. Dead `if (!res.ok) throw` blocks at call sites left in place pending a behaviour-preserving cleanup pass.
- **Test infrastructure scaffolded** (April 2026) — pytest + pytest-django on the Python side with a minimal Django bootstrap and DevAdminMiddleware smoke tests; vitest + jsdom on the TypeScript side with smoke tests for `createApiFetch`, `AUTH_ERROR_EVENT`, and the LocalSession provider helpers. 12 tests, both runners green in <1s. Future test additions are mechanical drop-ins.
- **`api-fetch.ts` migrated** (April 2026) — refactored into a `createApiFetch({ baseUrl })` factory in `packages/ccp4i2-auth/src/api-fetch.ts`. `client/renderer/api-fetch.ts` is now a 39-line thin binding to the factory, so all 21 import sites in the renderer keep working unchanged. The factory is consumer-agnostic; compounds-frontend can adopt the same canonical implementation by binding to its own baseUrl. Closes the JS-side v0 surface of the shared library.
- **Demo deploy validation** (April 2026) — `materia-demo-{web,server,worker}` deployed to `ccp4i2-demo-rg-uksouth` at image tag `20260429-141934`. Refactored AzureADAuthMiddleware accepts real Azure AD bearer tokens; sign-in flow works end-to-end. DDU and kawamura instances entirely untouched (separate `ccp4i2/*` lineage). Two small deploy-time fixes landed alongside: `@types/node` added to ccp4i2-auth devDependencies (Docker isolation surfaced what npm hoisting hid locally); `materia-demo-web.*` URLs added to the shared `ccp4i2-demo` AAD app registration's allowed redirect URIs (Azure portal). The shared AAD reg is a transition shortcut; revisit before broader cross-tenant materia user testing.
- **macOS Electron LocalSession smoke test** (April 2026) — `npm run start` on macOS, click "Launch CCP4i2", projects load via LocalSession-authenticated requests. Caught and fixed a layout-level bug: `app/layout.tsx` was gating MsalAuthProvider on `NEXT_PUBLIC_REQUIRE_AUTH`, which doesn't apply in Electron dev. Layout now also mounts AuthProvider when `CCP4I2_LOCAL_SESSION_TOKEN` is in env (visible during SSR via the Electron-spawned Next.js server's process.env).
- **CCP4i2 dependency on `apps/users/` dropped** (per LOCKED #1). `server/ccp4i2/api/admin_views.py` now uses `rest_framework.permissions.IsAdminUser`; the conditional `users` URL inclusion block in `api/urls.py` deleted. CCP4i2 cloud admin gating is now Django-built-in only; no role-system dependency.

After each next-action lands, update the "Status snapshot" table at the top of this document and the relevant workstream below.

## How to use this document

- **Starting a new session**: read the Status snapshot + Locked decisions + Currently in flight. Skim the rest.
- **Finishing a meaningful chunk of work**: update the relevant table row(s) here, in the same commit (or an immediately following one). Future you should be able to read this doc and know what's true.
- **Considering an architectural change**: read both proposal documents and confirm the change is consistent with their framing. If the change overrides a Locked decision, log the reason here and update the row.
