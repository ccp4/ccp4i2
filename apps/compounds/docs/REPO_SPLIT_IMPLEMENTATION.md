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
| **Lightweight desktop auth** | shipped — LocalSession end-to-end |
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
| `api-fetch` wrapper (bearer injection, 401 emit, error normalisation) | TS | ⏳ not yet migrated |
| `AUTH_ERROR_EVENT` protocol | TS | ⏳ not yet migrated |
| Identity primitives (`compound-identity`, `construct-identity`) | TS | ⏳ blocked on inventory step |
| Service-contract types for CCP4i2 REST API | TS | ⏳ blocked on step 4 |
| pytest test suite | Py | ⏳ not yet added |
| vitest test suite | TS | ⏳ not yet added |
| CI for the package | both | ⏳ not yet added |

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

Tracked because the technical work runs ahead of these and at some point waits on them.

| Ask | Stakeholder | State |
|---|---|---|
| Create `ccp4/ccp4i2-auth` repo on GitHub | `ccp4` org admin (Stuart McNicholas / Paul Bond contacted) | email sent (April 2026), awaiting reply |
| Create `@ccp4` npm organisation with publish rights for the auth lib | npm registry / `ccp4` org admin | downstream of repo creation |
| PyPI account / project for `ccp4i2-auth` | TBD | downstream of repo creation |
| New `materia` repo (host TBD; **not** under `newcastleuniversity` per memory note) | TBD | not yet raised |

None of these block the in-tree workspace work. They become load-bearing only at the cut moment.

## Workstreams not in the proposal but added by this work

- **Deploy script safety** (`--env` mandatory, no `:latest`, env-controlled `IMAGE_REPO_*` in build/deploy scripts and Bicep). Shipped on `django` (not just `django-sliced`) because it protects production deploys today.
- **`materia/*` image lineage + `.env.demo-materia`**. Wired but **not yet built or deployed**. First validation moment is the next demo deploy.
- **`DevAdminMiddleware` split + fail-closed settings**. Closes a misconfiguration backdoor that existed pre-migration (production with `CCP4I2_REQUIRE_AUTH=false` → auto-superuser). Strictly safer.
- **`BaseAuthMiddleware` contract + `AuthenticationFailed`/`AuthorizationFailed` exceptions**. Standardises 401/403 response shape and the `_ccp4i2_auth_middleware_ran` trust flag across all auth schemes. Single point of authority for the trust contract.

## Risk register

| Risk | Mitigation | State |
|---|---|---|
| DDU production user impact | Image-lineage separation, `--env` mandatory, pinned tags, deferred cutover | implemented |
| `django-sliced` drifts too far from `django` | Forward-port discipline (weekly target) | discipline-only; no automation |
| Materia's bilingual auth contract drifts between TS and Python | Single-repo bilingual design, single contract definition | implemented |
| Cloud auth chain (refactored) untested in cloud | Demo deploy of `django-sliced` via materia/* lineage | open — recommended next |
| Real Electron run with LocalSession untested | `npm run start:electron` smoke test on macOS, then Windows + Linux | open |
| Misconfigured cloud creates a superuser | `DevAdminMiddleware` DEBUG-gated; `settings.py` fail-closed | implemented |
| `:latest` tag re-pointed by a stray build re-introduces a moving target | `:latest` push removed from build script | implemented |
| Future renames of `apps/compounds/` make `git filter-repo` painful | Rename deferred to the moment of cut | discipline-only |

## Currently in flight

(Update this section frequently — it's the "what's on the desk right now" pointer.)

- Last commit on `django-sliced`: `594065ec3` — `Close dev_admin backdoor: split DevAdminMiddleware from AzureAD; fail-closed settings`.
- `origin/django-sliced` and `origin/django` are both up to date as of this writing.

## Recommended next actions

In rough priority order:

1. **Drop CCP4i2's dependency on `apps/users/`** (per [INVENTORY.md LOCKED #1](INVENTORY.md#locked-1-appsusers--materia)). Concrete:
   - In [server/ccp4i2/api/admin_views.py](../../../server/ccp4i2/api/admin_views.py): replace `from users.permissions import IsPlatformAdmin` with `from rest_framework.permissions import IsAdminUser`, and update the 4 `@permission_classes([IsPlatformAdmin])` decorators to `[IsAdminUser]`.
   - In [server/ccp4i2/api/urls.py](../../../server/ccp4i2/api/urls.py) lines 84-94: delete the conditional `users` URL inclusion block.
   - Verify ccp4-python imports cleanly with `users` removed from the path.
   - Lands on `django-sliced`. Doesn't need to forward-port to `django` immediately (DDU still has `users` overlaid in its image).
2. **Deploy `django-sliced` → ccp4i2-demo via `materia/*` image lineage.** Validates the AzureAD migration in cloud + exercises the materia/* lineage that has never been built. Catches Docker-side regressions before further migrations compound them. Wallclock ~30-60 min including build.
3. **Real Electron smoke test** — `npm run start:electron` on macOS, watch logs end-to-end. Then ideally same on Windows and Linux, but macOS is the immediate confidence-builder.
4. **Migrate `api-fetch.ts`** (457-line wrapper in client/renderer) into `packages/ccp4i2-auth/`. Completes the JS-side v0 surface and ends the next chunk of soon-to-be-duplicated code.
5. **Add test infrastructure** — pytest for the Python side, vitest for the TS side, both in `packages/ccp4i2-auth/tests/`.

After each next-action lands, update the "Status snapshot" table at the top of this document and the relevant workstream below.

## How to use this document

- **Starting a new session**: read the Status snapshot + Locked decisions + Currently in flight. Skim the rest.
- **Finishing a meaningful chunk of work**: update the relevant table row(s) here, in the same commit (or an immediately following one). Future you should be able to read this doc and know what's true.
- **Considering an architectural change**: read both proposal documents and confirm the change is consistent with their framing. If the change overrides a Locked decision, log the reason here and update the row.
