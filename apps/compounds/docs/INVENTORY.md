# CCP4i2 / Materia inventory

*The fault line. Where does each part of the codebase belong after the slice?*

*Companion to [REPO_SPLIT_IMPLEMENTATION.md](REPO_SPLIT_IMPLEMENTATION.md). This document is the **draft**: most cases are settled and verified, but a small number of genuinely ambiguous cases are flagged with **DECISION REQUIRED** and need stakeholder input before they can be locked.*

## Methodology

Every component falls into one of five categories:

| Category | Test |
|---|---|
| **CCP4i2** | Serves crystallography. Referenced only by crystallography code. |
| **Materia** | Serves compounds / assays / drug discovery. Not used by crystallography. |
| **Shared lib** | Divergence between sides would be a *bug* — single canonical contract required. |
| **Ambiguous** | Genuinely cross-domain or used by both sides; needs explicit decision. |
| **Deletable** | Unused, abandoned, or made redundant by migrations. |

**Tie-breaker for ambiguous cases:** *small identity primitive in shared lib + rich machinery in whichever side currently owns it.* Identity flows freely; rich machinery stays put. Per the proposal's "Cross-domain entities" section.

**Discipline:** every entry gets a one-line *why*. For ambiguous cases, the options and recommendation. Decisions get cross-linked back to [REPO_SPLIT_IMPLEMENTATION.md](REPO_SPLIT_IMPLEMENTATION.md)'s Locked decisions table once locked.

## Summary table

Verified by grep + visual inspection on `django-sliced` at commit `b5e4c3fd0`.

| Path | Category | Why |
|---|---|---|
| `server/ccp4i2/api/` | CCP4i2 | REST API for jobs, projects, files. Compounds URL inclusion is conditional (`COMPOUNDS_ENABLED`). |
| `server/ccp4i2/db/` | CCP4i2 | ORM models for projects, files, jobs. No compounds FKs. |
| `server/ccp4i2/middleware/` | CCP4i2 | Just `corp.py` (Moorhen WASM CORP) left after the auth migration. The folder stays for the corp middleware. |
| `core/` | CCP4i2 | Crystallography base classes, task manager. |
| `wrappers/`, `wrappers2/` | CCP4i2 | ~116 program wrappers (Refmac, Phaser, Coot, etc.). |
| `pipelines/` | CCP4i2 | Composite crystallographic workflows. |
| `pimple/` | CCP4i2 | Plugin/manifest system. |
| `smartie/` | CCP4i2 | Log parsing utilities. |
| `cli/` | CCP4i2 | i2run CLI. |
| `report/` | CCP4i2 | Crystallographic report generation. |
| `tipsOfTheDay/` | CCP4i2 | Crystallography tips. |
| `qticons/`, `svgicons/` | CCP4i2 | Task icons (Qt legacy + SVG). |
| `assets/` | CCP4i2 | (under `client/`) Crystallographic UI assets. |
| `client/main/` | CCP4i2 | Electron main process. |
| `client/preload/` | CCP4i2 | Electron preload (now also exposes LocalSession surface for shared auth). |
| `client/renderer/` | CCP4i2 | Next.js renderer (excluding compounds/* subtrees that live in `apps/compounds/frontend/`). |
| `Docker/server/` | CCP4i2 | Server image (now COPYs `packages/ccp4i2-auth/` too). |
| `Docker/client/Dockerfile.cross-platform` | **Mostly CCP4i2 + deletable lines** | Lines 80-120 are the per-file overlay of compounds frontend; deleted at the cut. The rest is the CCP4i2 web build. |
| `Docker/base/` | CCP4i2 | Base image (CCP4 + ARP/wARP). |
| `Docker/azure-uksouth/scripts/` | CCP4i2 | Generic deploy machinery (build-and-push, deploy-applications, etc.). All require `--env`; image lineage parameterised. |
| `Docker/azure-uksouth/.env.deployment` | CCP4i2 | DDU production env. Image lineage `ccp4i2/{web,server}`. |
| `Docker/azure-uksouth/.env.demo` | CCP4i2 | ccp4i2-demo (stable) env. |
| `Docker/azure-uksouth/.env.kawamura` | CCP4i2 | Kawamura instance env. |
| `Docker/azure-uksouth/.env.demo-materia` | Materia | Demo env for materia development. Moves to materia repo. |
| `Docker/azure-uksouth/infrastructure/` | CCP4i2 | Bicep templates. |
| `Docker/azure-uksouth/azure_extensions/` | CCP4i2 | Azure-specific Django extensions. |
| `Docker/azure-uksouth/legacy/` | likely **Deletable** | Pre-bicep deployment configs (per `azure_ad_auth.py` precedent — example/reference code). Verify before deletion. |
| `Docker/azure-uksouth/docs/` | CCP4i2 | Deploy operations docs. |
| `Docker/azure-uksouth/AZURE_AD_SETUP.md` etc. | CCP4i2 | Deployment runbooks. |
| `Docker/scripts/` | CCP4i2 | `entrypoint.sh`, `startup-server.sh`, `startup-worker.sh`. |
| `apps/compounds/registry/` | Materia | Compound + Batch + Lab Notebook + Supplier + Target. Compounds platform. |
| `apps/compounds/assays/` | Materia | Assay protocols, fitting, data series, analysis. |
| `apps/compounds/constructs/` | Materia (rich) **+ Shared lib (identity primitive)** | Plasmids, cassettes, sequencing results stay in Materia. Identity type extracted to shared lib so CCP4i2 can reference constructs. See **DECISION REQUIRED #2**. |
| `apps/compounds/nlp/` | Materia | NLP query, scaffold extension, substructure catalog. Compounds cheminformatics. |
| `apps/compounds/frontend/` | Materia | Next.js compounds app (the standalone build target). |
| `apps/compounds/scripts/` | Materia | Compounds-specific management commands. |
| `apps/compounds/tests/` | Materia | Compound-specific tests. |
| `apps/compounds/docs/` | Materia | Compounds docs (proposals, NLP, this inventory, the implementation plan). Move with materia. |
| `apps/compounds/admin_views.py`, `media_views.py`, `config_views.py`, etc. | Materia | Compounds-specific Django views. |
| `apps/compounds/settings.py` | Materia | Compounds Django app settings. |
| `apps/compounds/urls.py` | Materia | Compounds URL routing. |
| `apps/compounds/utils.py`, `formatting.py`, `validators.py` | Materia | Compounds-specific helpers. |
| `apps/users/` | **Materia** (locked) | The desktop path never touches it (`users` not in default `INSTALLED_APPS`). The 3 CCP4i2-side imports are cloud admin features that switch to Django built-in `IsAdminUser` as part of the cleanup. See locked decision below. |
| `packages/ccp4i2-auth/` | Shared lib | ✓ Already extracted. Bilingual TS+Python. CCP4-stewarded. |
| `tests/` (root, `ccp4i2/tests/`) | CCP4i2 | Existing test suite — covers crystallography. (Root `tests/` is empty; real tests live at `server/ccp4i2/tests/`.) |
| `docs/` | CCP4i2 | Crystallography documentation. |
| `mddocs/` | **AMBIGUOUS — DECISION REQUIRED #3** | Planning notes (API harmonisation, base class decisions, async changelog). Some clearly CCP4i2; some may be cross-cutting. |
| `migration/` | **Deletable** | `migration/CData/` only. No `from migration` imports anywhere. Stale. |
| `scripts/generate_task_interfaces.py` | CCP4i2 | Generates React task interfaces from def.xml. Crystallography. |
| `scripts/backfill_target_genes.py` | Materia (or **Deletable**) | One-off Materia operator tool. References `Target.genes`, talks to deployed REST API. Likely move-with-materia or delete after one final use. |

## Decisions required

The cases below are not unilaterally settleable from a code reading. They need stakeholder input.

### LOCKED #1: `apps/users/` → Materia

**Decision** (locked April 2026): `apps/users/` is Materia-resident. The CCP4i2 side has no dependency on it; cloud-CCP4i2 admin features use Django's built-in `IsAdminUser` (`is_staff` / `is_superuser`) for coarse admin gating.

**Rationale:** Verification showed the Electron desktop path never touches `apps/users/` — it's not in the default `INSTALLED_APPS` ([settings.py:44-53](../../../server/ccp4i2/config/settings.py#L44-L53)) and the URLs are conditionally included only when the app is. The 3 CCP4i2-side import sites are all cloud-admin features (legacy-import endpoints in `admin_views.py`); the 8+ import sites are all Materia-side (compounds admin, assays, constructs, registry). The role/permission system is genuinely a multi-user-platform concern, which is the Materia case, not the desktop case.

**Cleanup work tracked in REPO_SPLIT_IMPLEMENTATION.md:**

| Site | Action |
|---|---|
| [`server/ccp4i2/api/admin_views.py:20`](../../../server/ccp4i2/api/admin_views.py#L20) and 4 `@permission_classes([IsPlatformAdmin])` decorators | Replace with `rest_framework.permissions.IsAdminUser` |
| [`server/ccp4i2/api/urls.py:84-94`](../../../server/ccp4i2/api/urls.py#L84) | Delete the conditional `users` URL inclusion block |
| Comments referencing the role system in cloud-CCP4i2 settings | Update or remove |

After this cleanup, `apps/users/` migrates with Materia at the `git filter-repo` cut. No second shared library is created — `ccp4i2-users` does not exist.

**Trade-off:** any future pure-CCP4i2 cloud deployment that wants the CONTRIBUTOR role gets it back by either installing Materia's users module as a pip dependency or implementing its own. Realistically no such deployment exists; CCP4i2 cloud presence is the unified DDU-style deployment that keeps Materia in the build.

**Options considered and rejected:**
- Extract to a second shared library (`ccp4i2-users`) — ruled out because the desktop doesn't need users at all; making it a shared dependency burdens the desktop case unnecessarily.
- Keep in CCP4i2 with Materia as consumer — bilateral-stewardship risk; CCP4i2 changes could break Materia silently.
- Fork at the cut — divergence drift; users with both apps would see inconsistent roles.

---

### DECISION REQUIRED #2: location of identity primitives

Concrete primitives we'll likely want in shared form (per the "Cross-domain entities" section of the proposal):

- `ConstructIdentity` — `id`, `formatted_id` (e.g. `NCLCON-12345678`), `name`, optional `protein_uniprot_id`, optional cassette range. Lets CCP4i2 reference a construct by identity without depending on Materia's full registry.
- `CompoundIdentity` (perhaps later) — `id`, `compound_code`, optional `smiles`. Lets CCP4i2 record which compound a crystallographic experiment used, without dragging in the registry.
- `CampaignIdentity` (perhaps later) — `id`, `name`. Cross-domain campaign references.

**Where do these live?**

1. **In `packages/ccp4i2-auth/`** — same package, broader scope. Risk: name-mismatch ("auth" undersells what's there). Pro: one package to steward, one CI to watch. (The proposal noted this option already.)

2. **In a new `packages/ccp4i2-entities/`** — cleaner separation. Risk: another package, another publish/version cycle.

**My recommendation:** start by adding `ConstructIdentity` to `ccp4i2-auth/` — the surface is small, and we don't yet know how many other identity types will materialise. Split into a separate `entities` package later if (a) the entity-types section grows beyond ~5 types, or (b) third parties want entities without the auth dependency. Defer the rename/split decision until we see how it grows.

**Locking this:** if you agree, I'll update REPO_SPLIT_IMPLEMENTATION.md's Locked decisions table to add: *"Identity primitives go in `packages/ccp4i2-auth/` for v0; consider extracting later if catalogue exceeds ~5 types."*

---

### DECISION REQUIRED #3: `mddocs/`

`mddocs/` contains planning notes:
- `API_HARMONIZATION_PLAN.md`, `API_HARMONIZATION_PROGRESS.md`
- `BASE_CLASS_DECISION.md`
- `CDATAFILE_STUB_ANALYSIS.md`
- `CHANGELOG_ASYNC.md`
- (and more — full listing in the directory)

Some are clearly CCP4i2 historical (CDataFile, async work). Some may be cross-cutting (API harmonisation could be the cross-domain API work).

**Options:**

1. **Treat all as CCP4i2** — they live under `mddocs/` at the repo root, predate the materia split, are about the crystallography work.
2. **Audit per-file** — slow, but might surface a couple that should travel with materia.
3. **Mark deletable if older than the materia work** — these are historical decision logs; their relevance has often expired.

**My recommendation:** **Option 1 — keep all of them in CCP4i2 by default.** They're historical artefacts of CCP4i2's evolution; even cross-cutting ones (API harmonisation) are written from CCP4i2's perspective. Materia inherits a clean docs slate at the cut; if any of these are relevant they can be referenced or copied by hand at the time. Avoids per-file audit overhead.

Want me to spot-check one or two before locking?

---

## Deletable candidates (verified)

| Path | Evidence | Action |
|---|---|---|
| `migration/CData/` | No `from migration` imports anywhere in `server/` or `apps/`. Single subdirectory, looks pre-Django. | Delete. |
| `Docker/client/Dockerfile.cross-platform` lines 80-120 (compounds overlay) | Lines explicitly per-file COPY of compounds frontend onto client/renderer. | Delete **at the cut** (not now — the demo build still uses the unified image). |
| `Docker/client/Dockerfile` (non-cross-platform variant) lines copying compounds frontend | Same as above. | Delete at the cut. |
| `Docker/azure-uksouth/legacy/` | Naming + the existing `azure_ad_auth.py` precedent suggests this is reference/example code. | Verify imports before deletion. |
| `apps/compounds/dev_auth.py` | Already deleted in commit `392203a49`. | (Done.) |
| `server/ccp4i2/middleware/azure_ad_auth.py` | Already deleted in commit `0c014ccff`. | (Done.) |
| `server/ccp4i2/middleware/azure_auth.py` | Already migrated to `packages/ccp4i2-auth/`. | (Done.) |

## Cross-cutting concerns

### The Docker overlay machinery

[Docker/client/Dockerfile](../../../Docker/client/Dockerfile) and [Dockerfile.cross-platform](../../../Docker/client/Dockerfile.cross-platform) currently build a unified image that overlays compounds frontend (`apps/compounds/frontend/`) onto `client/renderer/` via per-file COPY. Post-cut:

- **CCP4i2 side:** the overlay block is deleted from `Docker/client/Dockerfile*`. The client image becomes a clean crystallography-only Next.js build.
- **Materia side:** Materia repo brings its own `Dockerfile` that builds `apps/compounds/frontend/` (or whatever it's called by then) directly. No overlay.
- **Transition:** during `django-sliced` development, the overlay still works — the `materia/*` image lineage is built from this same Dockerfile (for now).

### The `materia/*` image lineage

Locked safeguard: DDU and kawamura pull only `ccp4i2/{web,server}`; demo-materia pulls only `materia/{web,server}`. Image-repo names are env-controlled (`IMAGE_REPO_WEB`, `IMAGE_REPO_SERVER`), so a script aimed at DDU cannot push or pull a materia tag. Verified in [REPO_SPLIT_IMPLEMENTATION.md DDU safety contract](REPO_SPLIT_IMPLEMENTATION.md#ddu-safety-contract).

### UI components currently shared via overlay

Components that compounds-side currently borrows from `client/renderer/` (`DataTable`, `SearchField`, `RequireAuth`, theme tokens, snackbar styling): per the proposal, **fork at the cut, no shared library**. Divergence is a *feature* for UI; bug fixes do not auto-propagate, but the trade is worth it for design autonomy. No decision needed — this is locked in the proposal.

## Identity primitives — proposed shared-lib additions

(Subject to DECISION REQUIRED #2.)

| Type | Lives in | Minimal fields | Why shared |
|---|---|---|---|
| `ConstructIdentity` | `packages/ccp4i2-auth/src/entities/` (or similar) | `id`, `formatted_id`, `name`, `protein_uniprot_id?`, cassette range | CCP4i2 can reference a construct without depending on Materia's full registry. Both sides agree on the wire format. |
| `CompoundIdentity` (later) | same | `id`, `compound_code`, `smiles?` | If CCP4i2 ever records "which compound this experiment used", reference by identity. |
| `CampaignIdentity` (later) | same | `id`, `name` | If campaigns become cross-domain. |

Add these as the first non-auth contents of the shared library. Names can be tweaked at PR time.

## Open questions (further investigation)

1. Does `Docker/azure-uksouth/legacy/` actually contain reference code, or is something still consuming it? Quick `grep` would settle.
2. Does any compound code in `apps/compounds/` import from a `client/renderer/lib/compounds/` path at runtime (not just at Docker overlay time)? If yes, that's a real cross-imports we'd need to handle. The dynamic `require()` calls in compounds-frontend that hit `@ccp4/ccp4i2-auth` are already handled.
3. Are there any imports of `apps.compounds.X` from `server/ccp4i2/`? (Should be zero; the dependency direction is one-way.)

---

## Status

This document is a **draft** until DECISIONS REQUIRED #1, #2, and #3 are locked. Once locked, it moves into REPO_SPLIT_IMPLEMENTATION.md's Locked decisions table and the working plan's "Compounds repo cut" workstream gets its concrete next-step list.
