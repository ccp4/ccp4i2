# CCP4i2 / Materia inventory

*The fault line. Where does each part of the codebase belong after the slice?*

*Companion to [REPO_SPLIT_IMPLEMENTATION.md](REPO_SPLIT_IMPLEMENTATION.md). All categorisations locked April 2026. The **LOCKED #1**, **LOCKED #2**, **LOCKED #3** sections below explain the genuinely difficult cases and their decision rationale so future-you can re-read the reasoning.*

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
| `apps/compounds/constructs/` | Materia | Plasmids, cassettes, sequencing results all live in Materia. CCP4i2 has no construct concept currently (verified — see LOCKED #2), so no shared identity primitive needed for v0. |
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
| `mddocs/` | CCP4i2 (cleaned up) | 6 reference docs kept; 20 process artefacts deleted. See LOCKED #3 for the keep/delete split and rationale. |
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

### LOCKED #2: defer identity primitives

**Decision** (locked April 2026): no identity primitives in the shared library for v0. The shared lib stays auth-focused. Cross-domain identity types (`ConstructIdentity`, `CompoundIdentity`, `CampaignIdentity` etc.) live as Materia-internal types until empirical need arises.

**Rationale:** the proposal anticipated CCP4i2 wanting to reference Materia entities without dragging in Materia's full registry — and identity primitives in the shared lib were the response. **Empirical verification ([CCP4ModelData.py](../../../server/ccp4i2/core/CCP4ModelData.py)) shows CCP4i2 has no concept of constructs, plasmids, or expression vectors.** Protein-sequence handling is entirely `CAsuContent` (line 13), `CAsuContentSeq` (124), `CSequence` (2620), `CSeqDataFile` (2596), and `CSequenceAlignment` (2840) — crystallographic primitives. Zero "construct" or "plasmid" references in `core/`, `server/ccp4i2/`, `wrappers/`, `wrappers2/`, `pipelines/`, `pimple/`, `smartie/`, or `cli/`. The need the proposal anticipated has not manifested in code.

**Implication:** designing identity primitives speculatively, before a real reference exists, risks:
- A speculative API that doesn't match real usage when the need does arise.
- Turning the auth library into a "miscellaneous shared types" dumping ground.
- Carrying a maintenance burden for types nobody consumes.

**When to revisit:** if/when CCP4i2 grows a feature that wants to reference a Materia entity (e.g., a future `Project.construct_id` foreign key, or a crystallographic dataset annotated with the construct it came from), we add the identity primitive at that point — with the actual usage informing the shape rather than designing speculatively. The shared library's contract surface is one-way extensible; nothing about the current shape closes the door.

**This is an explicit divergence from the proposal**, recorded so future-you can re-read this and know why. The proposal's framing is forward-looking; the verification shows the forward-looking concern hasn't materialised yet.

---

### LOCKED #3: `mddocs/` cleanup — keep 6 reference docs, delete 20 process artefacts

**Decision** (locked April 2026): the `mddocs/` folder is mostly cruft — planning notes, progress reports, decision logs, completion summaries. Most have served their purpose; their content is either applied to the codebase or superseded. Triaged into 6 keepers (live reference material) and 20 deletions (process artefacts).

**Keep (6 files):**

| File | Why |
|---|---|
| `README.md` | Folder index (rewritten to drop refs to deleted files) |
| `PLUGIN_REGISTRY_README.md` | Live system documentation for the plugin registry |
| `STUBS_README.md` | Live documentation for PySide2/Qt stub modules used in plugin discovery |
| `STUB_IMPLEMENTATION_INHERITANCE_PATTERN.md` | **Drift-guard** for the multi-inheritance pattern (`class CObsDataFile(CObsDataFileStub, CMiniMtzDataFile)`). Prescriptive, captures non-obvious MRO concerns; deleting it would invite future contributors to re-derive the wrong way first. |
| `qt_task_gui_guide.md` | 1471-line comprehensive guide; recently maintained (2026-02-28) |
| `QUICK_REFERENCE.md` | Async execution infrastructure reference (recently maintained) |

**Delete (20 files):** all naming-pattern matches for `_PLAN`, `_PROGRESS`, `_DECISION`, `_ANALYSIS`, `_CHANGELOG`, `_COMPLETE`, `_SUCCESS`, `_SUMMARY`, `_APPLIED`, `_MILESTONE`, `_REFACTOR`, `_MIGRATION`, `_INTEGRATION`. Plus `MTZ_CONVERSION_SYSTEM.md` (status doc with `"Conversion logic pending"` marker — half-done refactor snapshot). Git history retains all of them; nothing is actually lost.

**Rationale:** Materia inherits a clean docs slate. CCP4i2 doesn't carry forward a folder of stale planning notes. If any specific deleted doc proves load-bearing later, it's recoverable from git. The keepers are exactly the docs that would be hard to re-derive (drift-guard, comprehensive guides, system READMEs).

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

## Identity primitives — deferred (not in v0)

Per LOCKED #2: no identity primitives in the shared library for v0. Cross-domain identity types live as Materia-internal types until empirical need arises. When CCP4i2 grows a feature that wants to reference a Materia entity (e.g., a future `Project.construct_id` foreign key), we add the primitive at that point with the actual usage informing its shape.

## Open questions (further investigation)

1. Does `Docker/azure-uksouth/legacy/` actually contain reference code, or is something still consuming it? Quick `grep` would settle.
2. Does any compound code in `apps/compounds/` import from a `client/renderer/lib/compounds/` path at runtime (not just at Docker overlay time)? If yes, that's a real cross-imports we'd need to handle. The dynamic `require()` calls in compounds-frontend that hit `@ccp4/ccp4i2-auth` are already handled.
3. Are there any imports of `apps.compounds.X` from `server/ccp4i2/`? (Should be zero; the dependency direction is one-way.)

---

## Status

All three LOCKED decisions are recorded above and reflected in [REPO_SPLIT_IMPLEMENTATION.md](REPO_SPLIT_IMPLEMENTATION.md)'s Locked decisions table. Proposal step 1 (Inventory) is now complete; the "Compounds repo cut" workstream can pick up next-step actions from the Recommended next actions list in the implementation plan.
