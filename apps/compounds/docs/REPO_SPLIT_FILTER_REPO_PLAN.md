# Repo split — `git filter-repo` execution plan

*Draft. Companion to [REPO_SPLIT_IMPLEMENTATION.md](REPO_SPLIT_IMPLEMENTATION.md). The mechanics of step 3 of [the proposal sequencing](CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md) — cutting the compounds-side carve-out into its own repo — written far enough in advance to be reviewed before execution.*

## Decisions locked into this plan

| Decision | Value | Why |
|---|---|---|
| **Materia repo home** | `newcastleuniversity/materia` | Institutional grounding aligns with the sustainability narrative; user has SSO + create access. CCP4i2-facing assets stay on `ccp4` org. |
| **History preservation** | Full history of carved paths preserved (no shallow squash) | `git blame` / `git log apps/compounds/` survives the cut for personal continuity. The carved paths are essentially single-author so there's no parity-of-attribution concern beyond personal traceability. |
| **Materia's relationship to the contract** | External-slash-internal canary | Post-cut, Materia consumes `@ccp4/ccp4i2-auth` from npm + `ccp4i2-auth` from PyPI like any third party. Its CI is the empirical proof the contract is sufficient. No in-tree shortcuts. |
| **CCP4i2-side carve-out copies** | Read-only after cut, deleted later | `apps/compounds/` and `apps/users/` continue to exist on `ccp4/ccp4i2:django-sliced` for the transition period (so DDU and demo deploys keep working from the same tree); deletion follows once Materia repo proves itself. |

## What carves out (goes to `newcastleuniversity/materia`)

Comprehensive list of paths to keep in the filter-repo invocation. Verified against current `django-sliced` tree.

### Server-side

- `apps/compounds/` — assays, registry, constructs, nlp, admin_views, settings, urls, etc. (the entire compounds Django app)
- `apps/users/` — Materia-resident per [LOCKED #1](INVENTORY.md#locked-1-appsusers--materia). Already removed from CCP4i2 server's URL conf (commit `fc5a48d09`).

### Frontend

- `apps/compounds/frontend/` (already nested under apps/compounds, included by the path above) — Next.js routes, shared components/, lib/, types/. Note: `frontend/node_modules/` is `.gitignore`d so it doesn't contribute to history weight.

### Docker / infra

- `Docker/azure-uksouth/.env.demo-materia` — the materia-demo deploy target, materia-only by design.
- *Selected* scripts under `Docker/azure-uksouth/scripts/` that exclusively target materia/* lineage. Most scripts are shared (`build-and-push.sh` accepts `--env`, used by both lineages); those stay shared in CCP4i2 and are *re-implemented* in Materia rather than carved (see "What forks at the cut" below).

### Docs

- `apps/compounds/docs/` — all of `CCP4I2_RELATIONSHIP_AND_SUSTAINABILITY.md`, `CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md`, `REPO_SPLIT_IMPLEMENTATION.md`, `CCP4I2_SERVICE_CONTRACT.md`, `INVENTORY.md`, this plan, plus the Materia-internal proposals (HELM, inventory, etc.). Note: the service contract doc is a Materia-side document about the *CCP4i2 contract Materia consumes* — it's correctly Materia-resident; CCP4i2's view of its own contract lives in its own README plus the `@ccp4/ccp4i2-auth` README.

### Strict exclusions (do **not** carve)

- `client/renderer/` — stays in CCP4i2 entirely.
- `core/`, `wrappers/`, `wrappers2/`, `pipelines/`, `pimple/`, `smartie/`, `report/` — CCP4i2 crystallographic logic, not Materia's concern.
- `server/` — Django backend, including the api/ contract guards (`tests/api/unit/test_contract.py`). The contract guards live with the *server they guard*, not the consumer.
- `packages/ccp4i2-auth/` — bound for its own repo (`ccp4/ccp4i2-auth`), not Materia.
- `Docker/azure-uksouth/.env.deployment` (template), `.env.demo`, `.env.kawamura` — CCP4i2-side cloud lineage.
- DDU-targeted infra (`ccp4i2-bicep-*`).

## What forks at the cut (duplicate-and-diverge)

Locked principle from [the proposal](CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md): bilateral stewardship is for things where divergence is a *bug* (auth, identity, service contracts); forking is fine where divergence is a *feature*. At the moment of the cut, Materia takes copies of the following from CCP4i2 and both repos own their own copy thereafter:

- **UI primitives** that compounds-frontend uses from `client/renderer/` (`DataTable`, `SearchField`, `RequireAuth`, theming utilities, etc.). Compounds-frontend has been importing these via relative `../../../client/renderer/` paths; at the cut, those copy into Materia and the imports rewrite to be Materia-relative.
- **Deploy scripts** under `Docker/azure-uksouth/scripts/` (`build-and-push.sh`, `deploy-applications.sh`, etc.) — shared today via the `--env <file>` indirection; at the cut, Materia takes its own copy and may simplify (it only needs the materia/* lineage codepath).
- **Bicep templates** for container-app provisioning — likewise duplicated, Materia owns its own from the cut moment.
- **Anything else surfaced during inventory** — the proposal puts the bar at "is divergence a bug or a feature here". Auth contract: bug. UI button styling: feature. Anything in between gets reasoned about explicitly.

## What replaces the in-tree dependency

After the cut, Materia's `package.json` and Python deps reference the published artifacts:

- `@ccp4/ccp4i2-auth` from npm — replaces all imports currently going through `packages/ccp4i2-auth/` workspace symlinks. Verified pre-cut by publishing a real version of the lib and pointing Materia's frontend at it; if everything still builds and runs, the workspace symlink wasn't hiding anything.
- `ccp4i2-auth` from PyPI — replaces server-side imports. Materia's compounds Django middleware stack pulls `LocalSessionAuthMiddleware` / `AzureADAuthMiddleware` / `DevAdminMiddleware` / `BaseAuthMiddleware` from there.

This dependency direction is the load-bearing test of the canary framing: if Materia's CI passes against the published lib, the contract is sufficient for any external consumer.

## Materia's overlay → standalone build

CCP4i2's web Docker image currently overlays the compounds frontend onto `client/renderer/` at Docker build time (see `Docker/client/Dockerfile` for the COPY-per-file pattern, source of recurring drift). After the cut:

- **Materia builds its own web image** standalone — no overlay step, just a normal Next.js build.
- **CCP4i2's web image stops including the overlay** — `Docker/client/Dockerfile` loses the compounds-related COPY block. The CCP4i2 web image becomes pure renderer (returns to its pre-overlay shape, modulo any other intentional shared additions).
- **DDU keeps pulling `ccp4i2/{web,server}` lineage** unchanged — what changes is what those images contain, not where DDU pulls from.
- **Materia gets `materia/{web,server}` lineage** (already used for development today) — published from Materia's own CI.

## `git filter-repo` invocation (sketch)

```bash
# In a fresh clone of ccp4/ccp4i2 dedicated to the cut (DO NOT do this in
# the working tree — filter-repo rewrites refs).
git clone https://github.com/ccp4/ccp4i2.git materia-extract
cd materia-extract

# Carve only the paths that move to Materia. Filter-repo preserves
# history of those paths; everything else (including the carve paths'
# ancestral commits that didn't touch them) drops out.
git filter-repo \
  --path apps/compounds/ \
  --path apps/users/ \
  --path Docker/azure-uksouth/.env.demo-materia

# At this point the repo contains only the carved paths' history.
# Push to the new home.
git remote add origin git@github.com:newcastleuniversity/materia.git
git push -u origin main
```

Notes on the invocation:

- `git filter-repo` (not `git filter-branch` — filter-repo is the modern, supported, ~100x faster tool).
- Path-based filtering only; no email rewrites, no message rewrites.
- The fresh clone is critical — filter-repo rewrites all refs, so it must run in a clone you don't care about losing.
- After filter-repo, `git log --all` should show only commits that touched the carved paths. Commits that purely touched (say) `client/renderer/` are dropped from the rewritten history.

## Pre-cut validation checklist

Before pulling the trigger:

- [ ] **Namespaces claimed**: `@ccp4` npm, `ccp4i2-auth` PyPI, `ccp4/ccp4i2-auth` GitHub repo. Without these, Materia can't `npm install` / `pip install` the shared lib post-cut.
- [ ] **`@ccp4/ccp4i2-auth` 1.0.0 published**: a real release on npm + PyPI, not just the workspace artefact. Materia's pre-cut canary build pulls it from registries.
- [ ] **Materia pre-cut canary build green**: in a scratch checkout of `django-sliced`, replace the workspace dependency with the registry version of `@ccp4/ccp4i2-auth`. Run the full compounds-frontend build + the Django test suite. If anything fails, the workspace was hiding a contract leak.
- [ ] **CCP4i2 still builds without the carved paths**: prove `client/renderer/` and `server/` build cleanly with `apps/compounds/` and `apps/users/` deleted from the tree (mock the deletion locally; do not commit). Catches any sneaky runtime imports that survived the inventory.
- [ ] **Service contract draft circulated**: `CCP4I2_SERVICE_CONTRACT.md` v0 reviewed by CCP4 dev team and merged into the eventual `ccp4/ccp4i2-auth` README. Materia's external-consumer status depends on a contract that exists.
- [ ] **DDU image-tag pinning verified**: confirm `ccp4i2-bicep-*` apps are pinned to a specific tag of `ccp4i2/{web,server}`, not a moving reference. They must be insulated from any subsequent CCP4i2 image rebuild that might inadvertently change behaviour.
- [ ] **Materia repo created on `newcastleuniversity`**: empty, with house rules set (LICENSE, branch protection, README sketch). Don't push the rewritten history into a default-config repo.

## Post-cut sequence

1. **Push rewritten history to `newcastleuniversity/materia`**.
2. **Materia first build from its own repo** — clone fresh, install deps from registries, verify build + tests green. *This is the canary.*
3. **CI configured on Materia repo** — GitHub Actions wired up.
4. **Materia first deploy from its own image build** — push `materia/web` and `materia/server` to the shared ACR from Materia's CI, deploy to `materia-demo-*` apps. Verify everything still works end-to-end.
5. **Transition period (carve paths in both places)** — `ccp4/ccp4i2:django-sliced` retains `apps/compounds/` and `apps/users/` for some weeks. Canonical copy is now Materia; sync direction is materia → django-sliced (one-way) for any updates that need to land in CCP4i2's deploy targets pre-deletion. Avoid bilateral edits during this period.
6. **Drop the CCP4i2 overlay** — once Materia's standalone build proves stable in the demo deploy, remove the compounds COPY block from `Docker/client/Dockerfile` and the `Dockerfile.cross-platform` companion. CCP4i2's web image stops carrying compounds.
7. **Delete `apps/compounds/` and `apps/users/` from `ccp4/ccp4i2`** — final step. After this, the only place the compounds work lives is Materia.

## Risks specific to the filter-repo step

| Risk | Mitigation |
|---|---|
| A path is silently shared (compounds imports a renderer module via a deep relative path) and the cut breaks Materia's build | Pre-cut canary build (checklist above) — exercises the import graph against registry versions of the shared lib. The dry run is the test. |
| `filter-repo` rewrites a ref we care about (`django-sliced`, `django`) | Always run filter-repo in a fresh, throwaway clone. Never on a working tree. |
| DDU production breaks because the CCP4i2 web image lost something it needed | Image-lineage separation (DDU pulls from `ccp4i2/*`, never `materia/*`) plus pinned tags mean DDU is decoupled from any post-cut rebuild *until* a deliberate redeploy. The overlay deletion is a deliberate step that can be tested on a demo CCP4i2 deploy first. |
| The transition-period dual location drifts (someone edits `apps/compounds/` in CCP4i2 and forgets to mirror) | Sync direction declared one-way (materia → django-sliced) and the CCP4i2 copy gets a CODEOWNERS / pre-commit guard that warns on edits. Probably overkill for a single-contributor situation; explicit instruction in the plan should suffice. |
| Published `@ccp4/ccp4i2-auth` on npm/PyPI doesn't quite match the workspace version (something in the publish pipeline filters or renames) | First-publish dry run before the cut. Pull the published artifact into a scratch project and verify it imports cleanly. |

## What this plan deliberately does *not* yet decide

- **License choice** for Materia. Open question per the strategy doc.
- **CI provider and concrete pipeline shape** for Materia. Mostly mechanical — settled when the repo lands.
- **First co-maintainer onboarding sequence**. Not gating the cut; happens after Materia is its own repo and there's something for them to maintain.
- **Cross-repo identity-resolution API surface**. Per the proposal, identity primitives are deferred from v0 of the shared lib (LOCKED #2). Revisit once Materia is standing on its own.

## Status

Draft. Not yet executed. Pending the pre-cut validation checklist (most items not yet started). When this is ready to execute, update [REPO_SPLIT_IMPLEMENTATION.md](REPO_SPLIT_IMPLEMENTATION.md) "Currently in flight" and walk through the checklist — every checkbox before the invocation runs.
