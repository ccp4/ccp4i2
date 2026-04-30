# Request: please create `ccp4/ccp4i2-auth` on GitHub

*A focused ask. The full proposal lives in [`apps/compounds/docs/CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md`](https://github.com/newcastleuniversity/materia/blob/main/apps/compounds/docs/CCP4I2_FACING_PROPOSAL_REPO_SPLIT.md) (in the Materia repo post-cut); this doc is the one-page version that the GitHub repo creation gates on.*

## The ask (one sentence)

**Create an empty repository `ccp4/ccp4i2-auth` on the `ccp4` GitHub org, with [@martinemnoble](https://github.com/martinemnoble) as repo admin / maintainer. Five-minute task; everything downstream is on us.**

## What goes there

A bilingual auth + API-contract library, **already shipping**:

| Artefact | Where it is now |
|---|---|
| npm package | [`@ccp4/ccp4i2-auth@0.1.0`](https://www.npmjs.com/package/@ccp4/ccp4i2-auth) (TypeScript) |
| PyPI package | [`ccp4i2-auth==0.1.0`](https://pypi.org/project/ccp4i2-auth/) (Python; Django middleware) |
| Source | [`ccp4/ccp4i2:django-sliced @ packages/ccp4i2-auth/`](https://github.com/ccp4/ccp4i2/tree/django-sliced/packages/ccp4i2-auth) (in the monorepo, awaiting extraction) |
| Tests | 31 Python tests + 13 TypeScript tests, all green |

After the GitHub repo exists, Martin extracts it from the monorepo via `git filter-repo` (preserving full history) and updates `ccp4/ccp4i2` to import the lib from the registry instead of the workspace. The library gets its own release cycle and stops cluttering `ccp4i2`.

## Why `ccp4` org (not Newcastle)

Three independent codebases consume this library:

1. **CCP4i2 desktop** — Electron app's auth uses it (LocalSession provider).
2. **CCP4i2 cloud / Materia** — Azure AD auth uses it (MSAL provider, AzureAD middleware).
3. **`i2remote` and any future external integrator** — bearer-token API calls use it.

It's structurally a CCP4 community asset, not a Newcastle product. The npm scope `@ccp4` and the PyPI name `ccp4i2-auth` were claimed from that reasoning (April 2026, Martin direct, ex-chair authority). The GitHub repo home should match.

If it sat under `newcastleuniversity/`, every external integrator would have to look outside the CCP4 org for a CCP4 contract — that's exactly the wrong signal.

## What this commits CCP4 to

- **Hosting the repo.** That's the load-bearing commitment.
- **Saying no to PRs that would break the contract.** Anyone (including Martin) can be asked to "wait, that's a major-version-bump-required change". The repo is the place where that conversation happens.
- **Co-ownership on the registries.** As CCP4 admins (Paul, Stuart, Dave) come online with npm and PyPI accounts, they get added as co-owners. Provenance binding (`@ccp4` npm ↔ `ccp4` GitHub) follows once the repo exists.

## What this does NOT commit CCP4 to

- Maintaining Materia or any compounds-side work — that lives at [`newcastleuniversity/materia`](https://github.com/newcastleuniversity/materia), under Newcastle's roof, not yours.
- Day-to-day maintenance of the auth library — Martin owns that. CCP4 admins are reviewers / safety-net, not authors.
- Learning MSAL / Azure AD — the provider abstraction hides those details. CCP4i2 desktop developers see `LocalSessionTokenProvider`; they don't need to understand the cloud side.
- Cloud-deploy responsibilities, Azure infrastructure, or any of the deployment-shaped concerns.
- Renaming, restructuring, or versioning of CCP4i2 itself — this is purely about giving the auth library a clean home.

## What's already in flight, gated on this

1. **`@ccp4` npm ↔ `ccp4` GitHub provenance binding** — one button-click on each side, but needs the GitHub repo to exist first.
2. **CCP4i2 monorepo cleanup** — `packages/ccp4i2-auth/` lives in `ccp4/ccp4i2:django-sliced` today as a "pending extraction" workspace; nothing changes for CCP4i2 contributors today, but every PR review that touches the auth library has to deliberate "is this a contract change?". Cleaner once it lives elsewhere.
3. **Service contract circulation** — the v0 contract document at [`docs/CCP4I2_SERVICE_CONTRACT.md`](CCP4I2_SERVICE_CONTRACT.md) describes what CCP4i2 promises through the library's TS types. ~37 endpoints, conservative, walked-through. Awaits dev-team review; the library's GitHub home is the natural place for that review to land.

## Concrete action

```
Repo name:   ccp4/ccp4i2-auth
Visibility:  public
Initial admin/maintainer: @martinemnoble
Initial content: empty (Martin pushes the filter-repo'd history)
License: LGPL-3.0-or-later (already in the published packages; matches CCP4 conventions)
```

Once the empty repo exists, Martin will:

1. `git filter-repo` extract `packages/ccp4i2-auth/` from `ccp4/ccp4i2` preserving full history.
2. Push to `ccp4/ccp4i2-auth:main`.
3. Update `ccp4/ccp4i2:django-sliced` to import the lib from npm + PyPI instead of the workspace path.
4. Open a PR (or push directly if you'd rather) on the new repo to wire up GitHub Actions CI.
5. Bind `@ccp4` npm org to the `ccp4` GitHub org for provenance attestation.

Steps 1–3 are mechanical. Step 4 is small. Step 5 is two settings clicks across two web UIs.

## What we need from you

- Create the repo.
- Add @martinemnoble as a maintainer / admin.
- (Optional, when ready) accept the npm + PyPI co-owner invitations as your accounts come online.

That's it.
