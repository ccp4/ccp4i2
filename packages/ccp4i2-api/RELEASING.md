# Releasing `ccp4i2-api` / `@ccp4/ccp4i2-api`

This runbook explains how to cut a release of the bilingual API library
that lives at `ccp4/ccp4i2:packages/ccp4i2-api/`. It is the operational
counterpart to [`README.md`](README.md) (which describes *what* the
package is). Read this before bumping the version, before publishing,
or when something goes wrong during a release.

The intended reader is a maintainer (human or LLM agent) who needs to
either:

1. ship a new version, or
2. understand the existing publication state and what level of
   automation underpins it.

---

## What gets published

One package, two language artifacts, two registries, one workflow.

| Side | Source | Registry artifact | Consumers |
|---|---|---|---|
| Python | `ccp4i2_api/` | [`ccp4i2-api` on PyPI](https://pypi.org/project/ccp4i2-api/) | CCP4i2 Django server, [Materia](https://github.com/martinemnoble1/materia) Django server, [Reinspect](https://github.com/martinemnoble1/pandda-inspect-api) Django server |
| TypeScript | `src/` (built to `lib/`) | [`@ccp4/ccp4i2-api` on npm](https://www.npmjs.com/package/@ccp4/ccp4i2-api) | CCP4i2 Electron client, Materia NextJS frontend, Reinspect Moorhen client |

Both halves implement the same wire contract (auth, api-fetch helpers,
typed payloads). They are versioned in **lockstep**: every release ships
both, with the same `x.y.z`. Skipping a version on one side (e.g. 0.3.2
on PyPI when the npm publish failed) creates a temporary gap; the next
release should re-align.

## Version lineage (current state at 2026-06-05)

| Version | PyPI | npm | Notes |
|---|---|---|---|
| 0.1.0 – 0.2.x | published as `ccp4i2-auth` | published as `@ccp4/ccp4i2-auth` | Old name. Yanked / unpublished after the rename at 0.3.0. |
| 0.3.0 | ✅ | ✅ | First release under the `ccp4i2-api` / `@ccp4/ccp4i2-api` name. Published manually 2026-05-02. |
| 0.3.1 | ✅ | ✅ | Relaxed Django bound from `<5.0` to `<6.0`. Published manually 2026-06-03. |
| 0.3.2 | — | — | Tag exists; never published. Skipped after the publish workflow landed. |
| 0.3.3 | ✅ | — | First *automated* release. npm publish failed due to npm CLI 10.8 vs OIDC requirement (≥11.5). |
| 0.3.4 | ✅ | ✅ | First *successful* end-to-end automated release. |

## Versioning policy

Semver, with the understanding that `0.x.y` indicates the public
contract has not been formally promised stable yet:

| Change | Bump | Example |
|---|---|---|
| Constraint relaxation, metadata-only fixes, doc-only changes | **patch** | `0.3.0` → `0.3.1` (Django bound widened) |
| New auth middleware, new exported type, new helper that doesn't change existing surface | **minor** | `0.3.x` → `0.4.0` |
| Breaking change in auth handshake, removed/renamed exports, request/response shape change | **major** (still inside 0.x — but document loudly) | `0.x.y` → `0.(x+1).0` with the breaking change called out |

If you're not sure, default to **patch** and explain in the commit
message why the change is not load-bearing for existing consumers. A
gratuitously-bumped minor is recoverable; an under-bumped patch that
breaks a consumer is not.

## Branch model

The package source lives in the `ccp4/ccp4i2` monorepo, which has
multiple long-lived branches. Releases follow this convention:

- **`django`** — canonical development and release branch. **All
  release commits land here first.** All version tags point at commits
  on this branch.
- **`moorhen-scenes`** — the branch Materia's submodule pin currently
  tracks. After landing a release commit on `django`, **cherry-pick it
  onto `moorhen-scenes`** so Materia picks the change up at the next
  submodule bump.
- **`main`** — legacy; not used for releases of this package.

The cherry-pick step exists because Materia tracks `moorhen-scenes` for
YAML scene-file support that hasn't yet merged back. When that merges,
Materia will move its pin back to `django` and the cherry-pick step
disappears.

## Tag pattern

```
ccp4i2-api-v<semver>
```

Examples: `ccp4i2-api-v0.3.4`, `ccp4i2-api-v0.4.0`.

The `ccp4i2-api-` prefix disambiguates from the **application** tags
(`v2.6.x`) that live in the same monorepo. The publish workflow
([`.github/workflows/publish-ccp4i2-api.yml`](../../.github/workflows/publish-ccp4i2-api.yml))
matches on this prefix; tags without it do nothing.

Tags are always created on **`django` HEAD**, never on
`moorhen-scenes`. The cherry-picked commit on `moorhen-scenes` has a
different SHA from the tagged commit; that is expected and not a
problem (Materia's submodule pin points at the cherry-pick).

---

## Cutting a release

The full recipe, assuming you're starting from a clean ccp4i2 checkout
with both branches up to date and the change you want to ship not yet
committed.

### 0. Decide the version

See "Versioning policy" above. Don't reuse a version that has been
published or tagged — both `pip` and `npm` refuse to overwrite, and
`git` will reject re-pushing a tag.

### 1. Make the change on `django`

```bash
git checkout django
# … edit code under packages/ccp4i2-api/ …
```

Bump the version in **three files**:

- `packages/ccp4i2-api/pyproject.toml` (`version = "x.y.z"`)
- `packages/ccp4i2-api/package.json` (`"version": "x.y.z"`)
- `packages/ccp4i2-api/package-lock.json` — regenerate via
  `(cd packages/ccp4i2-api && rm -rf node_modules && npm install)`.

Do **not** touch `ccp4i2_api/__init__.py` — `__version__` is derived
from package metadata via `importlib.metadata.version("ccp4i2-api")`,
so the literal can never drift again. (It did, between 0.3.1 and 0.3.2;
the fix landed at 0.3.2 and the runtime version has matched
`pyproject.toml` since.)

Run the tests locally as a sanity check:

```bash
cd packages/ccp4i2-api
python3 -m venv /tmp/release-check
/tmp/release-check/bin/pip install -e ".[test]"
/tmp/release-check/bin/python -m pytest tests/python -q
```

### 2. Commit + cherry-pick

```bash
git add packages/ccp4i2-api/{pyproject.toml,package.json,package-lock.json}
git commit -m "ccp4i2-api x.y.z: <one-line summary>"
DJANGO_SHA=$(git rev-parse HEAD)

git checkout moorhen-scenes
git cherry-pick "$DJANGO_SHA"
```

The cherry-pick should be clean (no conflicts), since the only files
involved are the three version files which rarely diverge between
branches. If there *is* a conflict, that means `moorhen-scenes` has
drifted in those files — resolve and document why.

### 3. Push both branches

```bash
git push origin django
git push origin moorhen-scenes
```

This triggers the test matrix (Python 3.9/3.11/3.12 × Django 4.2/5.2,
minus the impossible 3.9+5.2 cell) on both branches. The publish jobs
are **skipped** here because they gate on tag refs. Wait for both runs
to go green at https://github.com/ccp4/ccp4i2/actions/workflows/publish-ccp4i2-api.yml
before tagging — the tag will re-run the matrix anyway, but if it's
red on a plain branch push you've saved yourself a wasted tag.

### 4. Tag + push tag (on `django`)

```bash
git checkout django
git tag -a ccp4i2-api-vx.y.z -m "ccp4i2-api x.y.z: <summary>"
git push origin ccp4i2-api-vx.y.z
```

This is the moment the publish happens. The workflow re-runs the test
matrix and, on green, builds + publishes to both registries via OIDC.

### 5. Watch the workflow

```bash
gh run watch $(gh run list --workflow=publish-ccp4i2-api.yml --limit 1 --json databaseId --jq '.[0].databaseId') --exit-status
```

A successful run shows:

```
✅ Test (py3.9  / django4.2)
✅ Test (py3.11 / django4.2)
✅ Test (py3.11 / django5.2)
✅ Test (py3.12 / django4.2)
✅ Test (py3.12 / django5.2)
✅ Publish to PyPI
✅ Publish to npm
```

Verify the registries:

```bash
curl -s https://pypi.org/pypi/ccp4i2-api/json | python3 -c "import sys, json; print(json.load(sys.stdin)['info']['version'])"
npm view @ccp4/ccp4i2-api version
```

Both should print `x.y.z`.

### 6. Coordinate downstream

- **Materia**: bump the CCP4i2 submodule pin to pick up the
  cherry-picked release commit on `moorhen-scenes`. Materia's
  Dockerfile pin (`ccp4i2-api>=0.3`) accepts the new version
  automatically; no Dockerfile edit needed. From the Materia repo:
  ```bash
  cd external/ccp4i2 && git fetch && git checkout moorhen-scenes && git pull
  cd ../.. && git add external/ccp4i2 && git commit -m "Bump CCP4i2 submodule for ccp4i2-api x.y.z"
  ```
- **Reinspect**: if the release contains something Reinspect needs
  (e.g. new auth middleware), bump the requirement in Reinspect's
  `requirements.txt`. Otherwise no action — Reinspect's `>=0.3.1,<0.4`
  (or similar) range will pick it up at next install.
- **CCP4i2 itself**: ccp4i2-api is consumed in-tree via
  `pip install -e packages/ccp4i2-api/`; no version bump needed.

---

## What happens behind the scenes

The publish workflow at [`.github/workflows/publish-ccp4i2-api.yml`](../../.github/workflows/publish-ccp4i2-api.yml)
defines three jobs:

| Job | When it runs | What it does |
|---|---|---|
| `test` | every push/PR touching `packages/ccp4i2-api/` or the workflow itself; every tag push | 5-cell pytest matrix |
| `publish-pypi` | only on tag push matching `ccp4i2-api-v*` (or `workflow_dispatch` from such a tag) | `python -m build` + `pypa/gh-action-pypi-publish` |
| `publish-npm` | same trigger as PyPI | `npm install -g npm@latest` + `npm ci` + `npm run build` + `npm publish --access public --provenance` |

The two publish jobs use OIDC — no `PYPI_TOKEN` or `NPM_TOKEN` is
stored anywhere. Each registry's trusted-publisher configuration
points at the workflow file and an environment name (`pypi` / `npm`);
the OIDC token issued by GitHub Actions at job start is exchanged for
a short-lived registry credential.

### One-time trusted-publisher configuration (already done; record only)

If the trust relationship is ever lost (e.g. registry-side
de-configuration) and needs re-establishing:

#### PyPI

1. https://pypi.org/manage/project/ccp4i2-api/settings/publishing/
2. "Add a new publisher" with:
   ```
   Owner:               ccp4
   Repository name:     ccp4i2
   Workflow filename:   publish-ccp4i2-api.yml
   Environment name:    pypi
   ```

#### npm

1. https://www.npmjs.com/package/@ccp4/ccp4i2-api/access → "Trusted Publishers" tab
2. "Add publisher" → "GitHub Actions" with:
   ```
   Organization or user: ccp4
   Repository:           ccp4i2
   Workflow filename:    publish-ccp4i2-api.yml      (must end with .yml)
   Environment name:     npm
   Allow npm publish:    ✓
   Allow npm stage publish: leave unchecked
   ```

#### GitHub environments (optional belt-and-braces)

The trusted-publisher configs above reference environment names
(`pypi`, `npm`). Creating GitHub environments with those names under
https://github.com/ccp4/ccp4i2/settings/environments lets you add
required reviewers — gating publishes behind a one-click approve in
the GitHub UI. Without the environments existing, the workflow still
runs and publishes (the OIDC token carries the environment name from
the workflow file regardless), but a typo or mistake on the npm side
will not be caught by a human eye until after the publish happens.

The PyPI/npm trusted-publisher records are the load-bearing piece;
GitHub environments are reviewer-gating only.

---

## Troubleshooting

### "npm error code E404" / "'@ccp4/ccp4i2-api@x.y.z' is not in this registry"

**This bit us at 0.3.3.** The error message is misleading; it does not
mean the package or scope is missing. It means the registry PUT request
arrived without a valid OIDC token, so the registry treated the call
as unauthenticated and returned 404 (its standard not-authorised
response for write attempts).

**Root cause:** `actions/setup-node@v4` with Node 20 ships
**npm 10.8.x**. That version supports OIDC for sigstore provenance
signing but **not** for the registry PUT itself. OIDC-based registry
auth (the feature that lets you skip `NPM_TOKEN`) landed in
**npm 11.5**.

**Fix (already in the workflow):** a `npm install -g npm@latest` step
between `setup-node` and `npm ci`. Adds ~5s; brings npm up to a
version that knows how to use OIDC for the registry call.

**How to diagnose:** the failed run's `Publish to npm` job will show:

```
npm notice publish Signed provenance statement with source and build
            information from GitHub Actions
npm notice publish Provenance statement published to transparency log
npm error code E404
npm error 404 Not Found - PUT https://registry.npmjs.org/@ccp4%2fccp4i2-api
```

The fact that provenance was signed proves OIDC works for sigstore;
the 404 on the next call proves the registry isn't seeing a valid token.

### Workflow doesn't trigger on tag push

Confirm the tag name matches `ccp4i2-api-v*` exactly. A typo like
`ccp4i2-api-0.3.4` (no `v`) will be ignored by the workflow. The tag
itself still pushes; it just doesn't fire anything. Delete and re-tag:

```bash
git tag -d ccp4i2-api-x.y.z
git push origin :refs/tags/ccp4i2-api-x.y.z
# re-tag with the correct pattern
```

### Test matrix fails on a cell that locally passed

Most likely a Django version that's not in your local venv. The CI
matrix pins Django via `django~=${{ matrix.django }}.0` whereas a
local `pip install -e ".[test]"` resolves whatever floor the
`pyproject.toml` accepts. If your edit relies on a Django feature
newer than 4.2, the 4.2 cells will fail; choose the right floor in
`pyproject.toml` (and document the deprecation) rather than papering
over the failure.

### `npm publish` reports "version already exists"

Both registries refuse to re-publish a version. If you need to fix a
bad release, bump to the next patch and publish again — never try to
re-publish the same version. (Yanking a release on PyPI / unpublishing
on npm is a separate, last-resort operation; talk to the package owner
before attempting either.)

### `gh run watch` says the run succeeded but the registry doesn't show the new version

Wait 30–60 seconds. PyPI and npm both have read-replica propagation
delays that can leave `curl /pypi/.../json` or `npm view` stale for a
short window after the publish has actually landed. If you still don't
see it after a minute, click into the workflow run on github.com and
verify the publish job log shows the "Publishing to ..." success
message, not a swallowed error.

---

## Quick reference — release checklist

```
[ ] Decide version (semver; not previously published)
[ ] Bump 3 files on `django`: pyproject.toml, package.json, package-lock.json
[ ] Local test: pytest tests/python -q (against current Django)
[ ] git commit on django
[ ] git checkout moorhen-scenes && git cherry-pick <django-sha>
[ ] git push origin django moorhen-scenes
[ ] Verify branch test runs green at github.com/ccp4/ccp4i2/actions
[ ] git checkout django && git tag -a ccp4i2-api-vx.y.z -m "..."
[ ] git push origin ccp4i2-api-vx.y.z
[ ] gh run watch <id> --exit-status
[ ] Verify both registries report x.y.z (curl PyPI, npm view)
[ ] (Materia) bump submodule pin if needed
[ ] (Reinspect) bump requirement if release contains something needed
```
