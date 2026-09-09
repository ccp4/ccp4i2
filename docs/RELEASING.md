# Releasing CCP4i2

> For the workflows themselves — what runs on a PR, what is gated and what
> is not, and what to watch when editing them — see
> [GitHub Actions](github-actions.md).

One version tag produces a complete release: the **`ccp4i2` backend wheel on
PyPI** and the **desktop installers attached to a GitHub Release**. Driven by
[`.github/workflows/release.yml`](../.github/workflows/release.yml).

> This is the app/backend release. The shared `ccp4i2-api` library has its own,
> separate flow — see [`packages/ccp4i2-api/RELEASING.md`](../packages/ccp4i2-api/RELEASING.md).

## Single source of truth

The version lives in [`server/ccp4i2/__init__.py`](../server/ccp4i2/__init__.py)
(`MAJOR`/`MINOR`/`PATCH` + `PRERELEASE`). Everything else derives from it:

- the PyPI wheel version (setuptools `attr: ccp4i2.__version__`);
- the **exact** backend version the Electron app pins to (stamped into the build
  from the tag);
- the GitHub Release name.

The workflow **asserts the tag equals `ccp4i2.__version__`** and fails fast on a
mismatch, so the desktop app and PyPI can never drift apart at release time.

## Cutting a release

> **`django` is protected by a ruleset:** every change must go through a **pull
> request** with all required checks green, and **merge commits are disabled**
> (squash/rebase only). A direct `git push … django` is rejected with `GH006`.
> So a release is **two steps** — bump-via-PR, then tag-after-merge — not one
> push. `scripts/cut-alpha.sh` drives both.

**Step 1 — open the release PR.** `cut-alpha.sh` (no args) does the pre-flight:
syncs `django`, bumps `PRERELEASE` + the desktop exact-pin default
(`client/main/ccp4i2-server-version.ts`) in lockstep, updates the date, runs the
ccp4i2-api lock≥floor guard, refuses a version already on PyPI or an existing
tag, then commits on a `release-vX` branch, pushes it, and opens the PR into
`django`. It does **not** tag yet.

```bash
git checkout django
scripts/cut-alpha.sh              # aN -> a(N+1): bump, branch, push, open PR
scripts/cut-alpha.sh --dry-run    # show the plan, change nothing
scripts/cut-alpha.sh --version 3.1.0b1   # explicit version (e.g. move to beta)
scripts/cut-alpha.sh --no-push    # commit on the release branch locally only
```

**Review, wait for all checks green, and squash-merge the PR.**

**Step 2 — tag the merged bump.** Once the PR is on `django`:

```bash
scripts/cut-alpha.sh --tag        # tag the merged commit on django, push the tag
scripts/cut-alpha.sh --tag --dry-run   # show what it would tag
```

`--tag` reads the version back **from `django`** (so the tag can never disagree
with what merged), refuses if the tip isn't a fresh `release: …` commit or the
tag/PyPI version already exists, then pushes only the tag — which is what fires
the Release workflow.

<details>
<summary>Manual equivalent (if not using the script)</summary>

```bash
# Step 1 — bump on a release branch and PR it (NOT a direct push to django).
git checkout django && git pull
git checkout -b release-v3.1.0a1
#   edit server/ccp4i2/__init__.py -> PRERELEASE (+ date), and the exact-pin
#   default in client/main/ccp4i2-server-version.ts to match.
git commit -am "release: ccp4i2 3.1.0a1"
git push ccp4 release-v3.1.0a1
gh pr create --base django --head release-v3.1.0a1 --title "release: ccp4i2 3.1.0a1"
#   ... review, all checks green, squash-merge ...

# Step 2 — after merge, tag the merged commit on django and push the tag.
git fetch ccp4
git tag -a v3.1.0a1 ccp4/django -m "CCP4i2 3.1.0a1"
git push ccp4 v3.1.0a1          # <-- this fires the Release workflow
```

A **stacked** PR (one branched off another still-open PR) does not auto-rebase
when its base is squash-merged: `git rebase --onto ccp4/django <old-base>
<branch>` and force-push to drop the now-duplicated squashed commit.
</details>

The workflow then:

1. **verify-version** — tag `v3.1.0a1` must match `ccp4i2.__version__` (3.1.0a1);
   a pre-release version (`a`/`b`/`rc`) is flagged so the Release is marked
   **Pre-release**.
2. **publish-pypi** — builds `server/` sdist+wheel and publishes `ccp4i2` to
   PyPI via OIDC.
3. **build-desktop** — builds macOS `.dmg`, Windows `.exe`, Linux `.AppImage`
   with `CCP4I2_SERVER_VERSION_FLOOR=3.1.0a1` stamped in (the env var name is
   historical; it now carries the **exact** required version).
4. **release** — creates the GitHub Release and attaches all three installers as
   public, permanent, no-login assets.

A tag like `v3.0.2-rc1` is published as a **pre-release** (GitHub marks it as
such automatically; PyPI treats `rcN` as a pre-release version).

## One-time setup (before the first `v*` tag)

PyPI Trusted Publishers are **per project**. The existing `ccp4i2-api` publisher
does **not** cover `ccp4i2`. Add one for the main package:

- <https://pypi.org/manage/project/ccp4i2/settings/publishing/> →
  - **Owner:** `ccp4`
  - **Repository:** `ccp4i2`
  - **Workflow:** `release.yml`
  - **Environment:** `pypi` (reuses the existing GitHub environment — no new
    secret, no token)

No `PYPI_TOKEN` is required or used — publishing is OIDC trusted publishing, the
same as `ccp4i2-api`.

## macOS signing & notarisation (optional but recommended)

The release **degrades gracefully**: with no Apple secrets configured, the macOS
build ships **unsigned** (users clear quarantine with `xattr -cr`, exactly as
today). Add the secrets below and the Release workflow **signs and notarises**
the `.dmg` automatically — no `xattr` needed by users, no "app is damaged"
scare.

You need **two** Apple credentials — the API key notarises but does **not**
sign:

1. **Developer ID Application certificate** (`.p12`) — for code signing.
2. **App Store Connect API key** (`.p8` + Key ID + Issuer ID) — for notarising.

Add these as repository (or org) **secrets**:

| Secret | What it is |
|---|---|
| `MAC_CSC_LINK` | base64 of the Developer ID Application `.p12` (`base64 -i cert.p12 \| pbcopy`) |
| `MAC_CSC_KEY_PASSWORD` | the `.p12` export password |
| `APPLE_API_KEY_P8` | the **contents** of the App Store Connect `AuthKey_XXXX.p8` |
| `APPLE_API_KEY_ID` | the API Key ID |
| `APPLE_API_ISSUER` | the API Issuer ID |

The `build-desktop` job detects these: when `MAC_CSC_LINK` **and**
`APPLE_API_KEY_P8` are present it signs (`hardenedRuntime` + entitlements from
[`client/assets/entitlements.mac.plist`](../client/assets/entitlements.mac.plist))
and runs `electron-builder … -c.mac.notarize=true`; otherwise it builds unsigned
with `CSC_IDENTITY_AUTO_DISCOVERY=false`.

> **Validate with a pre-release tag first**, e.g. `v3.0.2-rc1` — signing +
> notarisation can only be truly confirmed by building on CI and opening the
> resulting `.dmg` on a clean Mac (no `xattr`). Once it's green, drop the `xattr`
> note from `docs/give-it-a-try.md`.

## Notes & gotchas

- **A PyPI version can be published only once.** If a release fails after the
  wheel is uploaded, bump the pre-release counter (`a1`→`a2`) and cut a new tag —
  you cannot re-upload the same version.
- **Pre-release versioning + exact pin (alpha discipline).** The version carries
  a PEP 440 pre-release segment (`PRERELEASE = "a1"` in `server/ccp4i2/__init__.py`).
  Pre-releases are invisible to a plain `pip install ccp4i2`, and the desktop app
  pins the backend **exactly** (`ccp4i2==<version>`, see
  [`client/main/ccp4i2-server-version.ts`](../client/main/ccp4i2-server-version.ts)),
  so an alpha app and backend are strictly bound and can't mix with the escaped
  `3.0.x` line. A pre-release tag (`vX.Y.Za1`) is auto-marked **Pre-release** on
  GitHub. To declare a stable line, set `PRERELEASE = ""`, tag `vX.Y.Z`, and (if
  desired) relax the app back to a floor/compatible-range check.
- **Optionally yank the superseded finals.** `3.0.0`/`3.0.1` can be *yanked* on
  PyPI (PEP 592) so normal resolution won't pick them, without deleting them
  (exact pins still work). Alpha apps are unaffected (they pin exact anyway).
- **Backend build needs no CCP4.** `server/ccp4i2/__init__.py` imports only the
  standard library, so the wheel builds on a bare runner.

## Troubleshooting a failed release run

- **Run fails, all jobs "cancelled", none ran a step (~15 min, no logs).**
  GitHub couldn't assign hosted runners — the run message says *"The job was not
  acquired by Runner of type hosted even after multiple attempts."* This is a
  transient GitHub Actions capacity/quota blip, **not your code**. Nothing
  publishes (jobs never run), so just **re-run the same run** — no new version:
  ```bash
  gh run rerun <run-id> --repo ccp4/ccp4i2 --failed
  ```
  Find the id with `gh run list --repo ccp4/ccp4i2 --workflow release.yml --branch vX.Y.ZaN --limit 1`.
  Only if the wheel already published (verify with
  `curl -s -o /dev/null -w '%{http_code}' https://pypi.org/pypi/ccp4i2/<version>/json`
  → 200) must you bump the version instead, because PyPI is immutable.

- **verify-version fails: "runtime lock ccp4i2-api==X < wheel floor >=Y".** The
  bundled lock (`server/ccp4i2/requirements-runtime.txt`) pins an older
  `ccp4i2-api` than the wheel requires. The desktop app installs the **lock**
  version (via `--no-deps`), so this would ship the wrong `ccp4i2-api`. Update the
  lock's `ccp4i2-api==` pin (regenerate with `server/scripts/gen_runtime_lock.py`
  if other deps changed too) and re-cut. `cut-alpha.sh` runs this same guard
  before tagging, so it catches the drift locally first.

- **macOS build fails: "MAC verification failed during PKCS12 import (wrong
  password?)".** `MAC_CSC_KEY_PASSWORD` doesn't match the `.p12` in
  `MAC_CSC_LINK`. Verify the pair **locally** before re-setting the secrets — see
  [macos-signing-setup.md](macos-signing-setup.md). Signing is opt-in
  (repo variable `ENABLE_MAC_SIGNING`); leave it off to ship unsigned and skip
  this class of failure entirely.

- **`gh workflow run` can't find a workflow.** `workflow_dispatch` only registers
  for workflows present on the **default branch** (`main`). All alpha work is on
  `django`, so branch-dispatch of a django-only workflow won't appear.

- **A build fails with `Failed to FinalizeArtifact: (403) Forbidden` after the
  installer already uploaded.** This is the org's **Actions artifact storage**
  at (or over) quota, not a code fault — a single re-run may squeak through, but
  the real fix is to prune. Each CI/release run uploads ~1.8 GB (three
  installers), so it fills fast. To measure: `gh api
  /repos/ccp4/ccp4i2/actions/artifacts --jq '.total_count'` and sum
  `size_in_bytes` of the non-`expired` ones (use `gh run list`, **not**
  `/actions/artifacts?per_page=100`, which times out expanding `workflow_run`).
  To prune safely: enumerate the artifact-producing workflows (`Electron
  Multiplatform Build`, `CI`) with `gh run list --json databaseId,headBranch`,
  keep the newest few runs, and `gh api -X DELETE
  /repos/ccp4/ccp4i2/actions/artifacts/{id}` the rest. These workflows only ever
  ran on 3.x branches, so pruning them cannot touch the 2.x release assets (those
  are **GitHub Release** assets — separate storage — not Actions artifacts).

- **"The desktop build shipped the wrong value" — verify the real packaged
  build, don't trust an installed app.** `next.config.ts` gates settings on
  `BUILD_TARGET === "electron"`, baked at `next build` time. To confirm what a
  release actually shipped, inspect the artifact itself, not `/Applications`:
  `gh release download vX --pattern '*arm64.dmg'`, `hdiutil attach` it, and
  `grep` the value in `…/ccp4i2-django.app/Contents/Resources/app.asar`. Installed
  apps and mounted dmgs are routinely **stale** — always check
  `CFBundleShortVersionString` (`/usr/libexec/PlistBuddy -c 'Print
  :CFBundleShortVersionString' …/Info.plist`) before trusting that a build is the
  version you think it is. A tester's "still broken in aN" is often a pre-fix
  build, but — as the 100 MB import cap showed — it can also be a real bug the
  packaged build genuinely carries; the dmg is the arbiter.
