# GitHub Actions

What runs, when, and what it would take to change it. Four workflows, in
`.github/workflows/`.

| Workflow | Fires on | Does | Cost |
|---|---|---|---|
| `ci.yml` | every pull request; pushes to `django` | Runs `tests/unit/` + `tests/api/unit/` on stock CPython 3.11 and 3.13 | ~1 min × 2 |
| `electron-multiplatform-build.yml` | pushes to `django`; PRs into `django` — both `paths-ignore` docs | Packages the desktop app on mac, Windows and Linux. Compiles, does not test | ~10 min × 3 |
| `publish-ccp4i2-api.yml` | changes under `packages/ccp4i2-api/`; tags `ccp4i2-api-v*` | Tests the shared package (Python 3.9/3.11/3.12 × Django 4.2/5.2); on a tag also publishes to PyPI and npm | minutes |
| `release.yml` | tags `v*` | Verifies the tag matches `ccp4i2.__version__`, publishes the backend wheel to PyPI, builds the three installers, creates the GitHub Release | ~20 min |

Nothing else is automated. There is no linting, type-checking, coverage or
scheduled job.

## What is and isn't gated

`ci.yml` is new (2026-08-24). Before it, **no workflow ran a server or client
test** — `release.yml` only checked that the tag matched the version, and the
Electron build proves the app compiles and packages, not that it works. The
`ccp4i2-api` package has always had its own test matrix, but that covers the
shared auth/API library, not CCP4i2 itself.

So the honest position today:

- **Server Python** — gated by `ci.yml` on every PR.
- **Client TypeScript** — not gated. `tsc --noEmit` is not run anywhere; the
  Electron build would fail on a hard compile error but is skipped entirely for
  docs-only changes.
- **A release** — not gated. `release.yml` publishes whatever the tag points at.
  A red `django` can still be cut.

Making `ci.yml` a **required check** is a repository setting, not something in
the workflow file: Settings → Branches → branch protection for `django` →
*Require status checks to pass*, selecting both matrix jobs. Until that is set
the job reports but does not block.

## `ci.yml` in detail

It runs the two suites that need nothing from the CCP4 suite — no binaries, no
`$CCP4` tree, no downloaded reference data — because `server/pyproject.toml`
declares every dependency they need. That is ~800 tests in about a minute.

Three decisions in it are load-bearing and easy to undo by accident:

**The two suites run as separate `pytest` invocations.** The conftest under
`tests/api/` builds a Django test database; once one exists, the report-fixture
tests in `tests/unit/lib/` take a database path and fail looking for jobs
belonging to another suite's fixtures. Each suite is green alone. Combining them
into one command reintroduces eight failures that have nothing to do with the
code under test.

**The job asserts CCP4 is absent.** If a runner image ever grew a CCP4 install,
the job would keep passing while quietly proving something weaker. It fails
instead.

**Anything needing an absent dependency skips, rather than failing.** The
`servalcat`/`ctruncate` converter tests, the `libtbx.phil` modules and one
`scipy` test are guarded, so the same command is green with or without CCP4 —
run it under `ccp4-python` locally and you get strictly more coverage, not a
different answer. Adding a test that hard-requires a CCP4 binary without a guard
turns the job red for everyone.

What is deliberately *not* in it: `tests/db/` and `tests/async/` do not collect
without more setup, and `tests/i2run/` and `tests/api/e2e/` run real
crystallographic jobs and download from PDBe/RCSB. Those need a machine with
CCP4; see [Testing](../CLAUDE.md) for running them locally.

## Maintaining these

**A `pull_request` run uses the workflow file from the PR's own head, not from
the base branch.** So a change to a workflow takes effect on the PR that makes
it — convenient — but a PR branched *before* that change still runs the old
version. This bites when stacking: PR B based on PR A does not see A's workflow
edits until B is rebased. (Widening `ci.yml`'s trigger and then wondering why
the stacked PR had no checks is exactly how this entry got written.)

**`paths-ignore` is evaluated over the whole changeset, not per file.** A commit
mixing a doc edit with a source edit still triggers the Electron build. Do not
rely on it to skip work for a mixed PR.

**Action versions drift.** The tree currently pins `actions/checkout@v4`,
`actions/setup-python@v5`, `actions/setup-node@v4`,
`actions/upload-artifact@v4`, `actions/download-artifact@v4`,
`softprops/action-gh-release@v2` and `pypa/gh-action-pypi-publish@release/v1`.
Release runs currently log a Node 20 deprecation warning, forced onto Node 24 by
the runner; bumping the `actions/*` majors clears it. Keep versions consistent
across workflows so one upgrade is one decision.

**Publishing uses OIDC, not stored tokens.** Both PyPI publishes and the npm
publish authenticate via GitHub's identity token against a Trusted Publisher
configuration, so there is no secret to rotate — but the trust is configured
per project on PyPI, and it is keyed to the *workflow filename*. Renaming
`release.yml` or `publish-ccp4i2-api.yml` breaks publishing until the
corresponding PyPI setting is updated. See [RELEASING.md](RELEASING.md).

**Testing a workflow change** without merging it: push the branch and open a PR
(the head's copy runs, per above), or use `workflow_dispatch` where the workflow
declares it — `ci.yml`, `release.yml` and `publish-ccp4i2-api.yml` all do.

## See also

- [RELEASING.md](RELEASING.md) — cutting a release, `scripts/cut-alpha.sh`, the
  one-time PyPI Trusted Publisher setup, and troubleshooting a failed run.
- [macos-signing-setup.md](macos-signing-setup.md) — the signing and
  notarisation secrets the desktop build consumes.
- `packages/ccp4i2-api/RELEASING.md` — releasing the shared package.
