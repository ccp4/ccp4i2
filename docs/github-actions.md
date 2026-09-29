# GitHub Actions

What runs, when, and what it would take to change it. Four workflows, in
`.github/workflows/`.

| Workflow | Fires on | Does | Cost |
|---|---|---|---|
| `ci.yml` | every pull request; pushes to `django` | Runs `tests/unit/` + `tests/api/unit/` on stock CPython 3.11 and 3.13 | ~1 min × 2 |
| `electron-multiplatform-build.yml` | pushes to `django`; PRs into `django` — a `changes` job skips the three builds for docs-only changesets | Packages the desktop app on mac, Windows and Linux. Compiles, does not test | ~10 min × 3 |
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

What is deliberately *not* in it: `tests/i2run/` and `tests/api/e2e/` run real
crystallographic jobs and download from PDBe/RCSB. Those need a machine with
CCP4; see [Testing](../CLAUDE.md) for running them locally.

## Maintaining these

**A `pull_request` run uses the workflow file from the PR's own head, not from
the base branch.** So a change to a workflow takes effect on the PR that makes
it — convenient — but a PR branched *before* that change still runs the old
version. This bites when stacking: PR B based on PR A does not see A's workflow
edits until B is rebased. (Widening `ci.yml`'s trigger and then wondering why
the stacked PR had no checks is exactly how this entry got written.)

**Jobs are skipped by a job-level `if:`, never by a path filter, and that is
deliberate.** `build-mac`, `build-win` and `build-linux` are required checks on
`django`. A workflow skipped by `paths-ignore` never reports them: they sit at
"Expected" for ever, so until 2026-09-18 a Markdown-only PR could not be merged
at all (the first project skill had to ride on a code PR to land). A *job*
skipped by an `if:` reports success. So the workflows always start, decide what
the changeset needs, and gate each job on that.

`.github/scripts/changes-need.sh` is the one rule, shared by both workflows. It
prints three decisions:

```
backend=true|false     the Python unit + API suites      (ci.yml)
frontend=true|false    the client's vitest + tsc         (ci.yml)
build=true|false       the mac/win/linux packages        (electron-multiplatform-build.yml)
```

- **Docs-only** — every changed file is `*.md`, under `docs/`, `LICENSE`, or
  under `server/.test-baselines/` (committed i2run results and summaries:
  evidence about a build, never an input to one) — needs nothing.
- **`server/`** needs the backend suites only. The desktop builds run `npm` in
  `client/` and `packages/` and read nothing under `server/`, so they cannot say
  anything about a backend change; on one measured server-only PR they were 87%
  of the billable minutes and none of the signal (macOS bills at 10x, Windows at
  2x). And the client's vitest tests use a mocked fetch, so they never reach the
  server either.
- **`client/`**, **`packages/`** or a root lockfile needs the frontend tests and
  the desktop builds. `packages/` is the shared API contract, so it needs the
  backend suites too.
- **`.github/`** needs everything: those files decide what runs, so a mistake in
  them is invisible to a narrowed run.
- A mixed changeset needs the union, and anything undecidable — an unknown base,
  a failing git, a path no rule claims — needs everything.

A pull request is judged against its merge base, so a base branch that has moved
on with code does not make a docs PR build. The script is plain bash and can be
run locally:
`.github/scripts/changes-need.sh <base-sha> <head-sha> [merge-base|direct]`.

It is tempting to make a change to `server/ccp4i2/__init__.py` force a desktop
build, because the app pins the exact backend version it requires. Don't: a
version bump always touches the client pin as well (`scripts/cut-alpha.sh` edits
both files and greps to confirm), so the `client/` rule already claims it, and
the invariant is checked directly by the **Version lockstep** job below.

**Never skip a MATRIX job at job level.** A skipped job publishes the
unexpanded template as its check name — `Unit + API tests (Python ${{
matrix.python-version }})` — because there is no matrix to expand for a job that
never ran. The required contexts `(Python 3.11)` and `(Python 3.13)` then never
report, sit at "Expected", and the pull request cannot be merged. That is the
same breakage `paths-ignore` used to cause, arriving by a different route: it was
introduced by gating `unit-tests` at job level and a docs-only PR hit it
immediately (#659). So `unit-tests` runs always and skips every *step*, which
keeps the names right and costs seconds. `frontend-tests` and the three builds
are not matrix jobs, so their names are static and they are safely gated at job
level.

Three ways to break the skipping, all tempting: putting `paths-ignore` back,
gating a matrix job as above, and adding a second workflow that reports the same
job names for docs changes. The last is worse than it looks: a PR touching docs
and code runs both, and a passing stub can mask a failing build.

**The version lockstep is checked on every pull request, in seconds.** An alpha
app and its backend are strictly bound: the app pins the EXACT version it
requires and refuses any other, so `CCP4I2_REQUIRED_SERVER_VERSION` in
`client/main/ccp4i2-server-version.ts` must equal `ccp4i2.__version__`. Nothing
checked that — `release.yml`'s `verify-version` compares the *tag* with
`__version__` and the `ccp4i2-api` lock with the wheel floor, but never the
client pin, and `cut-alpha.sh` only holds the invariant for bumps that go
through it. A hand-edit to either file would have shipped an app that refuses
its own backend. `.github/scripts/check-version-lockstep.sh` runs unconditionally
in `ci.yml` and again in `release.yml`; it takes no arguments and can be run
locally.

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
