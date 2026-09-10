# Draining the alpha-testing Google Doc into GitHub Issues

The alpha testers' feedback used to live in a Google Doc
([CCP4i2 Django Alpha Testing Feedback][doc]). Now that issues are tracked on
GitHub, this workflow **drains** that doc's *Feedback* section into issues so a
tester who still writes in the doc doesn't get lost.

## What it does

`.github/workflows/alpha-doc-drain.yml` runs `.github/scripts/alpha_doc_drain.py`,
which:

1. Fetches the doc's **HTML export** — `…/export?format=html`. The doc is
   link-shared, so this needs no credentials; the script uses only the Python
   standard library.
2. Parses the `Feedback` `<h1>` section and its per-tester `<h2>` blocks (the
   HTML export keeps the heading levels the plain-text export throws away).
3. Skips the scaffolding blocks (`Anonymous` / `Your Name Here` while they still
   hold their placeholder text) and everything under `Appendix`.
4. Opens one issue per **new** block, labelled `alpha-feedback` +
   `from-tester-doc`, attributing the tester and linking the doc.

### Deduplication

Each issue it creates carries a hidden marker `<!-- alpha-doc-fp:HASH -->`,
where `HASH` fingerprints the block (tester name + body). On every run the
script reads the markers off existing `from-tester-doc` issues and skips those
fingerprints — so re-running is safe and needs **no state file**.

Caveat: the fingerprint covers the block text, so if a tester substantially
edits an already-imported block a fresh issue can appear. That is rare and
low-cost (triage and close the duplicate); the doc is being wound down in favour
of GitHub anyway.

## Running it from a checkout (works today)

```bash
DRY_RUN=true python3 .github/scripts/alpha_doc_drain.py   # report only
python3 .github/scripts/alpha_doc_drain.py                # create issues
```

Needs `gh` authenticated (or `GH_TOKEN` set). This is the whole drain — the
GitHub Actions wrapper below just runs this same script on a schedule.

## Running it as a GitHub Action — needs the default branch

GitHub only registers a workflow that exists on the repository's **default
branch** (`main`). A workflow living *only* on `django` cannot be triggered at
all — not on a schedule, and **not even by `workflow_dispatch`** (the API
returns *"workflow not found on the default branch"*). So there is no way to run
this as an Action without a file on `main`, and CronCreate-style schedulers do
not survive without a live session.

The minimum-footprint way to keep the **logic on `django`** while satisfying
that rule is a small scheduled workflow on `main` that checks out `django` and
runs this script:

```yaml
# .github/workflows/alpha-doc-drain.yml  — on main
name: Alpha-doc feedback drain
on:
  schedule:
    - cron: "17 7 * * *"    # daily, ~07:17 UTC
  workflow_dispatch:
    inputs:
      dry_run: { type: boolean, default: false }
permissions:
  issues: write
jobs:
  drain:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with: { ref: django }          # take the script from django
      - env:
          GH_TOKEN: ${{ github.token }}
          GITHUB_REPOSITORY: ${{ github.repository }}
          DRY_RUN: ${{ inputs.dry_run }}
        run: python3 .github/scripts/alpha_doc_drain.py
```

`main` then carries only this ~15-line scheduler; the parser and issue logic
stay on `django`, maintained there. The `alpha-doc-drain.yml` in *this* branch
is the reference/definition and the `workflow_dispatch` copy for once the
default-branch scheduler exists.

[doc]: https://docs.google.com/document/d/1mLbsfvJV0JHdHbwOogGfytc4fWh3M-YKtGq9ulXv81M/edit
