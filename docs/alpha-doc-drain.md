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

## Running it

- **Manually:** Actions → *Alpha-doc feedback drain* → **Run workflow**. Tick
  *dry_run* first to see what it would import without creating anything.
- **From a checkout:** `DRY_RUN=true python3 .github/scripts/alpha_doc_drain.py`
  (needs `gh` authenticated).

## Making it periodic

A GitHub **`schedule:`** trigger only fires from the repository's **default
branch** (currently `main`). This workflow lives on `django` with a manual
trigger only, so it does not yet run on its own. To make it periodic, carry a
copy with a `schedule:` block on the default branch — e.g.

```yaml
on:
  schedule:
    - cron: "17 7 * * *"   # daily, ~07:17 UTC
  workflow_dispatch:
```

That is the only durable option: a scheduled Action needs the default branch,
and there is no off-branch scheduler that survives without it.

[doc]: https://docs.google.com/document/d/1mLbsfvJV0JHdHbwOogGfytc4fWh3M-YKtGq9ulXv81M/edit
