# i2run baseline — glean-fix-clean

- **Date:** 2026-09-22
- **Commit:** branch `fix/glean-status-decoupling`, based on b162055ed (django)
- **Result:** 211 tests · **191 passed · 1 failed · 19 skipped** · 63.3 min
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh
- **Deselected:** `test_moorhen.py` (see below)
- **Machine-readable:** results.xml (JUnit)

Run for the gleaning/terminal-status ordering change: gleaning now precedes the
terminal status write, so a job is FINISHED only once its outputs are in the
database. That sits on the completion path of every task and pipeline, which is
why the whole tier was run rather than a subset.

Compare with `ccp4-20251105` (148 passed / 4 failed / 18 skipped, 2026-06-12):
more tests, and the four ML-dependent failures recorded there are gone.

## The one failure is a suite-context flake, not a regression

`test_xia2.py::test_xia2_dials_file` — `FileNotFoundError` on
`program.xml_tmp -> program.xml` inside `xia2_dials.processOutputFiles()`
(`flushXML`). Attributed by running the 2x2:

| | isolated | in the full suite |
|---|---|---|
| clean b162055ed | pass (145.6 s) | — |
| this branch | pass (136.1 s) | **fail** |

`test_xia2_dials_directory` passed in the same run through the same
`processOutputFiles()` path, and the failure occurs inside `process()`, which
the change never touches. Worth fixing on its own account: `flushXML` writes a
temp file and renames it, and the temp file was missing.

**Method note.** An earlier run of this tier reported `test_dm_vs_parrot.py::
test_parrot_recovery` failing. That was *contamination* — unit suites were being
run concurrently against the same worktree. Both trees pass it isolated and
across an identical 27-file prefix (35 passed / 7 skipped; 390.9 s vs 389.7 s).
Do not run anything else against a tree while its i2run suite is executing.

## test_moorhen.py is deselected

It drives the only `interactive=True` task by faking the window with a
`daemon=True` thread that gives up after 120 s, while the plugin waits on the
session row with no timeout by design. A lost race hangs the whole suite
indefinitely instead of failing (observed: 1.63 s alone, infinite hang 89
job-creating tests deep). The i2run tier should not exercise interactive tools;
the coot_find_ligand / coot_find_waters / coot_rsr_morph tests are headless
wrappers and remain. Its session/dispatch coverage is to be rewritten against
`interactive.drop_file` / `finish_session` directly, CCP4-free.
