# i2run baseline — post-cdatafile-annotations

- **Date:** 2026-08-31
- **Commit:** 687d844b0 (annotations-are-the-declaration, PR #325)
- **Result:** 182 tests · **169 passed · 0 failed · 13 skipped** · 68.7 min
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh
- **Runner:** server/run_i2run_baseline.sh
- **Machine-readable:** results.xml (JUnit)

## The four "known baseline failures" are gone

`ccp4-20251105/SUMMARY.md` recorded four failures and judged them
"likely environment/tooling gaps, not code bugs":

1. `test_modelcraft.py::test_8xfm` — FileNotFoundError: XYZOUT.pdb
2. `test_modelcraft.py::test_gamma_ep` — FileNotFoundError: XYZOUT.pdb
3. `test_nucleofind.py::test_1hr2` — SystemExit: 2
4. `test_nucleofind.py::test_1hr2_raw` — SystemExit: 2

All four **pass** on `ccp4-20260702`, with no code change addressing them. That
confirms the judgement: they were gaps in the older CCP4 tree (NucleoFind's
downloaded models, ModelCraft's nautilus/buccaneer chain), not defects.

**So the protected set is now 169, and the expected failure count is zero.**
A single failure in a later run is a regression, with no known-bad list to
discount it against. Earlier notes in this repo describing "4 known failures"
are describing `ccp4-20251105` and are stale for this tree.

## Provenance

Supersedes a discarded run labelled `post-annotations`, which was started
before the `CDataFile` commit and then had the tree changed underneath it
six minutes in — pytest imports from the live filesystem, so part of that
run tested one tree and the rest another. **A worktree with a baseline
running in it is read-only until it finishes**; branch a second worktree
to keep working.
