# i2run run — post-C1, before the fixes (ccp4-20260702)

- **Date:** 2026-08-25
- **Commit:** 0eb1287da (`c1-process-output-files`)
- **Result:** 179 tests · **144 passed · 16 failed · 19 skipped** · 50.5 min
- **Against:** [pre-C1](../pre-C1/SUMMARY.md) — `NEWLY FAILING: 6`,
  `NEWLY PASSING: 0`, `OTHER CHANGES: 0`

Kept as the evidence of what C1 revealed. This is not a regression: a test that
flipped here was green while the wrapper it ran had already concluded its own
output was unusable, and said so into a `print()`.

**All six newly-failing carried error code 993** — a Python exception raised
inside `processOutputFiles()`, previously caught and printed to the server
console. None was a wrapper's deliberate `return FAILED`, which is why four of
the six fell outside the predicted at-risk list: `--predict-red-list` models the
`FAILED` return and cannot model exceptions.

| Test | Exception | Consequence, every run, for months |
|---|---|---|
| `test_chainsaw` | `TypeError: float() argument must be … not 'CInt'` | no sequence identity |
| `test_cmapcoeff` | `TypeError: tostring() got an unexpected keyword argument 'pretty_print'` | no `program.xml` written |
| `test_phaser_ensembler` | `AttributeError: 'CPdbData' object has no attribute 'mmdbManager'` | remarked ensemble never written |
| `test_sculptor` | `AttributeError: … no attribute 'setFromPdbDataFile'` | no atom-count KPI |
| `test_xia2` ×2 | `FileNotFoundError: … cmtzsplit.log` | no xia2 log file copied out of `LogFiles/` |

Five distinct causes; `chainsaw` had two. All fixed on the same branch, one
commit each. Triage was a single `grep` for the code, with no reading of
tracebacks to classify — which is the argument for reserving a distinct code per
core-detected failure mode.

The other ten failures are the diagnosed pre-existing set, unchanged.
