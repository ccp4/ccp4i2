# i2run confirmation — post-C1 with fixes (ccp4-20260702)

- **Date:** 2026-08-25
- **Commit:** 8e68bdcbc (`c1-process-output-files`)
- **Result:** 179 tests · **150 passed · 10 failed · 19 skipped** · 49.3 min
- **Against:** [pre-C1](../pre-C1/SUMMARY.md) — `NEWLY FAILING: 0`,
  `NEWLY PASSING: 0`, `OTHER CHANGES: 0`

Identical to the baseline, test for test. C1 is in, the six defects it revealed
are fixed, and no test changed outcome in either direction. The ten remaining
failures are the pre-existing set diagnosed in the pre-C1 summary; none is C1's.

No `xfail` markers were needed. The strict-xfail machinery in the remediation
document exists for revealed defects too large to fix in the same pass, and none
of the six was.

Unit (860) and api/unit (79) pass under `ccp4-python` and under the CCP4-free
venv — the latter matters because two of the five fixes are in core.

**This is the comparator from here.** The next work to measure against it is the
coordinate-format audit, where `substitute_ligand` ×2 are expected to go green.

Compare with `python3 scripts/diff_i2run_baselines.py pre-C1 post-C1-fixed`.
