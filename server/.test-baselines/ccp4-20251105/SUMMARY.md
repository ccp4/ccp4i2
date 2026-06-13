# i2run baseline — ccp4-20251105

- **Date:** 2026-06-12
- **Commit:** f2e0cf700 (django)
- **Result:** 170 tests · **148 passed · 4 failed · 18 skipped** · 73.1 min
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20251105/bin/ccp4.setup-sh
- **Runner:** server/run_i2run_baseline.sh
- **Machine-readable:** results.xml (JUnit) — diff target for the ccp4-10 run

## Known baseline failures (pre-existing on old CCP4 — NOT migration regressions)

1. `test_modelcraft.py::test_8xfm` — FileNotFoundError: XYZOUT.pdb not produced
2. `test_modelcraft.py::test_gamma_ep` — FileNotFoundError: XYZOUT.pdb not produced
3. `test_nucleofind.py::test_1hr2` — SystemExit: 2
4. `test_nucleofind.py::test_1hr2_raw` — SystemExit: 2

Both tools are ML-dependent (NucleoFind needs downloaded models; ModelCraft
chains nautilus/buccaneer). Likely environment/tooling gaps, not code bugs.

## Spotting regressions after migrating to ccp4-10

A regression = a test **green here but red on ccp4-10**. The 148 passes are the
protected set. To diff:

```bash
# After installing ccp4-10 (e.g. /Users/nmemn/Developer/ccp4-10):
#   edit CCP4_SETUP / CCP4_LABEL in run_i2run_baseline.sh, then:
bash server/run_i2run_baseline.sh
# compare server/.test-baselines/ccp4-20251105/results.xml
#    vs  server/.test-baselines/ccp4-10/results.xml
```
