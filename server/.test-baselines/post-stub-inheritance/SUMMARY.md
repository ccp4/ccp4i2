# i2run baseline — post-stub-inheritance

- **Date:** 2026-08-31
- **Commit:** 1b87061f1 (contents-order-from-declaration, PR #326)
- **Raw result:** 182 tests · 167 passed · 0 failed · 13 skipped · **2 errors** · 63.8 min
- **Effective result:** **169 passed · 0 failed · 13 skipped** — see below
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh
- **Machine-readable:** results.xml (JUnit)

Covers three changes: the declaration-order fallback, `CDataFile.project`
becoming a `CProjectId`, and the removal of 625 lines of inherited
restatements from the stubs.

## The two errors are a download timeout, not a regression

    ERROR test_xia2.py::test_xia2_dials_file      - TimeoutError
    ERROR test_xia2.py::test_xia2_dials_directory - TimeoutError

Both are at **fixture setup**, before any CCP4i2 code runs, and the traceback
is `urllib` → HTTPS → `TimeoutError: The read operation timed out`. The xia2
fixture pulls a 335 MB archive to use twenty frames of it, and files over
100 MB are deliberately **not** cached, so it re-downloads on every run and is
the one fixture exposed to a slow network.

Re-run on the same commit immediately afterwards: **2 passed, 1 skipped in
5:06**. So the protected set of 169 is intact.

## Why this needed checking rather than assuming

167 + 2 accounting for the usual 169 is suggestive, not evidence, and the
static evidence (snapshot byte-identical on all six surfaces) predicts the
same answer — which is exactly the situation in which it is tempting to wave a
failure through. The i2run tier exists to catch what the container snapshot
cannot see, so a run that does not reach two of its tests has not done its
job for those two. Re-running them is cheap; assuming is not.

**If this recurs**, `CCP4I2_TEST_CACHE_MAX_MB=400` brings the xia2 archive
inside the cache and removes the exposure.
