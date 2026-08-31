# i2run baseline — post-declaration-order

- **Date:** 2026-08-31
- **Commit:** 1fe2cf6cd (contents-order-from-declaration, PR #326)
- **Result:** 182 tests · **169 passed · 0 failed · 13 skipped** · 68.1 min
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh
- **Runner:** server/run_i2run_baseline.sh
- **Machine-readable:** results.xml (JUnit)

Identical to `post-cdatafile-annotations` (169/0/13), which is the expected
result: the change alters the order children are *rendered* in, and the
snapshot showed `params.xml`, i2run addressing and validity all byte-identical
across 171 tasks. This run confirms that at the level of actually executing
the tasks.
