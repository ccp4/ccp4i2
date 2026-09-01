# i2run baseline — post-hoist-error-codes

- **Date:** 2026-09-01
- **Commit:** e4dff149b (hoist-shared-error-codes, PR #330)
- **Result:** 182 tests · **169 passed · 0 failed · 13 skipped** · exit 0
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh

Covers moving the 47 error codes every CData shares onto CData itself: 1,710
copies and ~5,600 lines removed, verified beforehand as strictly additive
(0 codes lost across 221 classes, 7,565 gained by classes that had never
declared them).

Run with `CCP4I2_TEST_CACHE_MAX_MB=400` so the 335 MB xia2 archive stays
cached.
