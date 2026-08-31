# i2run baseline — post-error-code-inheritance

- **Date:** 2026-08-31
- **Commit:** 146e176e2 (stubs-inherit-error-codes)
- **Result:** 182 tests · **169 passed · 0 failed · 13 skipped** · 65.8 min
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh
- **Machine-readable:** results.xml (JUnit)

Covers six changes: the declaration-order fallback, `CDataFile.project` as a
`CProjectId`, removal of the inherited restatements, removal of `attributes=`,
error codes inheriting along the MRO, and the deletion of
`qualifiers_definition`.

Run with `CCP4I2_TEST_CACHE_MAX_MB=400`, which brings the 335 MB xia2 archive
inside the download cache. The previous baseline lost both xia2 tests to an
HTTPS read timeout while re-fetching it; with the archive cached there was no
network exposure and both passed. Worth keeping for any run on a slow link.
