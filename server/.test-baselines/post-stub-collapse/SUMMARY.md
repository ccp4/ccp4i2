# i2run baseline — post-stub-collapse

- **Date:** 2026-09-01
- **Commit:** d2b936a61 (collapse-stub-classes, PR #329, stacked on #328)
- **Result:** 182 tests · **169 passed · 0 failed · 13 skipped** · exit 0
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh
- **Machine-readable:** results.xml (JUnit)

Covers the merge of all 209 stub classes into their implementations, and the
`ERROR_CODES` fix that made it possible (#328).

Run with `CCP4I2_TEST_CACHE_MAX_MB=400`, which keeps the 335 MB xia2 archive in
the download cache. Without it that fixture re-fetches every run and cost the
`post-stub-inheritance` baseline both xia2 tests to an HTTPS read timeout.

## What this run tests that the snapshot cannot

The container snapshot was byte-identical on all six surfaces, so it says the
*shape* is unchanged. It cannot say whether code still resolves against the
merged classes --- every `isinstance`, every `custom_class="..."` lookup, every
plugin importing `CObsDataFile` or `CPdbDataFile`. Those run in plugin bodies,
and 169 tasks exercise them.

`aimless` and `chltofom` matter most here: they work the MTZ hierarchy that
carried the worst of the MRO bypass, where `isinstance(obs, CMtzDataFile)` was
False and `lib/utils/files/digest.py` imports stub names to work around it.
