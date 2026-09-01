# i2run baseline — post-retire-dead-args

- **Date:** 2026-09-01
- **Commit:** cbb3d2225 (retire-dead-decorator-args)
- **Result:** 182 tests · **169 passed · 0 failed · 13 skipped** · exit 0
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh

Covers the removal of `qualifiers_order` (1,537 lines, stored and never read),
154 empty `error_codes={}`, and the conversion of the last two classes using
`attributes=`.

The sixth consecutive 169/0/13 for this strand.
