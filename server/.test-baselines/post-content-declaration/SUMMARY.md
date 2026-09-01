# i2run baseline — post-content-declaration

- **Date:** 2026-09-01
- **Commit:** 38626df19 (content-declaration)
- **Result:** 182 tests · **169 passed · 0 failed · 13 skipped** · exit 0
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh

Covers `content()` and the conversion of 307 fields across 81 classes, on top
of the stub collapse and the error-code hoist.

The fifth consecutive 169/0/13 for this strand. Over the week `server/ccp4i2/core`
went from 32,913 code lines to 23,075 --- 29% --- with `params.xml`, i2run
addressing and validity byte-identical at every step.
