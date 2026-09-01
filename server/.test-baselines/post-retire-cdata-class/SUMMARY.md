# i2run baseline — post-retire-cdata-class

- **Date:** 2026-09-01
- **Commit:** 01088a566 (retire-cdata-class)
- **Result:** 182 tests · **169 passed · 0 failed · 13 skipped** · exit 0
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh

Covers the removal of `gui_label` and, with it, the last `@cdata_class` in the
tree. Class-level metadata is now declared one way: an in-class `Meta` block.

Also covers decoupling `content()` declarations from `_metadata` --- they were
read only inside `if metadata:`, which worked because a class with no `Meta`
inherited CData's. Removing CData's decorator exposed it: `CExePath` lost both
its fields.

The eighth consecutive 169/0/13, and the end of the CData simplification
strand. `server/ccp4i2/core` finished the week at 21,012 code lines from
32,913 --- 36% --- with `params.xml`, i2run addressing and validity
byte-identical throughout.
