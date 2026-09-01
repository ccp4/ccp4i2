# i2run baseline — post-meta-migration

- **Date:** 2026-09-01
- **Commit:** ca84582a9 (meta-migration)
- **Result:** 182 tests · **169 passed · 0 failed · 13 skipped** · exit 0
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh

Covers moving class-level metadata into an in-class `Meta` block for 219
classes across 16 modules, leaving one `@cdata_class` in the tree --- on
`CData` itself, which cannot use the hook it defines.

The seventh consecutive 169/0/13. Over the week `server/ccp4i2/core` went from
32,913 code lines to 21,012 --- 36% --- with `params.xml`, i2run addressing and
validity byte-identical at every step.
