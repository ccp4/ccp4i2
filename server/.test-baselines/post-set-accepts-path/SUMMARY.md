# i2run baseline — post-set-accepts-path

- **Date:** 2026-09-02
- **Commit:** 849962d25 (cdatafile-set-accepts-cstring, PRs #341 + #342)
- **Result:** 182 tests · **163 passed · 0 failed · 19 skipped** · exit 0 · 56.3 min
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh

Covers `CFilePath.__fspath__` (#341) and `CDataFile.set()` accepting a path
that is not a `str` (#342), on a tree verified as 0 commits behind `django`.

## 163/19, where earlier baselines were 169/13

Same 182 tests. Six moved from passed to skipped, none to failed, and every
one names a program that is genuinely absent:

    test_shelx.py       not installed: shelxc, shelxd
    test_arcimboldo.py  not installed: shelxe
    test_arpwarp.py     ARP/wARP not installed ($warpbin unset)
    test_xia2.py:58     not installed: xds_par

That is the optional-program resolution working: separately-licensed software
that is not present now skips with a reason naming where it was looked for,
instead of running and failing. Not attributable to either PR.
