# i2run baseline — ccp4-20260702 (July 2026 CCP4 10 build, macOS arm64)

- **Date:** 2026-09-26 (run after the ccp4-20260904 baseline, same commit)
- **Commit:** 82b1a45fa (django, release 3.1.0a79) plus the Coot
  molecules-container fix that became #621
- **Result:** 224 tests · **203 passed · 0 failed · 21 skipped** · 59.9 min
- **Purpose:** the like-for-like control for `ccp4-20260904`: same code, same
  day, same machine. Results and skip set are identical to that label's; see
  its SUMMARY.md for the skip table and the build differences.
- Confirms the Coot fix in #621 is safe on this build too (it has both
  container names and the helper resolves to `molecules_container_t`).
