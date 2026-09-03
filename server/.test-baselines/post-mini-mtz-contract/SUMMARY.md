# i2run baseline — post-mini-mtz-contract

- **Date:** 2026-09-03
- **Commit:** dc406f40a (django, v3.1.0a39) **plus the uncommitted mini-MTZ
  column-contract change in the working tree** — see *What this covers*.
- **Result:** 183 tests · **164 passed · 0 failed · 19 skipped** · exit 0 · 71.5 min
- **CCP4 setup:** /Users/nmemn/Developer/ccp4-20260702/bin/ccp4.setup-sh
- **Diffed against:** post-set-accepts-path (163 / 0 / 19) with
  `scripts/diff_i2run_baselines.py` — newly failing 0, newly passing 0, other
  changes 0, only-in-after 1 (the new test below).
- **Machine-readable:** results.xml (JUnit)

Run with `CCP4I2_TEST_CACHE_MAX_MB=400` (keeps the 335 MB xia2 archive cached).

## What this covers

A mini-MTZ's columns are how every consumer decides what it holds; the
contentFlag on the object is not consulted. `splitMtz` kept source labels for
everything but Phs, so a dropped `beta_blip_P3221.mtz` (`Fobs,Sigma`) became a
file phaser_MR_AUTO's `makeHklin` could not read ("Cannot determine contentFlag
from input file", ObsDataConverter code 2).

In the tree when this ran:

- `CMiniMtzDataFile.canonical_column_mapping()` — the one relabelling, driven by
  `CONTENT_SIGNATURE_LIST`; `splitMtz` and `import_common.canonical_mapping`
  delegate to it.
- `CMiniMtzDataFile.column_contract_violation()` and
  `CPluginScript.checkMiniMtzOutputs()`, run in `postProcess()` between
  `processOutputFiles()` and gleaning: a task that writes a mini-MTZ without
  canonical labels now fails (error 990) instead of persisting a file the next
  job cannot read.

The point of this run is the guard: it fires on *every* producer, so any
wrapper that had been writing a non-canonical mini-MTZ unnoticed would have
turned red here. None did — every outcome matches the previous baseline.

The 19 skips are the same absent-licensed-software set as post-set-accepts-path
(shelx, arp_warp, buster, …) plus the unimplemented xia2_ssx_reduce body.

## Not covered by this run

The phaser `F_OR_I` resolution (`phaser_MR.resolve_f_or_i`, the pipeline-level
call, the validity warning, and the always-visible widget) landed *after* this
run started. It was exercised separately: `test_phaser_expert.py`,
`test_phaser_simple.py`, `test_phaser_rnp.py`, `test_phaser_rnp_pipeline.py`,
`test_phasertng_picard.py` — 13 passed, including the new
`test_beta_blip_amplitudes_with_f_or_i_untouched`.

## New test

- `test_split_mtz::test_auto_detect_beta_blip_relabels_obs_to_canonical` —
  the file that failed in the wild, asserting the split comes out as `F,SIGF`.
