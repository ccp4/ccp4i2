# i2run baseline — post-docs-stack (ccp4-20260904, macOS arm64)

- **Date:** 2026-10-01
- **Commit:** ec4a9d8d2: django 639e11ad2 (#705) with the user-help stack
  (#706–#709) and #710 (fileOut= inside a list item) merged in. #711 was
  added afterwards; see *The one failure*.
- **Result:** 243 tests · **227 passed · 1 failed · 15 skipped** · 70.5 min
- **Against `ccp4-20260904`** (82b1a45fa, 2026-09-25/26: 203 passed, 0
  failed, 21 skipped): no test that passed there fails here because of the
  new code.
- **Environment difference:** `SHELXDIR=/Applications/ccp4-9/bin` was set,
  so the SHELX tests ran (the reference skipped them).

## Test by test against ccp4-20260904

| Change | Tests |
|---|---|
| passed → failed | `test_failure_surfacing::test_a_failure_in_a_subjob_names_the_subjob` |
| skipped → passed | `test_shelx` (3), `test_shelxe_mr` (SHELXDIR set); `test_arcimboldo`, `test_crank2` |
| new, passed | 19: Aimless dataset names, EDSTATS resolution default, ProSMART reference-model import, MakeLink (4), PanDDA dispatch (6), parrot ASU copies, servalcat Platonyzer, SubstituteLigand report, Phaser EP inverted hand, import_merged without a free set, fileOut= inside a list item |

The 15 skips are the reference's remainder: xds_par not installed, the
stub tests with no body, and the opt-in PanDDA end-to-end tests.

## The one failure

`test_a_failure_in_a_subjob_names_the_subjob` fails on django itself, not
because of the new code. `git bisect` over the 56 commits since the
reference names #638. It made UNSATISFACTORY a completion in
`recordCauses` and, with it, downgraded the errors such a job inherits from
its subjobs. `aimless_pipe` turns a failed step into UNSATISFACTORY, so its
pointless step's error reached the pipeline as a warning and the pipeline
"failed and reported nothing". Fixed by #711 (downgrade only when the job
succeeded); with #711 merged in, `test_failure_surfacing.py` passes 3/3.

## How this was run

In a separate worktree of the combined branch, from `server/`:
`SHELXDIR=/Applications/ccp4-9/bin CCP4_SETUP=<ccp4-20260904>/bin/ccp4.setup-sh
CCP4_LABEL=post-docs-stack bash run_i2run_baseline.sh`, with
`DJANGO_SETTINGS_MODULE` unset. pytest imported the worktree's own
`ccp4i2`, not the copy installed in the CCP4 tree.
