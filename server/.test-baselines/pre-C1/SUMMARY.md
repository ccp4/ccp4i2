# i2run baseline — pre-C1 (ccp4-20260702)

- **Date:** 2026-08-25
- **Commit:** 59579c932 (django, 3.1.0a28)
- **Result:** 179 tests · **150 passed · 10 failed · 19 skipped** · 51.3 min
- **Purpose:** the comparator for C1 of
  [docs/error-handling-remediation.md](../../../docs/error-handling-remediation.md).
  Captured before any change, so the effect of C1 is a diff rather than an
  impression.

All ten failures were diagnosed **before** C1 was started, so that any red in
the post-C1 run either matched a known cause or belonged to C1. Full analysis
is in the remediation document under *The pre-C1 baseline*; in brief:

| Failures | Cause |
|---|---|
| `crank2`, `shelx` ×3 | `shelxc` absent from this CCP4 build. SHELX ships separately under Sheldrick's licence (see the ccp4-20260520 summary), so this is an environment gap, not a defect. Two other tests already skip for it. |
| `nucleofind` ×2 | **See the correction below.** |
| `substitute_ligand` ×2 | `f8087b070` (2026-06-15) replaced mmdb2's content-based `ReadCIFASCII` with `gemmi.read_structure`, which detects format by extension. `SubstituteLigand` copies an mmCIF input verbatim into a file it names `selected_atoms.pdb`, so dimple's reader fails. Audited in [docs/coordinate-format-fidelity.md](../../../docs/coordinate-format-fidelity.md). |
| `dm_multidomain` | Test asserts on two mask arrays that stopped being the same shape when `c100c6f53` made masks tight sub-boxes. The masks are disjoint — verified in the unit-cell frame. The test is wrong. |
| `import_merged::test_2ceu_cif` | `64cd3bb3d` (2026-06-15) changed this path's free set from a 20-bin partition to a binary flag. Verified by running the test at `64cd3bb3d^`, where it passes with 21 bins at ~5% each. Needs a decision, not a fix. |

## Correction to the ccp4-20260520 summary

That baseline recorded `test_nucleofind::test_1hr2` / `test_1hr2_raw` as
*"pre-existing env/tooling gap, not a ccp4i2 regression"*, and as "likely a
NucleoFind ML-model/tooling gap". **That is wrong, and the misdiagnosis is why
the task stayed unusable for months.**

`wrappers/nucleofind` declares `TASKCOMMAND` but no `TASKNAME`, and
`CPluginScript.__init__` loads a task's def.xml only `if self.TASKNAME`. The
def.xml was therefore silently never loaded, the container came up empty,
`i2run` built a parser with no task arguments, and argparse rejected `--FPHIIN`
with `SystemExit: 2` — which reads exactly like a tooling problem and is not
one. The `nucleofind` binary is present in this build; the models were never
reached.

A sweep of all 173 registered tasks found `nucleofind` is the only one missing
`TASKNAME`, and also found `SIMBAD` and `pisa`, whose plugin classes do not load
at all. All three are one loop over `TASKS` away from being caught in CI (G1),
and the class of failure — no job, so nowhere to write a diagnostic — is now C8
in the remediation document.

## Note on comparability

Serial (no xdist), as `run_i2run_baseline.sh` requires: per-test results must be
attributable, and the programs under test are themselves parallel. Compare only
against runs made the same way — `post-C1` and `post-C1-fixed` both were.
