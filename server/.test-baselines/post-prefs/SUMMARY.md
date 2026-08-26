# i2run run — tests and jobs agree about "installed" (ccp4-20260702)

- **Date:** 2026-08-26
- **Commit:** d659bd212 (`tests-honour-program-preferences`)
- **Result:** 179 tests · **158 passed · 4 failed · 17 skipped** · 61.2 min
- **Against:** [post-format](../post-format/SUMMARY.md) — `NEWLY PASSING: 4`,
  `NEWLY FAILING: 0`, two skips became passes

```
NEWLY PASSING: test_crank2, test_shelx::{substrdet, gamma_sad, gamma_siras}
OTHER CHANGES: test_arcimboldo   skipped -> passed
               test_shelxe_mr    skipped -> passed
```

Six tests now run that did not, and the skip count falls from 19 to 17.
`test_arcimboldo` had been marked `skip` unconditionally since June; a probe
self-heals where a hand-written skip never does.

The four SHELX/crank2 failures were recorded in the pre-C1 summary as an
environment gap. They were a defect — see the correction there. SHELX resolves
on this machine through the `SHELXDIR` preference; the jobs could not find it
because a task's children inherited an unmodified `PATH`.

**Header stamp.** Every run now records how each gated program resolved:

```
program resolution (as a job would resolve it):
  shelxc     suite_dir   /Applications/ccp4-9/bin/shelxc
  ...
  xds        missing
```

so two baselines that disagree say why in their own headers. Read it when
comparing runs from different machines: a program resolving through a *user
preference* rather than the CCP4 installation means the run is not testing what
a fresh install would do.

## The four that remain

| Failures | Status |
|---|---|
| `nucleofind` ×2 | No `TASKNAME`, so the def.xml has never been loaded. One line, plus the G1 loop that also catches `SIMBAD` and `pisa`. C8 in the remediation document. |
| `dm_multidomain` | The test ANDs two mask arrays that stopped being the same shape when masks became tight sub-boxes. The masks are disjoint. Fix the test. |
| `import_merged::test_2ceu_cif` | The free set changed from a 20-bin partition to a binary flag in `64cd3bb3d`. The test was detecting a real divergence from CCP4 freerflag; fixed on `freerflag-conformance`. |
