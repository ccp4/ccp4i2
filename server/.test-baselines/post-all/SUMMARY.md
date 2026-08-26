# i2run — the first clean run of this sequence (ccp4-20260702)

- **Date:** 2026-08-26
- **Commit:** `101c65a5a` (`tasks-that-never-loaded`, rebased on `django`
  at `8240732c2`, so it carries #284, #285, #286 and #287)
- **Result:** 179 tests · **162 passed · 0 failed · 17 skipped** · 65.8 min
- **Against:** [pre-C1](../pre-C1/SUMMARY.md) — `NEWLY PASSING: 10`,
  `NEWLY FAILING: 0`, and two skips became passes

```
test_crank2                              test_shelx::test_substrdet
test_dm_multidomain                      test_shelx::test_gamma_sad
test_import_merged::test_2ceu_cif        test_shelx::test_gamma_siras
test_nucleofind::test_1hr2               test_substitute_ligand_no_ligand
test_nucleofind::test_1hr2_raw           test_substitute_ligand_with_smiles

skipped -> passed:  test_arcimboldo, test_shelxe_mr::test_gamma
```

## The chain

| Run | Failed | What moved |
|---|---:|---|
| `pre-C1` | 10 | the baseline; all ten diagnosed before anything changed |
| `post-C1` | 16 | C1 revealed six latent defects — kept as evidence |
| `post-C1-fixed` | 10 | those six fixed; identical to the baseline |
| `post-format` | 8 | coordinate format decided by content |
| `post-prefs` | 4 | a job's children see the configured programs |
| `post-registry` | 1 | three tasks that could not load; a mask test |
| `post-all` | **0** | with the freerflag contract restored |

Every run moved only what its own change was responsible for. No test failed
unexpectedly at any point.

## What "zero" depends on

Four of these passes need SHELX, which ships separately under Sheldrick's
licence and is not in this CCP4 build. They pass here because `SHELXDIR` points
at a ccp4-9 installation. On a machine without it, `test_crank2` and
`test_shelx` ×3 will **skip**, and the run header will say why:

```
program resolution (as a job would resolve it):
  shelxc     suite_dir   /Applications/ccp4-9/bin/shelxc
  xds        missing
```

Read that header before comparing this baseline with one from another machine.
A program resolving through a user preference rather than the CCP4 installation
means the run is not testing what a fresh install would do.

## The 17 remaining skips

Not yet audited, and worth it: two of the 19 in `pre-C1` turned out to be
stale, and both now run. Three more look stale on inspection —
`test_find_ligand` skips "until task is in registry" while the task it runs,
`coot_find_ligand`, is registered; `test_coot_script_lines` skips because coot
"may not be available in CI" when it resolves fine here; `test_merge_mtz`
describes a CLI parsing bug that may have been fixed. Three others claim bugs
that the intervening months may have closed (`test_ctruncate`,
`test_editbfac::test_alphafold_pae`, `test_refmac` ×2 — the last blaming a
clipper crash, and clipper has since been factored out).

An unconditional skip records a belief at a moment in time and never re-checks
it. `test_arcimboldo` was skipped "temporarily" in June and passes today.
