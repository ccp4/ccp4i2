# i2run run — coordinate format fidelity (ccp4-20260702)

- **Date:** 2026-08-25
- **Commit:** 40d00a9fd (`coordinate-format-fidelity`)
- **Result:** 179 tests · **152 passed · 8 failed · 19 skipped** · 49.7 min
- **Against:** [post-C1-fixed](../post-C1-fixed/SUMMARY.md) —
  `NEWLY PASSING: 2`, `NEWLY FAILING: 0`, `OTHER CHANGES: 0`

The two that went green are `test_substitute_ligand_no_ligand` and
`test_substitute_ligand_with_smiles`, red since 2026-06-15. The core change
touches every task that loads a coordinate file, so a full run rather than the
eight targeted tests is what makes "nothing else moved" mean anything.

Work described in
[docs/coordinate-format-fidelity.md](../../../docs/coordinate-format-fidelity.md):
format decided by content rather than by filename, and each of the two selection
APIs made to keep the promise its name makes, on the no-selection path as well
as the selection path.

Unit (887) and api/unit (79) pass under `ccp4-python` and under the CCP4-free
venv.

## The eight that remain

| Failures | Status |
|---|---|
| `crank2`, `shelx` ×3 | `shelxc` absent from this CCP4 build. SHELX ships separately under Sheldrick's licence. Wants one shared `skipif` on a binary probe, replacing four hand-written skips elsewhere and four hard failures here. |
| `nucleofind` ×2 | `wrappers/nucleofind` declares no `TASKNAME`, so its def.xml has never been loaded. A one-line fix, and the class of defect is C8 in the remediation document. |
| `dm_multidomain` | The test ANDs two mask arrays that stopped being the same shape when masks became tight sub-boxes. The masks are disjoint, verified in the unit-cell frame. Fix the test. |
| `import_merged::test_2ceu_cif` | `64cd3bb3d` changed this path's free set from a 20-bin partition to a binary flag. Needs a decision on which is intended, not a fix. |

Three of the four are small. The fourth is a question for a person.
