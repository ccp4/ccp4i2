# i2run — a stricter CLI, and four skips lifted (ccp4-20260702)

- **Date:** 2026-08-26
- **Commit:** `fd27a43f9` (`honest-skips`)
- **Result:** 179 tests · **166 passed · 0 failed · 13 skipped** · 63 min
- **Against:** [post-all](../post-all/SUMMARY.md) — `NEWLY FAILING: 0`,
  `NEWLY PASSING: 0`, four skips became passes

```
skipped -> passed:  test_find_ligand::test_8xfm
                    test_merge_mtz::test_gamma_merge_two_mtz
                    test_refmac::test_8xfm
                    test_refmac::test_gamma
```

## Why this run was taken

The change refuses a `key=value` sub-field the target object does not declare,
where the parser previously accepted it and created a junk attribute. That can
only *reveal* callers who were relying on the silence, so it needed the same
treatment as C1: measure against a known-clean baseline, and read any red as
this change's doing.

Nothing went red. The paths that had looked risky were checked first and each
turns out to be handled before the guarded branch or to use fields that really
are declared:

| Form | Path | Outcome |
|---|---|---|
| `fullPath=`, `columnLabels=`, `annotation=`, `file=`, `seqFile=` | special-syntax handlers, earlier | never reach the check |
| `--DOMAINS segments= mode=`, `--UNMERGEDFILES file=`, `--MINIMTZINLIST fileName=`, `--ENSEMBLES use=` | the guarded branch | all declared |
| `--ASSEMBLY m=A` | guarded branch, item declares nothing | check skips by design |
| `imageFile/baseName=`, `selection/text=`, `columnList[0]/columnLabel=` | nested-path branch | not guarded |

## What the run cannot say

The suite exercises the parser only where a test uses it. A user's script
passing a mistyped keyword to a task no test drives will now stop — correctly,
since the job it used to run had that input unset — but this run is not evidence
about that. It is a **breaking change for such scripts**, and belongs in the
release notes rather than being met in the wild.

## The four skips

Three were stale (see the pre-C1 summary's audit): `test_find_ligand` waited on
a task registration that had already happened, and the `test_refmac` pair blamed
a clipper crash from before clipper was factored out. `test_merge_mtz` is the
one this branch actually fixed, four layers deep.

Skips across the sequence: 19 (pre-C1) → 17 (post-prefs) → 13 here. Six of the
thirteen remaining are tests with no body, and now say so.
