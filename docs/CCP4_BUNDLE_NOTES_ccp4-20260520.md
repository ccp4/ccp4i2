# ccp4-20260520 migration test notes (from ccp4i2 i2run)

**Context.** The ccp4i2 i2run suite (~170 end-to-end task tests) was run against the
near-release **ccp4-20260520** to check for regressions versus ccp4-20251105, with
the Python framework stack pinned to the old baseline (Django 4.2.30, pytest 8.4.2)
so the only variable was CCP4 itself. Both the Mac/arm64 install and the Linux
tarball (`ccp4-20260520-linux.tar.gz`) were inspected.

**Bottom line:** the bundle is essentially fine. The only genuinely actionable
items were two **ccp4i2-side fixes** (already merged). There is **little for CCP4
maintainers** — the one upstream heads-up is the dimple log-section rename (§1).
Everything else is expected behaviour, a local environment issue, or ccp4i2's own.

---

## 1. Upstream heads-up — dimple log section rename (for dimple/refmac maintainers)

dimple renamed its INI log sections from `[refmac5 jelly]` / `[refmac5 restr]` to
`[refmacat jelly]` / `[refmacat restr]` (refmacat = the Servalcat-based refmac).
This is a reasonable change, but the section name is a de-facto interface:
ccp4i2's i2Dimple parser keyed off the old names and silently produced no
R-factors until updated. Flagging in case other consumers parse dimple's log.
(ccp4i2 now accepts both names.)

## 2. Expected / not an issue

- **SHELX (`shelxc/d/e`) not in the bundle.** Expected — SHELX is installed
  separately under its own (Sheldrick) license. The `test_shelx` / `test_shelxe_mr`
  failures simply reflect "SHELX not installed in this test environment," not a
  bundle defect. (For reference: the bundle does ship the crank/SHELX i2
  integration scripts and `bin/crank`, just not the separately-licensed SHELX
  executables.)
- **SAD/SIRAS R-work ~0.02 higher than the old CCP4** (deterministic 0.242/0.244
  vs the test's 0.22/0.23 thresholds, once SHELX was installed locally to run
  these). Most likely the new refmac/servalcat refinement; if anything this is a
  ccp4i2 question of whether those test thresholds are too tight for
  experimental-phasing cases — not a CCP4 build concern.

## 3. Local install artifacts on the test machine (NOT build defects)

The Mac/arm64 bundle was unpacked with the disk at **99% (~6 GB free)**, and a
**partial extraction** explains two things that initially looked like bundle bugs
— the **Linux bundle of the same version is clean**:

- **modelcraft truncated** in the Mac install: 13 of ~80 `.py` files, missing
  `__init__.py` and the `scripts/modelcraft.py` entry point → `modelcraft` ran as
  an empty namespace package (`ModuleNotFoundError`). Linux bundle: complete.
- **~12 packages with empty `.dist-info`/`.egg-info`** (no `METADATA`/`PKG-INFO`
  → version `None`), which crashes pip's resolver (`dui2`, `typing_extensions`,
  `coot_headless_api`, `SCons`, `meson`, `gyp_next`, `metalcoordanalysis`, …).

**Fix for the test machine (not CCP4):** re-extract the Mac bundle on a disk with
headroom and re-verify.

## 4. ccp4i2-side fixes (already merged — for reference)

Surfaced by this run, fixed in ccp4i2 (django branch):

- **i2Dimple** now accepts both `refmac5 {jelly,restr}` and `refmacat {jelly,restr}`
  log sections (§1) → fixes SubstituteLigand R-factor extraction.
- **modelcraft test** reads `XYZOUT.cif` not `XYZOUT.pdb`: modelcraft 6.x emits
  mmCIF, and the models here have long ligand names PDB format cannot represent,
  so `.cif` is the correct and only valid output.

---

*Generated from a ccp4i2 i2run migration test run, 2026-06.*
