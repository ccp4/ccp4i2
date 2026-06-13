# i2run candidate run — ccp4-20260520 (near-release ccp4-10)

- **Date:** 2026-06-12
- **Commit:** f2e0cf700 (django)
- **Result:** 170 tests · **143 passed · 8 failed · 19 skipped** · 65.1 min
- **Baseline:** ccp4-20251105 = 148 pass / 4 fail / 18 skip
- **Stack pinned to baseline** (CCP4-only variable): Django 4.2.30, pytest 8.4.2,
  pytest-django 4.11.1, pytest-asyncio 1.3.0, xdist 3.8.0, django-filter 24.3.
  (Unavoidable diffs: Python 3.11.11→3.11.14 patch, gemmi → 0.7.5, both part of
  the CCP4 distribution itself.)

## Diff vs baseline

**Regressions (green on old → red on new): 4 — in two clusters**

### Cluster 1 — SHELX/crank pipelines (NOT a ccp4i2 bug)
- test_shelx::test_gamma_sad   (best R 0.2429, threshold 0.22)
- test_shelx::test_gamma_siras (FileNotFoundError: n_REFMAC5.pdb — pipeline broke)
- test_shelxe_mr::test_gamma   (pass → **skip** "SHELXE not installed")

Root cause: **SHELX binaries (shelxc/d/e) are absent from the ccp4-20260520
bundle** — confirmed: `bin/shelx*` all missing.

**Follow-up (binaries restored):** copied shelxc/d/e from /Applications/ccp4-9/bin
into the candidate bin/. SHELX is **closed-source x86_64-only** (no arm64 build),
so on this Apple-Silicon machine it runs via **Rosetta 2**. After restoring:
- test_shelxe_mr::test_gamma -> PASS
- test_shelx::test_substrdet  -> PASS
- test_shelx::test_gamma_sad / test_gamma_siras -> still fail, but now on a
  *quality threshold* (R 0.242/0.244 vs 0.22/0.23), no longer missing-binary.
  Likely SHELXD stochastic variance near a tight threshold; variance check run.

**Strategic risk for CCP4 (bigger than the tests):** a native-arm64 bundle that
drops the x86_64 shelx breaks SHELX/crank phasing on Apple Silicon, and the
Rosetta 2 fallback is on Apple's deprecation path (available through macOS 26/27,
limited after). The question for CCP4 is the arm64 SHELX story, not just
re-bundling. No ccp4i2 change needed for this cluster.

### Cluster 2 — SubstituteLigand / DIMPLE (genuine ccp4i2 regression)
- test_substitute_ligand::test_substitute_ligand_no_ligand   (IndexError)
- test_substitute_ligand::test_substitute_ligand_with_smiles (IndexError)

`IndexError: list index out of range` on `assert rworks[-1] < 0.23`, because
`rworks = xml.iter("r_factor")` is **empty**. The pipeline actually SUCCEEDS —
all file outputs (XYZOUT.pdb, MTZs, dict) are produced, and dimple refines fine
(`dimple.log: R/Rfree 0.2145/0.2376`, which would PASS the thresholds). The
break is purely in **result extraction**: the new CCP4's dimple/refmac no longer
emits `<r_factor>`/`<r_free>` into the program.xml that ccp4i2's i2Dimple scraper
aggregates. So this is a ccp4i2-vs-new-dimple **XML-schema compatibility** issue,
not a science regression.
- Root cause (found): newer dimple renamed its INI log sections from
  `[refmac5 jelly]`/`[refmac5 restr]` to `[refmacat jelly]`/`[refmacat restr]`
  (refmacat = the Servalcat-based refmac). i2Dimple's hardcoded section-name
  checks missed them, so no `<r_factor>` was written. The R-factor data is all
  present under the new section names.
- **FIX APPLIED:** wrappers/i2Dimple/script/i2Dimple.py now matches both
  `refmac5 {jelly,restr}` and `refmacat {jelly,restr}` sections (in log order),
  with a `has_option` guard. Backward-compatible. Verification re-run in progress.

**Fixed (red→green): 0**

**Still red on both (pre-existing, unchanged): 4**
- test_modelcraft::test_8xfm, test_modelcraft::test_gamma_ep (ML/tooling)
- test_nucleofind::test_1hr2, test_1hr2_raw (ML model / SystemExit 2)

## Net actionable for the migration
- **CCP4 build side:** re-bundle SHELX (shelxc/d/e) — also fix the 9 packages that
  shipped empty metadata (no Version:) which crash pip's resolver.
- **ccp4i2 side:** one real regression — the i2Dimple R-factor scraper vs the new
  dimple/refmac XML. Two tests, one root cause.

---

## FINAL STATUS (after fixes/repairs)

ccp4i2 code/test fixes (real, verified green):
- **i2Dimple.py** — accept `refmacat {jelly,restr}` sections → substitute_ligand x2 PASS.
- **test_modelcraft.py** — read `XYZOUT.cif` (modelcraft 6.x emits mmCIF; PDB can't
  hold the long ligand names) instead of `XYZOUT.pdb` → modelcraft x2 PASS.

CCP4 bundle repairs (env-side, NOT ccp4i2 — for the CCP4 build team):
1. ~12 packages shipped empty metadata (no Version:) → crashes pip. (synthesized)
2. SHELX binaries absent → copied x86_64 shelx{c,d,e} from ccp4-9 (runs via Rosetta).
   Strategic: SHELX is x86_64-only, Rosetta 2 is on Apple's deprecation path.
3. modelcraft package truncated (134/173 files) → reinstalled modelcraft==6.1.1.

Open item (a decision, not a clear bug):
- test_shelx::test_gamma_sad / test_gamma_siras — **deterministic** R-work
  0.242 / 0.244 across 3 runs vs thresholds 0.22 / 0.23. Not flaky; a consistent
  ~0.02 offset vs old CCP4, most likely the new refmac/servalcat refinement.
  Decide: investigate the refinement delta, or relax thresholds (R~0.24 is
  reasonable for these SAD/SIRAS experimental-phasing cases).

Not investigated:
- test_nucleofind::test_1hr2 / test_1hr2_raw (SystemExit 2) — pre-existing,
  likely a NucleoFind ML-model/tooling gap (another env-side item).

---

## CORRECTION (2026-06-13) — bundle "defects" are mostly a local bad unpack

Inspecting `ccp4-20260520-linux.tar.gz` (the Linux build of the SAME version)
changes the conclusion for two of the three "bundle" items:

- **modelcraft is COMPLETE in the Linux bundle** (80 .py files, all key files
  present, at `lib/python3.11/site-packages/modelcraft/`). The truncation
  (13/80 files) is specific to the local **Mac/arm64 install**. The disk was at
  **99% (≈6 GB free)** during this session — a partially-extracted package is the
  classic signature of an extraction that ran out of disk. The ~12 empty-metadata
  packages are consistent with the same cause.
  → NOT a CCP4 build defect. Action: re-extract the Mac bundle on a disk with
    headroom and re-verify, rather than reporting it as a build bug.

- **SHELX**: absent from the bin/ of this bundle (see note for the Linux check).
  This is the one genuinely build/distribution-level item, separate from the
  unpack issue. SHELX is x86_64-only (no arm64 build) and runs via Rosetta 2,
  which is on Apple's deprecation path — that's the real question for CCP4.

So the only thing worth raising with CCP4 maintainers as a build issue is SHELX
(+ the arm64/Rosetta strategy). The modelcraft/metadata breakage was almost
certainly this machine's full-disk unpack.

**Further correction:** SHELX absence is also NOT a CCP4 build issue — SHELX is
installed separately under its own (Sheldrick) license. So nothing here is
actionable for CCP4 maintainers; the only real changes were the two ccp4i2 fixes.
