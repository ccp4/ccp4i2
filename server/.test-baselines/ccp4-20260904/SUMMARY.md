# i2run baseline — ccp4-20260904 (September 2026 CCP4 10 build, macOS arm64)

- **Date:** 2026-09-25/26
- **Commit:** 82b1a45fa (django, release 3.1.0a79) plus the Coot
  molecules-container fix that became #621
- **Result:** 224 tests · **203 passed · 0 failed · 21 skipped** · 60.3 min
- **Like-for-like:** `ccp4-20260702` run on the same commit, same day:
  **203 passed · 0 failed · 21 skipped**, identical skip set (see that
  label's SUMMARY.md). The September build introduces no i2run regression.
- **Stack:** Python 3.11.14, gemmi 0.7.5, Django 5.2.15 / DRF 3.17.1 from
  `requirements-runtime.txt`, pytest 7.3.2, pytest-django 4.14.0.

## What changed in this build that ccp4i2 had to absorb

- `coot_headless_api` keeps only `molecules_container_t`; the July build had
  `molecules_container_py` as well. Four call sites named the old class and
  failed minutes into a run. Fixed in `server/ccp4i2/lib/coot_api.py` (#621),
  which resolves the name and works on both builds.
- `bin/pandda2.analyse` is present (PanDDA 2 in its own micromamba env).

## Skips (21) — all environment gaps, none are failures

| Reason | Tests |
|---|---|
| shelxc / shelxd / shelxe not in the bundle | test_shelx (3), test_shelxe_mr (1) |
| xds_par not installed | test_xia2 (1) |
| Documented stubs with no body | test_xia2_ssx_reduce, test_xia2_multiplex, test_import_serial_pipe, test_pdb_redo_api |
| Historic tests needing test101 project zips | test_i2run |
| Opt-in PanDDA end-to-end (volume / `CCP4I2_PANDDA_E2E=1`) | test_pandda_events, test_pandda_campaign |
| Remaining | per-test `skipif` markers, identical on both builds |

## Slowest

servalcat neutron 316 s · zanuda 187 s · lorestr mmCIF 122 s · refmac 119 s ·
simbad lattice 113 s.

## How this was run

Fresh tarball → `./BINARY.setup` (creates `bin/ccp4.setup-sh`) → install the
ccp4i2 stack the way the desktop app does (`pip install --no-deps -e server`,
`--no-deps -e packages/ccp4i2-api`, `--no-deps --ignore-installed -r
server/ccp4i2/requirements-runtime.txt`, `--no-deps pytest-django`) → from
`server/`, `CCP4_SETUP=<tree>/bin/ccp4.setup-sh bash run_i2run_baseline.sh`.
pip's resolver crashes on this distribution's corrupt `typing_extensions`
dist-info, which is why every install step is `--no-deps`.
