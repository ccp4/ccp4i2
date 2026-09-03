# Testing Guide

All tests run from the `server/` directory with `ccp4-python -m pytest`
(cross-platform: mac, Linux, Windows). Test paths are relative to the
`ccp4i2` package. `./run_test.sh` remains as a mac/linux convenience wrapper
that sources the environment first.

The authoritative description of the test tiers — what each layer tests, how
fast it is, and whether it needs CCP4 — is the **Tests** section of
[CLAUDE.md](../../CLAUDE.md). The short version:

```
ccp4i2/tests/
├── unit/          # fast, no CCP4 binaries (containers, mtz, pdb, phil, ...)
├── async/         # async execution infrastructure
├── db/            # database, project import/export
├── api/unit/      # REST endpoints via the Django test client (no CCP4)
├── api/e2e/       # full pipelines via the REST API (slow, needs CCP4)
├── parity/        # native gemmi/numpy ports vs the CCP4 binary they replaced
└── i2run/         # full task execution via CLI (slow, needs CCP4)
```

## Running tests

```bash
cd server
source /path/to/ccp4-<build>/bin/ccp4.setup-sh

# Fast unit tests (no CCP4 binaries needed)
ccp4-python -m pytest ccp4i2/tests/unit/ -v

# One end-to-end task test / one test function
ccp4-python -m pytest ccp4i2/tests/i2run/test_parrot.py -v
ccp4-python -m pytest ccp4i2/tests/i2run/test_servalcat.py::test_servalcat_basic -v

# In parallel
ccp4-python -m pytest ccp4i2/tests/i2run/ -n 4
```

## Reproducing CI locally (no CCP4)

CI runs `tests/unit/` and `tests/api/unit/` on stock Python with **no CCP4**,
as two separate pytest invocations:

```bash
cd server
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pytest ccp4i2/tests/unit/ -q
.venv/bin/python -m pytest ccp4i2/tests/api/unit/ -q
```

Running these suites under `ccp4-python` gives strictly *more* coverage, so it
cannot prove they are CCP4-free — use the venv for that.

## Test output

Failed-test projects are preserved in `~/.cache/ccp4i2-tests/` (never pruned —
delete old ones when disk gets tight). Downloaded fixtures are cached in
`~/.cache/ccp4i2-tests/downloads/`; files over 100 MB are not cached unless
`CCP4I2_TEST_CACHE_MAX_MB` is raised.

## Environment variables

- `DJANGO_SETTINGS_MODULE`: set by pytest.ini to `ccp4i2.config.test_settings`
- `CCP4I2_TEST_CACHE_MAX_MB` / `CCP4I2_TEST_REFETCH`: fixture-cache controls
- `DEBUG_MERGE=1`: debug output for container merging

## Compounds app tests

The compounds app was extracted to the `newcastleuniversity/materia`
repository (2026-04-29 cut); run its tests from a Materia checkout.

## Troubleshooting

Import errors: source the CCP4 environment
(`source /path/to/ccp4-<build>/bin/ccp4.setup-sh`) and install editable
(`cd server && ccp4-python -m pip install -e .`).

`fixture 'django_db_blocker' not found`:
```bash
ccp4-python -m pip install pytest-django pytest-xdist
```
