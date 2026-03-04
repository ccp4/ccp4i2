# Testing Guide

## Running i2run Tests

### Using the Test Runner

The `run_test.sh` script sets up the environment and runs pytest with proper configuration.

#### Run all tests in a file:
```bash
./run_test.sh tests/i2run/test_parrot.py -v
```

#### Run a specific test:
```bash
./run_test.sh tests/i2run/test_refmac.py::test_8xfm_basic -v
```

#### Run tests in parallel:
```bash
./run_test.sh tests/i2run/ -n 4
```

### Manual pytest Invocation

If you prefer to run pytest directly, set up the environment first:

```bash
source /path/to/ccp4-20251105/bin/ccp4.setup-sh
export CCP4I2_ROOT=$(pwd)
export DJANGO_SETTINGS_MODULE=ccp4i2.config.test_settings

# Then run pytest
ccp4-python -m pytest tests/i2run/test_parrot.py -v -s
```

### Test Output Location

Test projects are created in `~/.cache/ccp4i2-tests/` with timestamped names:
```
~/.cache/ccp4i2-tests/
├── 20251204_113645_968f_refmac_test_8xfm_basic/   # Failed test (preserved)
├── 20251204_114012_a2b3_parrot_test_parrot/       # Passed test (cleaned up)
└── ...
```

- **Passed tests**: Automatically cleaned up
- **Failed tests**: Preserved for debugging
- **Cleanup**: `rm -rf ~/.cache/ccp4i2-tests/*`

### Available Test Suites

- **i2run tests**: `tests/i2run/` - Integration tests that run crystallographic tasks
- **Inheritance tests**: `tests/test_prosmart_refmac_inheritance.py`
- **Hash collision fix**: `tests/test_hash_collision_fix.py`

### Test Categories

| Test | Description | Duration |
|------|-------------|----------|
| `test_parrot.py` | Density modification | ~3s |
| `test_sheetbend.py` | Model rebuilding | ~3s |
| `test_coordinate_selector.py` | Coordinate selection | ~1s |
| `test_servalcat.py` | Refinement pipeline | ~7s |
| `test_refmac.py` | Refmac refinement | ~10s |
| Full suite | All i2run tests | ~40min |

### Environment Variables

- `CCP4I2_ROOT`: Project root directory (set by `run_test.sh`)
- `DJANGO_SETTINGS_MODULE`: Django settings (set by `run_test.sh` to `ccp4i2.config.test_settings`)
- `DEBUG_MERGE`: Set to `1` to enable debug output for container merging

### Compounds App Tests

The compounds app has its own Django settings:

```bash
cd server
source /path/to/ccp4-20251105/bin/ccp4.setup-sh
PYTHONPATH="$PWD:$PWD/../apps" DJANGO_SETTINGS_MODULE=compounds.settings \
  ccp4-python -m pytest ../apps/compounds/assays/tests/test_aggregation.py -v
```

### Troubleshooting

If tests fail with import errors:
1. Ensure CCP4 environment is sourced: `source /path/to/ccp4-20251105/bin/ccp4.setup-sh`
2. Ensure ccp4i2 is pip-installed: `cd server && ccp4-python -m pip install -e .`
3. Check CCP4I2_ROOT is set: `echo $CCP4I2_ROOT`

If Django tests fail:
1. Ensure DJANGO_SETTINGS_MODULE is set
2. Check database migrations are up to date

If `fixture 'django_db_blocker' not found`:
```bash
ccp4-python -m pip install pytest-django pytest-xdist
```
