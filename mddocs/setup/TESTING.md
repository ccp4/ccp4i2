# Testing Guide

## Running i2run Tests

### Using the Test Dispatcher

The `run_test.sh` script sets up the environment and runs pytest with proper configuration.

#### Run all tests in a file:
```bash
./run_test.sh server/ccp4x/tests/i2run/test_i2run.py
```

#### Run a specific test:
```bash
./run_test.sh server/ccp4x/tests/i2run/test_i2run.py test_case1
```

#### Run tests from other directories:
```bash
./run_test.sh tests/test_prosmart_refmac_inheritance.py
./run_test.sh tests/test_hash_collision_fix.py
```

### Manual pytest Invocation

If you prefer to run pytest directly, set up the environment first:

```bash
export CCP4I2_ROOT=$(pwd)
export PYTHONPATH=$(pwd):$(pwd)/server:$PYTHONPATH
export DJANGO_SETTINGS_MODULE=ccp4x.settings

# Then run pytest
python -m pytest server/ccp4x/tests/i2run/test_i2run.py -v -s
python -m pytest server/ccp4x/tests/i2run/test_i2run.py::TestI2Run::test_case1 -v -s
```

### Available Test Suites

- **i2run tests**: `server/ccp4x/tests/i2run/test_i2run.py`
  - `test_shlex` - Test argument parsing
  - `test_case1` - Test aimless_pipe
  - `test_case2` - Test aimless_pipe with reference
  - `test_case3` - Test complex pipeline

- **Inheritance tests**: `tests/test_prosmart_refmac_inheritance.py`
  - Verifies DEF XML inheritance works correctly
  - Tests that NCYCLES is inherited from refmac_i2

- **Hash collision fix**: `tests/test_hash_collision_fix.py`
  - Verifies CInt/CFloat objects with same value work correctly
  - Tests identity-based hashing

### Example: Testing prosmart_refmac with NCYCLES

To verify the hash collision fix works for prosmart_refmac:

```bash
# Test the inheritance mechanism
./run_test.sh tests/test_prosmart_refmac_inheritance.py

# Test the hash collision fix
./run_test.sh tests/test_hash_collision_fix.py
```

### Environment Variables

- `CCP4I2_ROOT`: Project root directory (required)
- `PYTHONPATH`: Must include project root and server directory
- `DJANGO_SETTINGS_MODULE`: Django settings (for i2run tests)
- `DEBUG_MERGE`: Set to `1` to enable debug output for container merging

### Troubleshooting

If tests fail with import errors:
1. Ensure virtual environment is activated: `source .venv/bin/activate`
2. Check CCP4I2_ROOT is set: `echo $CCP4I2_ROOT`
3. Verify PYTHONPATH includes server: `echo $PYTHONPATH`

If Django tests fail:
1. Ensure DJANGO_SETTINGS_MODULE is set
2. Check database migrations are up to date
3. Verify test database can be created

If i2run tests fail with "NOT NULL constraint failed: ccp4x_job.title":
- This has been fixed in `server/ccp4x/lib/utils/jobs/create.py`
- Jobs now fallback to using `taskName` if plugin doesn't define TASKTITLE
- All jobs will have a valid title set
