# MolProbity Validation Test Debugging

## Issue Summary

Validation tests with MolProbity are failing with:
```
WARNING: Failed to run MolProbity; continuing without MolProbity analyses
AssertionError: MolProbity tag missing in program.xml
```

## Root Cause Identified

The `cctbx-base==2025.10` package from PyPI installed 21 broken MolProbity executables in `.venv/bin/` with hardcoded conda build paths. These shadowed the working CCP4 MolProbity tools.

**Fix Applied**: Removed broken executables with `rm .venv/bin/molprobity.*`

## Current Status

After removing broken executables:
- ✅ MolProbity tools point to correct CCP4 location
- ✅ `which molprobity.ramalyze` returns `/Users/nmemn/Developer/ccp4-20251105/bin/molprobity.ramalyze`
- ✅ Manual execution works: `molprobity.ramalyze <file>` succeeds
- ❌ Tests still fail: iris_validation unable to run MolProbity

## Hypothesis

There appears to be a deeper issue with how `iris_validation` invokes MolProbity that is unrelated to the PATH shadowing. Possible causes:

1. **Python subprocess environment**: iris_validation may spawn subprocesses that don't inherit the CCP4 environment
2. **Missing dependencies**: MolProbity tools may require additional Python packages or environment variables
3. **CCP4 version mismatch**: iris_validation from CCP4-20251105 may expect different MolProbity behavior
4. **Silent failures**: MolProbity may be failing for a different reason and iris_validation is catching the exception

## Testing Evidence

```bash
# Manual test - WORKS
source /Users/nmemn/Developer/ccp4-20251105/bin/ccp4.setup-sh
source .venv/bin/activate
cd tests/i2run/test_projects/.../CCP4_IMPORTED_FILES
molprobity.ramalyze tmpv9usefsg_8ola_final.cif
# Output: Ramachandran analysis completes successfully

# Pytest test - FAILS
./run_test.sh tests/i2run/test_validate.py::test_8ola_8olf_with_molprobity
# Output: WARNING: Failed to run MolProbity
```

## Comparison with Working Laptop

The tests work on your laptop, suggesting an environment difference. Possibilities:
1. Different CCP4 version
2. Different Python version
3. Different iris_validation version
4. Additional environment variables set
5. Different cctbx-base installation state

## Actions Taken

1. ✅ Documented cctbx-base MolProbity conflict in DEVELOPMENT_SETUP_SUMMARY.md
2. ✅ Added warning to requirements.txt
3. ✅ Created setup_venv.sh script to automate cleanup
4. ✅ Removed broken MolProbity executables from .venv/bin/

## Remaining Work

To fully resolve the MolProbity test failures, need to investigate:

1. **Compare environments**: Check differences between laptop and desktop:
   ```bash
   # On both machines
   pip show iris-validation
   ls -la .venv/lib/python3.11/site-packages/iris_validation/
   python -c "import iris_validation; print(iris_validation.__version__)"
   env | grep -E "(CCP4|PATH|PYTHON)" | sort
   ```

2. **Enable iris_validation debug logging**: Modify iris_validation calls to capture stderr/stdout from MolProbity subprocess calls

3. **Test MolProbity with same arguments**: Find out exactly what command iris_validation is running and replicate it manually

4. **Check for Python module conflicts**: iris_validation may be trying to import a Python module that conflicts with venv packages

## Workaround

Tests can pass the basic validation checks (without MolProbity) by using:
```python
check_program_xml_basic(xml)  # Instead of check_program_xml_with_molprobity(xml)
```

## Documentation Impact

- DEVELOPMENT_SETUP_SUMMARY.md now includes comprehensive troubleshooting for cctbx-base MolProbity conflict
- requirements.txt warns about the issue
- setup_venv.sh automates the cleanup

## Next Steps

1. Compare laptop vs desktop iris_validation versions and environments
2. Add debug logging to understand why iris_validation's MolProbity calls fail
3. Test hypothesis: Try running tests with laptop's iris_validation version
