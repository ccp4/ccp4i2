# Package Management Summary

This document provides a quick reference for the hybrid package management strategy used in ccp4i2.

## Environment
- **Python**: 3.11.11
- **CCP4 Distribution**: CCP4-20251105 (`$CCP4`)
- **Virtual Environment**: `.venv/` (Python 3.11)

## Pip-Installed Packages (56 total)

Installed via `pip install -r requirements.txt`:

### Web Framework & API
- `Django==5.2.8`
- `djangorestframework==3.16.1`
- `django-cors-headers==4.9.0`
- `django-filter==25.2`
- `whitenoise==6.11.0`

### Scientific Computing
- `numpy==1.26.4` (pinned < 2.0 for CCP4 compatibility)
- `scipy==1.16.3`
- `pandas==2.3.3`
- `matplotlib==3.10.7`

### Structural Biology (PyPI)
- `biopython==1.86`
- `gemmi==0.7.3`
- `mrcfile==1.5.4`
- `cctbx-base==2025.10` (CCTBX utilities, not full stack)
  - ⚠️ **WARNING**: Installs broken MolProbity executables that must be removed (see below)

### Testing
- `pytest==8.4.2`
- `pytest-django==4.11.1`
- `pytest-asyncio==1.2.0`
- `pytest-xdist==3.8.0` (parallel testing)
- `pytest-ordering==0.6`

### Code Quality
- `autopep8==2.3.2`
- `mypy==1.18.2`

### Utilities
- `lxml==6.0.2`
- `PyYAML==6.0.3`
- `requests==2.32.5`
- `svgwrite==1.4.3`
- `psutil==7.1.3`

**Full list**: Run `.venv/bin/pip list --format=freeze`

## Symlinked CCP4 Modules (17 required)

Symlinked from CCP4-20251105 distribution at:
`$CCP4/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/`

**IMPORTANT**: For MolProbity validation to work, you must also symlink `share/cctbx` (see section 9 below).

### 1. Clipper (Crystallography)
```bash
ln -sf $CCP4_SITE_PACKAGES/clipper.py .venv/lib/python3.11/site-packages/
ln -sf $CCP4_SITE_PACKAGES/_clipper.so .venv/lib/python3.11/site-packages/
```
**Purpose**: Crystallographic calculations (FFT, map operations)
**Used by**: refmac, ctruncate, freerflag tests

### 2. CCP4MG (Molecular Graphics + MMDB2)
```bash
ln -sf $CCP4_SITE_PACKAGES/ccp4mg .venv/lib/python3.11/site-packages/
```
**Purpose**: Macromolecular database (mmdb2), coordinate manipulation
**Used by**: Model validation, structure manipulation tests

### 3. PyRVAPI (Report Viewing)
```bash
ln -sf $CCP4_SITE_PACKAGES/pyrvapi.so .venv/lib/python3.11/site-packages/
ln -sf $CCP4_SITE_PACKAGES/pyrvapi_ext .venv/lib/python3.11/site-packages/
```
**Purpose**: Interactive HTML report generation
**Used by**: All plugins with GUI reports

### 4. RDKit (Chemistry Toolkit) ⚠️ CRITICAL
```bash
ln -sf $CCP4_SITE_PACKAGES/rdkit .venv/lib/python3.11/site-packages/
```
**Purpose**: Chemical structure manipulation, SMILES parsing
**Version**: v2023.03.3 (CCP4-bundled)
**Why symlink?**:
- PyPI version v2025.09.1 incompatible with phaser's pickle implementation
- Phaser crashes with: `RuntimeError: Pickling of "rdkit.rdBase._vectd" instances is not enabled`
**Used by**: acedrg tests (8 tests), phaser tests (18 tests)

### 5. Phaser (Molecular Replacement)
```bash
ln -sf $CCP4_SITE_PACKAGES/phaser .venv/lib/python3.11/site-packages/
```
**Purpose**: Crystallographic phasing, molecular replacement
**Used by**: phaser_simple (4 tests), phaser_expert (2 tests), phaser_ep (1 test)

### 6. MrBUMP (MR Pipeline)
```bash
ln -sf $CCP4_SITE_PACKAGES/mrbump .venv/lib/python3.11/site-packages/
```
**Purpose**: Automated molecular replacement pipeline
**Used by**: mrbump test (1 test)

### 7. Iris Validation (Structure Validation)
```bash
ln -sf $CCP4_SITE_PACKAGES/iris_validation .venv/lib/python3.11/site-packages/
```
**Purpose**: Structure quality assessment
**Used by**: validate_protein tests (7 tests)

### 8. Chem_data (Top8000 Database)
```bash
ln -sf $CCP4_SITE_PACKAGES/chem_data .venv/lib/python3.11/site-packages/
```
**Purpose**: Ramachandran and rotamer reference database (Top8000)
**Used by**: MolProbity validation, iris_validation

### 9. CCTBX Suite (Required for MolProbity) ⚠️ CRITICAL
```bash
# Symlink CCTBX Python modules
ln -s $CCP4_SITE_PACKAGES/cctbx .venv/lib/python3.11/site-packages/
ln -s $CCP4_SITE_PACKAGES/mmtbx .venv/lib/python3.11/site-packages/
ln -s $CCP4_SITE_PACKAGES/iotbx .venv/lib/python3.11/site-packages/
ln -s $CCP4_SITE_PACKAGES/scitbx .venv/lib/python3.11/site-packages/
ln -s $CCP4_SITE_PACKAGES/smtbx .venv/lib/python3.11/site-packages/
ln -s $CCP4_SITE_PACKAGES/boost_adaptbx .venv/lib/python3.11/site-packages/

# Symlink CCTBX build environment (contains libtbx_env)
ln -sf $CCP4/Frameworks/Python.framework/Versions/3.11/share/cctbx .venv/share/cctbx
```
**Purpose**: Computational Crystallography Toolbox - core libraries for MolProbity
**Why critical?**:
- MolProbity validation requires `libtbx.env.has_module("probe")` to return True
- `libtbx.env` is loaded via `import libtbx.load_env` which requires:
  - Full CCTBX suite modules (cctbx, mmtbx, iotbx, scitbx, libtbx, smtbx, boost_adaptbx)
  - Build environment file at `share/cctbx/libtbx_env`
- Without these, iris_validation falls back to basic validation (no MolProbity data)
**Used by**: validate_protein tests with MolProbity (3 tests)

## Verification

### Check Pip Packages
```bash
source .venv/bin/activate
pip list | wc -l  # Should show ~56 packages
```

### Check Symlinked Modules
```bash
ls -la .venv/lib/python3.11/site-packages/ | grep "^l" | wc -l  # Should show 17 symlinks
```

### Test Imports
```bash
# Core crystallography
python -c "import clipper; from ccp4mg import mmdb2; print('✓ clipper, mmdb2')"

# Critical RDKit version
python -c "import rdkit; print(f'RDKit: {rdkit.__version__}')"  # Must show 2023.03.3

# Molecular replacement
python -c "import phaser; import mrbump; print('✓ phaser, mrbump')"

# Validation
python -c "import iris_validation; print('✓ iris_validation')"

# CCTBX suite and MolProbity probe detection
python -c "import libtbx.load_env; print(f'✓ libtbx.env loaded, has probe: {libtbx.env.has_module(\"probe\")}')"  # Must show True
```

## Test Results

With all 56 pip packages + 17 symlinked modules:
- **53 passed** out of 69 tests
- **77% pass rate**
- Runtime: ~38 minutes
- Commit: `9187a41`

## Quick Setup

```bash
# 1. Create venv with CCP4's Python 3.11
$CCP4/Frameworks/Python.framework/Versions/3.11/bin/python3.11 -m venv .venv

# 2. Install pip packages
source .venv/bin/activate
pip install -r requirements.txt

# 3. Run post-install cleanup (removes broken MolProbity executables)
./setup_venv.sh

# 4. Symlink CCP4 modules
CCP4_SITE_PACKAGES="$CCP4/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages"
cd .venv/lib/python3.11/site-packages

# Core modules (10 symlinks)
ln -sf "$CCP4_SITE_PACKAGES/clipper.py" .
ln -sf "$CCP4_SITE_PACKAGES/_clipper.so" .
ln -sf "$CCP4_SITE_PACKAGES/ccp4mg" .
ln -sf "$CCP4_SITE_PACKAGES/pyrvapi.so" .
ln -sf "$CCP4_SITE_PACKAGES/pyrvapi_ext" .
ln -sf "$CCP4_SITE_PACKAGES/rdkit" .
ln -sf "$CCP4_SITE_PACKAGES/phaser" .
ln -sf "$CCP4_SITE_PACKAGES/mrbump" .
ln -sf "$CCP4_SITE_PACKAGES/iris_validation" .
ln -sf "$CCP4_SITE_PACKAGES/chem_data" .

# CCTBX suite (6 symlinks - required for MolProbity validation)
ln -s "$CCP4_SITE_PACKAGES/cctbx" .
ln -s "$CCP4_SITE_PACKAGES/mmtbx" .
ln -s "$CCP4_SITE_PACKAGES/iotbx" .
ln -s "$CCP4_SITE_PACKAGES/scitbx" .
ln -s "$CCP4_SITE_PACKAGES/smtbx" .
ln -s "$CCP4_SITE_PACKAGES/boost_adaptbx" .

# CCTBX build environment (1 symlink - required for libtbx.env)
cd $CCP4I2_ROOT
ln -sf $CCP4/Frameworks/Python.framework/Versions/3.11/share/cctbx .venv/share/cctbx

# 5. Verify
python -c "import clipper, rdkit, phaser; print('✅ Core modules OK')"
python -c "import libtbx.load_env; print(f'✅ CCTBX OK, has probe: {libtbx.env.has_module(\"probe\")}')"
```

## Common Issues

### RDKit Version Conflict
**Symptom**: `RuntimeError: Pickling of "rdkit.rdBase._vectd" instances is not enabled`

**Cause**: Pip-installed RDKit v2025.09.1 instead of CCP4's v2023.03.3

**Fix**:
```bash
rm -rf .venv/lib/python3.11/site-packages/rdkit*
ln -sf "$CCP4_SITE_PACKAGES/rdkit" .venv/lib/python3.11/site-packages/
```

### Module Not Found
**Symptom**: `ModuleNotFoundError: No module named 'phaser'`

**Fix**: Check symlink exists and points to correct location:
```bash
ls -la .venv/lib/python3.11/site-packages/phaser
```

If broken, re-create symlink as shown above.

### MolProbity Validation Tests Fail (⚠️ CRITICAL)
**Symptom**:
- Validation tests with MolProbity fail
- Error: `WARNING: Failed to run MolProbity; continuing without MolProbity analyses`
- Test assertion: `MolProbity tag missing in program.xml`

**Root Cause #1 - Broken MolProbity Executables**:
The `cctbx-base==2025.10` package from PyPI installs 21 broken MolProbity executables in `.venv/bin/` with hardcoded paths to conda build environments (`/Users/runner/miniforge3/conda-bld/...`). These broken executables shadow the working CCP4 MolProbity tools in `$CCP4/bin/`.

**Root Cause #2 - Missing CCTBX Suite** (ACTUAL BLOCKER):
MolProbity validation requires `libtbx.env.has_module("probe")` to return True. This requires:
1. Full CCTBX suite modules symlinked (cctbx, mmtbx, iotbx, scitbx, smtbx, boost_adaptbx)
2. CCTBX build environment symlinked at `.venv/share/cctbx` (contains `libtbx_env` file)

**Fix**:
```bash
# 1. Remove broken MolProbity executables from virtualenv
rm .venv/bin/molprobity.*

# 2. Symlink CCTBX suite modules (see section 9 under "Symlinked CCP4 Modules")
CCP4_SITE="$CCP4/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages"
cd .venv/lib/python3.11/site-packages
ln -s "$CCP4_SITE/cctbx" .
ln -s "$CCP4_SITE/mmtbx" .
ln -s "$CCP4_SITE/iotbx" .
ln -s "$CCP4_SITE/scitbx" .
ln -s "$CCP4_SITE/smtbx" .
ln -s "$CCP4_SITE/boost_adaptbx" .

# 3. Symlink CCTBX build environment
cd $CCP4I2_ROOT
ln -sf $CCP4/Frameworks/Python.framework/Versions/3.11/share/cctbx .venv/share/cctbx

# 4. Verify MolProbity detection works
source /path/to/ccp4-20251105/bin/ccp4.setup-sh
source .venv/bin/activate
python -c "import libtbx.load_env; print(f'has probe: {libtbx.env.has_module(\"probe\")}')"
# Should output: has probe: True
```

**Prevention**:
1. Run the automated cleanup script after pip install: `./setup_venv.sh`
2. Follow the complete symlink setup in "Quick Setup" section above (includes CCTBX suite)
