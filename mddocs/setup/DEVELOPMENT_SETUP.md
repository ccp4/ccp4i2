# Development Environment Setup

This guide walks through setting up a development environment for ccp4i2-django with full CCP4 integration.

## Prerequisites

- **CCP4 Suite 2024/2025** - Download from [CCP4 Downloads](https://ccp4serv6.rc-harwell.ac.uk/10/)
- **Git** for version control
- **macOS, Linux, or WSL** (Unix-like environment)

## Quick Start (Recommended)

The simplest approach uses `ccp4-python` directly, which includes all crystallographic libraries.

### Step 1: Install CCP4

Download and install CCP4 from https://ccp4serv6.rc-harwell.ac.uk/10/

Common installation locations:
- **macOS**: `/Applications/ccp4-20251105` or `~/Developer/ccp4-20251105`
- **Linux**: `/opt/ccp4-20251105` or `~/ccp4-20251105`

### Step 2: Clone the Repository

```bash
git clone https://github.com/ccp4/ccp4i2.git
cd ccp4i2
git checkout ccp4i2-django
```

### Step 3: Source CCP4 Environment

```bash
# macOS/Linux
source /path/to/ccp4-20251105/bin/ccp4.setup-sh

# Verify ccp4-python is available
which ccp4-python
ccp4-python --version  # Should show Python 3.11.x
```

### Step 4: Install Django Dependencies

CCP4 includes Django, but we need a few additional packages:

```bash
# Core Django extensions
ccp4-python -m pip install django-cors-headers django-filter djangorestframework whitenoise

# For running the Electron client (ASGI server)
ccp4-python -m pip install uvicorn==0.20.0
```

### Step 5: Verify Setup

```bash
# Test Django
ccp4-python -c "import django; print(f'Django {django.__version__}')"

# Test baselayer (Qt-free compatibility)
ccp4-python -c "from baselayer import DJANGO; print(f'DJANGO mode: {DJANGO()}')"

# Test crystallographic libraries
ccp4-python -c "import clipper; import phaser; print('CCP4 libraries OK')"
```

### Step 6: Run Tests

```bash
# Run a quick test
./run_test.sh tests/i2run/test_parrot.py -v

# Run multiple tests in parallel
./run_test.sh tests/i2run/ -n 4
```

Test projects are stored in `~/.cache/ccp4i2-tests/`. Cleanup instructions are printed after each test run.

---

## Alternative: Virtual Environment with Symlinks

For development requiring packages not in CCP4, you can create a Python virtual environment and symlink CCP4 modules into it.

### When to Use This Approach

- You need specific package versions different from CCP4
- You're developing features that require additional Python packages
- You want isolation from CCP4's Python environment

### Setup Steps

1. **Create virtual environment using CCP4's Python**:
   ```bash
   CCP4_PYTHON="/path/to/ccp4-20251105/Frameworks/Python.framework/Versions/3.11/bin/python3.11"
   $CCP4_PYTHON -m venv .venv
   source .venv/bin/activate
   ```

2. **Install requirements**:
   ```bash
   pip install -r requirements.txt
   pip install "numpy<2"  # CCP4 modules require NumPy 1.x
   ```

3. **Symlink CCP4 modules**:
   ```bash
   CCP4_SITE="/path/to/ccp4-20251105/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages"
   VENV_SITE=".venv/lib/python3.11/site-packages"

   # Core modules
   ln -sf $CCP4_SITE/{clipper.py,_clipper.so,ccp4mg,pyrvapi.so,pyrvapi_ext} $VENV_SITE/

   # Crystallography
   ln -sf $CCP4_SITE/{phaser,mrbump,rdkit,iris_validation,chem_data} $VENV_SITE/

   # CCTBX suite (for MolProbity validation)
   ln -sf $CCP4_SITE/{cctbx,mmtbx,iotbx,scitbx,libtbx,smtbx,boost_adaptbx} $VENV_SITE/

   # CCTBX build environment
   mkdir -p .venv/share
   ln -sf /path/to/ccp4-20251105/Frameworks/Python.framework/Versions/3.11/share/cctbx .venv/share/
   ```

4. **Verify**:
   ```bash
   python -c "import clipper; import phaser; print('OK')"
   ```

---

## Environment Configuration

### Environment Variables

The test runner (`run_test.sh`) automatically sets these, but for manual testing:

```bash
export CCP4I2_ROOT=/path/to/ccp4i2           # Project root
export DJANGO_SETTINGS_MODULE=ccp4x.config.test_settings
export PYTHONPATH=$CCP4I2_ROOT:$CCP4I2_ROOT/server:$PYTHONPATH
```

### Optional .env File

Create a `.env` file in the project root for local configuration:

```bash
CCP4_ROOT=/path/to/ccp4-20251105
CCP4_VERSION=ccp4-20251105
PYTHON_VERSION=3.11
```

---

## Running Tests

### Test Runner

The `run_test.sh` script handles environment setup:

```bash
# Single test
./run_test.sh tests/i2run/test_parrot.py -v

# Test file
./run_test.sh tests/i2run/test_servalcat.py

# All i2run tests (parallel)
./run_test.sh tests/i2run/ -n 4

# Specific test function
./run_test.sh tests/i2run/test_refmac.py::test_8xfm_basic -v
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

### Test Categories

| Test | Description | Duration |
|------|-------------|----------|
| `test_parrot.py` | Density modification | ~3s |
| `test_sheetbend.py` | Model rebuilding | ~3s |
| `test_coordinate_selector.py` | Coordinate selection | ~1s |
| `test_servalcat.py` | Refinement pipeline | ~7s |
| `test_refmac.py` | Refmac refinement | ~10s |
| Full suite | All i2run tests | ~40min |

---

## Troubleshooting

### ModuleNotFoundError: No module named 'corsheaders'

Install missing Django packages:
```bash
ccp4-python -m pip install django-cors-headers django-filter djangorestframework whitenoise
```

### No module named 'uvicorn' (Electron client)

The Electron client uses uvicorn as the ASGI server:
```bash
ccp4-python -m pip install uvicorn==0.20.0
```

### NumPy Version Conflicts

CCP4 modules require NumPy 1.x:
```bash
pip install "numpy<2"
```

### Segmentation Fault on Import

Python version mismatch. Ensure you're using Python 3.11 matching CCP4:
```bash
ccp4-python --version  # Must be 3.11.x
```

### libtbx.env FileNotFoundError

Missing CCTBX share directory symlink:
```bash
mkdir -p .venv/share
ln -sf /path/to/ccp4/share/cctbx .venv/share/
```

### Django Database Errors

Ensure settings module is configured:
```bash
export DJANGO_SETTINGS_MODULE=ccp4x.config.test_settings
```

---

## Platform Notes

### macOS

- CCP4 Python: `$CCP4_ROOT/Frameworks/Python.framework/Versions/3.11/bin/python3.11`
- Site-packages: `$CCP4_ROOT/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages`

### Linux

- CCP4 Python: `$CCP4_ROOT/bin/ccp4-python`
- Site-packages: `$CCP4_ROOT/lib/python3.11/site-packages`

### WSL (Windows)

Follow Linux instructions. CCP4 paths are Linux-style within WSL.

---

## Next Steps

After setup:

1. **Run quick tests** to verify everything works:
   ```bash
   ./run_test.sh tests/i2run/test_parrot.py tests/i2run/test_sheetbend.py -v
   ```

2. **Read the architecture docs**:
   - [Migration Strategy](../MIGRATION_STRATEGY.md)
   - [Quick Reference](../QUICK_REFERENCE.md)

3. **Understand the baselayer**:
   - Qt-free compatibility layer at `baselayer/`
   - Automatic Qt vs Django mode detection
   - Use `from baselayer import QtCore` instead of `from PySide2 import QtCore`
