# Development Environment Setup

This guide walks through setting up a development environment for ccp4i2-django with full CCP4 integration.

## Prerequisites

- **CCP4 Suite 2024/2025** - Download from [CCP4 Downloads](https://ccp4serv6.rc-harwell.ac.uk/10/)
- **Node.js 18+** and **npm 9+** - For the Electron client (download from [nodejs.org](https://nodejs.org/))
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

### Step 4: Replace CCP4's ccp4i2 with This Branch

CCP4 ships with its own version of ccp4i2 in site-packages. To use this development branch instead:

```bash
# First, remove CCP4's bundled ccp4i2
CCP4_SITE=$(ccp4-python -c "import site; print(site.getsitepackages()[0])")
rm -rf "$CCP4_SITE/ccp4i2"

# Install this branch in editable mode with Django dependencies
ccp4-python -m pip install -e ".[django]"

# Or install everything (Django + test dependencies)
ccp4-python -m pip install -e ".[full]"

# Note: pip may show an error about "checking for conflicts" - this is a known
# issue with some CCP4 packages and can be ignored if the install succeeds.
```

This replaces CCP4's bundled ccp4i2 with a link to your development checkout, so changes you make are immediately reflected. It also installs/upgrades Django and related packages (djangorestframework, django-cors-headers, django-filter, whitenoise, uvicorn).

**To revert to CCP4's original version**: Re-install CCP4 or copy the ccp4i2 directory from a fresh CCP4 installation.

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

## CLI Tools (Django Management Commands)

CLI tools are implemented as Django management commands. Run them via `ccp4-python` with the `manage.py` script:

```bash
cd server
ccp4-python manage.py <command> [options]
```

### Available Commands

**Job Execution:**
- `i2run <task> [options]` - Run a task (equivalent to legacy i2run)
- `create_job <task>` - Create a job without running it
- `run_job <job_id>` - Run an existing job
- `execute_job <job_id>` - Execute a job (internal)

**Project Management:**
- `list_projects` - List all projects
- `list_project <project_id>` - Show project details
- `create_project <name>` - Create a new project
- `tree_project <project_id>` - Show project job tree

**Job Management:**
- `list_jobs [project_id]` - List jobs in a project
- `tree_job <job_id>` - Show job file tree
- `clone_job <job_id>` - Clone a job
- `set_job_status <job_id> <status>` - Update job status
- `set_job_parameter <job_id> <param> <value>` - Set job parameter
- `validate_job <job_id>` - Validate job configuration

**File Operations:**
- `list_files <job_id>` - List files in a job
- `list_filetypes` - List available file types
- `cat_job_file <job_id> <filename>` - Display job file contents
- `cat_project_file <project_id> <filename>` - Display project file contents
- `preview_file <file_id>` - Preview file contents

**Import/Export:**
- `export_job <job_id>` - Export job to archive
- `export_project <project_id>` - Export project to archive
- `import_ccp4_project_zip <zipfile>` - Import legacy CCP4 project
- `import_i2xml <xmlfile>` - Import i2 XML file

**Reports:**
- `get_job_report <job_id>` - Generate job report

### Examples

```bash
# Run a refinement task
cd server
ccp4-python manage.py i2run refmac --hklin data.mtz --xyzin model.pdb

# List all projects
ccp4-python manage.py list_projects

# Show job tree
ccp4-python manage.py tree_job 42
```

---

## Running the Electron Client

The Electron client provides a modern GUI for CCP4i2 with a React frontend and Django backend.

### Prerequisites

- **Node.js 18+** and **npm 9+** (check with `node --version` and `npm --version`)
- **CCP4 with ccp4-python** installed and accessible
- **ccp4i2 installed** with Django dependencies (Step 4: `ccp4-python -m pip install -e ".[django]"`)

### CCP4 Detection

The client automatically detects your CCP4 installation. In development mode, it searches:

1. **Sibling directories** - Scans `../` for `ccp4-*` folders (e.g., `../ccp4-20251105`)
2. **Standard locations** - `/Applications/ccp4-9` (macOS), `C:\CCP4\ccp4-9` (Windows), `/opt/ccp4` (Linux)

The first directory containing `bin/ccp4-python` is used. Newer versions are preferred (sorted by name descending).

**Note**: You do NOT need to source `ccp4.setup-sh` before running the client - the Electron app sets up the CCP4 environment internally.

### Running in Development Mode

```bash
cd client
rm -rf node_modules  # Clean install (first time or after package.json changes)
npm install
npm run start:electron
```

This will:
1. Install npm dependencies
2. Build the Electron main process
3. Start Next.js development server
4. Launch the Electron window

### Configuration

On first launch, the client opens a configuration page where you can:
- Verify/change the CCP4 installation path
- Set the projects directory
- Configure other settings

Settings are persisted in `electron-store` (platform-specific location).

### Environment File

The client uses `client/renderer/.env.local` for configuration (not tracked in git).
Copy from template if needed:

```bash
cp client/renderer/.env.local.template client/renderer/.env.local
```

Default configuration runs without authentication. See `.env.local.auth-example` for Azure AD setup.

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

### Electron client can't find CCP4

The client looks for `ccp4-python` in sibling directories first, then standard locations.
Verify your CCP4 installation:
```bash
ls ../ccp4-*/bin/ccp4-python  # Should show your CCP4 installation
```

If CCP4 is elsewhere, you can configure the path in the client's config page on first launch.

### npm install fails

Ensure you have Node.js 18+ and npm 9+:
```bash
node --version  # Should be v18.x or higher
npm --version   # Should be 9.x or higher
```

If using an older version, update Node.js from [nodejs.org](https://nodejs.org/).

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
