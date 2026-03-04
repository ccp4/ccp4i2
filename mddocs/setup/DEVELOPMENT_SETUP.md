# Development Environment Setup

This guide walks through setting up a development environment for CCP4i2 with full CCP4 integration.

## Prerequisites

- **CCP4 Development Build** - Download from [CCP4 Development Builds](https://ccp4serv6.rc-harwell.ac.uk/10/downloads/packages/)
- **Node.js 18+** and **npm 9+** - For the Electron client (download from [nodejs.org](https://nodejs.org/))
- **Git** for version control
- **macOS, Linux, or WSL** (Unix-like environment)

## Quick Start (Recommended)

CCP4i2 is installed as a pip package into `ccp4-python`, which already includes all crystallographic libraries.

### Step 1: Install CCP4

Download a development build from https://ccp4serv6.rc-harwell.ac.uk/10/downloads/packages/

Common installation locations:
- **macOS**: `/Applications/ccp4-20251105` or `~/Developer/ccp4-20251105`
- **Linux**: `/opt/ccp4-20251105` or `~/ccp4-20251105`

### Step 2: Clone the Repository

```bash
git clone https://github.com/ccp4/ccp4i2.git
cd ccp4i2
```

### Step 3: Source CCP4 Environment

```bash
# macOS/Linux
source /path/to/ccp4-20251105/bin/ccp4.setup-sh

# Verify ccp4-python is available
which ccp4-python
ccp4-python --version  # Should show Python 3.11.x
```

### Step 4: Install ccp4i2 into ccp4-python

CCP4i2 is defined as a pip-installable package in `server/pyproject.toml`. Install it in editable mode so your source changes are immediately reflected:

```bash
cd server
ccp4-python -m pip install -e .
cd ..
```

This installs ccp4i2 and its dependencies (Django, djangorestframework, uvicorn, etc.) into `ccp4-python`'s site-packages. If CCP4 ships with an older bundled ccp4i2, the editable install takes precedence.

> **Note**: pip may show warnings about dependency conflicts with existing CCP4 packages — these can generally be ignored if the install succeeds.

### Step 5: Verify Setup

```bash
# Test ccp4i2 is importable
ccp4-python -c "import ccp4i2; print('ccp4i2 OK')"

# Test Django
ccp4-python -c "import django; print(f'Django {django.__version__}')"

# Test crystallographic libraries
ccp4-python -c "import clipper; import phaser; print('CCP4 libraries OK')"
```

### Step 6: Run Tests

```bash
# Run a quick test
./run_test.sh tests/i2run/test_parrot.py -v

# Run multiple tests in parallel (requires pytest-xdist, installed by pip)
./run_test.sh tests/i2run/ -n 4
```

Test projects are stored in `~/.cache/ccp4i2-tests/`. Cleanup instructions are printed after each test run.

---

## CLI Tools

The `i2` command provides a clean interface for CCP4i2 operations. It wraps Django management commands with intuitive syntax.

### Quick Reference

```bash
# Projects
i2 projects                     # List all projects
i2 projects create myproject    # Create a project
i2 projects show <id>           # Show project details
i2 projects tree <id>           # Show project job tree

# Jobs
i2 jobs <project>               # List jobs in a project
i2 jobs create <project> <task> # Create a job
i2 jobs run <job_id>            # Run an existing job
i2 jobs tree <job_id>           # Show job file tree
i2 jobs clone <job_id>          # Clone a job

# Run a task directly (create + run)
i2 run <task> [options]         # e.g., i2 run refmac --hklin data.mtz

# Files
i2 files <job_id>               # List files in a job
i2 files cat <job_id> <name>    # Display file contents

# Reports
i2 report <job_id>              # Generate job report

# Import/Export
i2 export job <job_id>          # Export job to archive
i2 export project <project_id>  # Export project
i2 import <zipfile>             # Import legacy CCP4 project
```

### Django Management Commands

For advanced usage, Django management commands are available via `python -m django`:

```bash
# From anywhere (with DJANGO_SETTINGS_MODULE set)
ccp4-python -m django <command> [options]

# Or from the server/ directory
cd server
ccp4-python manage.py <command> [options]
```

Commands include: `i2run`, `list_projects`, `list_jobs`, `create_job`, `run_job`, `tree_job`, `export_job`, `import_ccp4_project_zip`, and more. Run `ccp4-python -m django --help` for the full list.

---

## Running the Electron Client

The Electron client provides a modern desktop GUI for CCP4i2 with a React frontend and Django backend.

### Prerequisites

- **Node.js 18+** and **npm 9+** (check with `node --version` and `npm --version`)
- **CCP4 with ccp4-python** installed and accessible
- **ccp4i2 pip-installed** into ccp4-python (Step 4 above)

### How it Works

The Electron app:
- Detects your CCP4 installation automatically
- Uses `ccp4-python` to run Django migrations and start the Uvicorn ASGI server
- Loads the ASGI app via the module path `ccp4i2.config.asgi:application` — no Python source code is bundled in the Electron app

### CCP4 Detection

In development mode, the client searches:
1. **Sibling directories** - Scans `../` for `ccp4-*` folders (e.g., `../ccp4-20251105`)
2. **Standard locations** - `/Applications/ccp4-9` (macOS), `C:\CCP4\ccp4-9` (Windows), `/opt/ccp4` (Linux)

The first directory containing `bin/ccp4-python` is used. Newer versions are preferred.

**Note**: You do NOT need to source `ccp4.setup-sh` before running the client — the Electron app sets up the CCP4 environment internally.

### Running in Development Mode

```bash
cd client
npm install
npm run start:electron
```

This will:
1. Build the Electron main process
2. Start Next.js development server
3. Launch the Electron window
4. Start Django/Uvicorn via `ccp4-python`

### Building a Packaged App

```bash
cd client
npm run package-mac        # macOS .dmg
npm run package-win        # Windows .exe
npm run package-linux-x64  # Linux .AppImage
```

Packaged apps are written to `client/release/`.

**Important**: The packaged app requires `ccp4i2` to be pip-installed in the target machine's CCP4 installation. Without it, the Django backend will fail to start with `ModuleNotFoundError: No module named 'ccp4i2'`.

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

## Environment Configuration

### Environment Variables

After pip installation (`ccp4-python -m pip install -e .`), Python module discovery is handled automatically. The key environment variables are:

```bash
export CCP4I2_ROOT=/path/to/ccp4i2           # For resource files (set by run_test.sh)
export DJANGO_SETTINGS_MODULE=ccp4i2.config.settings  # Or test_settings for tests
```

The test runner (`run_test.sh`) and the Electron app set these automatically.

### Optional .env File

Create a `.env` file in the project root for local CCP4 configuration:

```bash
CCP4_ROOT=/path/to/ccp4-20251105
```

This is used by `run_test.sh` to find CCP4 if not already in the environment.

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

2. **Install ccp4i2**:
   ```bash
   cd server
   pip install -e .
   pip install "numpy<2"  # CCP4 modules require NumPy 1.x
   cd ..
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

## Troubleshooting

### ModuleNotFoundError: No module named 'ccp4i2'

ccp4i2 needs to be pip-installed into ccp4-python:
```bash
cd server
ccp4-python -m pip install -e .
```

### ModuleNotFoundError: No module named 'corsheaders'

This is installed automatically by `pip install -e .`. If missing, reinstall:
```bash
cd server
ccp4-python -m pip install -e .
```

### Electron client can't find CCP4

The client looks for `ccp4-python` in sibling directories first, then standard locations.
Verify your CCP4 installation:
```bash
ls ../ccp4-*/bin/ccp4-python  # Should show your CCP4 installation
```

If CCP4 is elsewhere, you can configure the path in the client's config page on first launch.

### Electron app starts but Django fails

Ensure ccp4i2 is pip-installed into the CCP4 installation that Electron detects:
```bash
source /path/to/ccp4-20251105/bin/ccp4.setup-sh
cd server
ccp4-python -m pip install -e .
```

### npm install fails

Ensure you have Node.js 18+ and npm 9+:
```bash
node --version  # Should be v18.x or higher
npm --version   # Should be 9.x or higher
```

### NumPy Version Conflicts

CCP4 modules require NumPy 1.x:
```bash
ccp4-python -m pip install "numpy<2"
```

### Coot Python 2 Not Found

Some CCP4i2 wrappers (e.g., `coot_find_waters`, `coot_script_lines`) require the legacy Python 2 version of Coot. Recent CCP4 builds may not include `coot_py2`. If you see errors like:

```
exec: .../coot_py2/bin/coot: No such file or directory
```

You need to symlink `coot_py2` from an older CCP4 installation (e.g., CCP4 9):

```bash
# macOS example
ln -s /Applications/ccp4-9/coot_py2 /path/to/ccp4-20251105/coot_py2

# Linux example
ln -s /opt/ccp4-9/coot_py2 /path/to/ccp4-20251105/coot_py2
```

### Segmentation Fault on Import

Python version mismatch. Ensure you're using Python 3.11 matching CCP4:
```bash
ccp4-python --version  # Must be 3.11.x
```

### Django Database Errors

Ensure settings module is configured:
```bash
export DJANGO_SETTINGS_MODULE=ccp4i2.config.test_settings
```

### fixture 'django_db_blocker' not found

This error means `pytest-django` is not installed. It should be installed automatically by `pip install -e .`, but can be installed explicitly:
```bash
ccp4-python -m pip install pytest-django pytest-xdist
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
   - [Quick Reference](../QUICK_REFERENCE.md)

3. **Frontend development** (if working on the UI):
   - [Frontend README](../../client/README.md) - Quick start
   - [Frontend Development Guide](../../client/FRONTEND_DEVELOPMENT.md) - Full documentation
