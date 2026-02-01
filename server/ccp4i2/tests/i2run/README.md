# CCP4i2 i2run Testing Guide

This document describes the i2run testing infrastructure for CCP4i2, covering test execution, the command-line job runner mechanism, and how to write new tests.

## Quick Start

```bash
# From the server/ directory
./run_test.sh ccp4i2/tests/i2run/ -v

# Run a specific quick test
./run_test.sh ccp4i2/tests/i2run/test_acorn.py -v

# Run multiple quick tests in parallel
./run_test.sh ccp4i2/tests/i2run/ -n 4 -k "acorn or aimless or parrot"
```

## What are i2run Tests?

The i2run tests exercise CCP4i2's **command-line job execution** system. They test:

1. **Keyword Discovery** - Automatic extraction of parameters from `.def.xml` files
2. **Argument Parsing** - Translation of command-line arguments to plugin parameters
3. **Job Execution** - Running crystallographic tasks through the Django backend
4. **Output Validation** - Verifying job outputs (MTZ files, PDB files, XML reports)

Each test defines a complete CCP4i2 job on a single command line, runs it, and validates the results.

## Prerequisites

### 1. CCP4 Environment

The tests require the CCP4 suite with `ccp4-python`. Before running tests:

```bash
# Option A: Source CCP4 setup directly
source /path/to/ccp4-20251105/bin/ccp4.setup-sh

# Option B: Place CCP4 in a sibling directory of the project (auto-detected)
# ~/Developer/ccp4-20251105/  (sibling of ~/Developer/ccp4i2/)

# Option C: Set CCP4_ROOT in .env file (in server/ or project root)
echo "CCP4_ROOT=/path/to/ccp4-20251105" > .env
```

### 2. Python Packages

Install ccp4i2 with test dependencies:

```bash
ccp4-python -m pip install -e ".[full]"
```

Key packages required:
- `pytest` - Test framework
- `pytest-django` - Django integration
- `pytest-xdist` - Parallel test execution
- `gemmi` - Crystallographic file validation

### 3. Test Data

Tests use data from two sources:

**Local demo_data** (included in repository):
```
server/ccp4i2/demo_data/
├── gamma/           # Gamma-B crystallin data
├── beta_blip/       # Beta-lactamase/BLIP complex
├── baz2b/           # BAZ2B bromodomain
└── ...
```

**Downloaded fixtures** (fetched automatically):
- PDBe FASTA sequences
- PDB-REDO MTZ/mmCIF files
- RCSB mmCIF structures

Downloaded files are cached in temporary directories and cleaned up after tests.

## Running Tests

### Basic Usage

```bash
# Run all i2run tests (default when no path specified)
./run_test.sh

# Run all i2run tests explicitly
./run_test.sh ccp4i2/tests/i2run/

# Run a specific test file
./run_test.sh ccp4i2/tests/i2run/test_acorn.py

# Run a specific test function
./run_test.sh ccp4i2/tests/i2run/test_aimless.py::test_gamma
```

### Quick Tests

Some tests are notably fast (< 10 seconds):

```bash
# Fast tests for quick validation
./run_test.sh ccp4i2/tests/i2run/test_acorn.py -v
./run_test.sh ccp4i2/tests/i2run/test_aimless.py -v
./run_test.sh ccp4i2/tests/i2run/test_parrot.py -v
./run_test.sh ccp4i2/tests/i2run/test_freerflag.py -v
```

### Slow Tests

Some tests run full crystallographic pipelines and take longer:

| Test | Approximate Time | Notes |
|------|------------------|-------|
| `test_crank2.py` | 5-15 minutes | Full experimental phasing pipeline |
| `test_arcimboldo.py` | 10+ minutes | Ab initio phasing |
| `test_xia2.py` | 10+ minutes | Data processing pipeline |
| `test_arpwarp.py` | 5-10 minutes | Automated model building |
| `test_mrbump.py` | 5-10 minutes | Molecular replacement pipeline |
| `test_buster.py` | 5+ minutes | Refinement (requires license) |

To exclude slow tests:

```bash
./run_test.sh ccp4i2/tests/i2run/ \
    --ignore=ccp4i2/tests/i2run/test_crank2.py \
    --ignore=ccp4i2/tests/i2run/test_arcimboldo.py \
    --ignore=ccp4i2/tests/i2run/test_xia2.py \
    --ignore=ccp4i2/tests/i2run/test_mrbump.py \
    --ignore=ccp4i2/tests/i2run/test_arpwarp.py \
    --ignore=ccp4i2/tests/i2run/test_buster.py \
    -v
```

### Parallel Execution

Tests support parallel execution via pytest-xdist:

```bash
# Run with 4 parallel workers
./run_test.sh ccp4i2/tests/i2run/ -n 4

# Auto-detect CPU count
./run_test.sh ccp4i2/tests/i2run/ -n auto

# Exclude slow tests and run in parallel
./run_test.sh ccp4i2/tests/i2run/ -n 4 \
    --ignore=ccp4i2/tests/i2run/test_crank2.py \
    --ignore=ccp4i2/tests/i2run/test_arcimboldo.py
```

### Verbose Output

```bash
# Show test names as they run
./run_test.sh ccp4i2/tests/i2run/ -v

# Show print statements and program output
./run_test.sh ccp4i2/tests/i2run/ -v -s

# Short traceback on failures
./run_test.sh ccp4i2/tests/i2run/ --tb=short
```

## Test Mechanism

### The i2run Context Manager

Tests use the `i2run()` context manager from `utils.py`:

```python
from .utils import i2run, demoData

def test_acorn():
    args = ["acorn"]  # Task name
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]

    with i2run(args) as job:
        # job is a Path to the job directory (CCP4_JOBS/job_N)
        assert (job / "PHSOUT.mtz").exists()
        gemmi.read_mtz_file(str(job / "PHSOUT.mtz"))
```

The context manager:
1. Creates an isolated project directory with unique name
2. Invokes the Django `i2run` management command in-process
3. Waits for job completion
4. Checks `diagnostic.xml` for errors
5. Yields the job directory for output validation
6. Cleans up on success; preserves on failure

### Keyword Discovery from .def.xml

The i2run system automatically discovers parameters from plugin `.def.xml` files:

```
wrappers/acorn/acorn.def.xml
    └── Defines: F_SIGF, XYZIN, PHSOUT, etc.
         └── Parsed by KeywordExtractor
              └── Creates argparse arguments: --F_SIGF, --XYZIN, etc.
```

This allows defining any CCP4i2 job entirely on the command line:

```bash
# What the test does internally:
ccp4-python manage.py i2run acorn \
    --project_name test_acorn \
    --F_SIGF /path/to/data.mtz \
    --XYZIN /path/to/model.pdb
```

### Database Isolation

Each test gets its own isolated environment:

1. **Unique project directory**: `~/.cache/ccp4i2-tests/YYYYMMDD_HHMMSS_XXXX_{test_name}/`
2. **Isolated SQLite database**: `{project_dir}/project.sqlite`
3. **Job outputs**: `{project_dir}/CCP4_JOBS/job_1/`

This enables safe parallel execution without conflicts.

### Key Components

| Component | File | Purpose |
|-----------|------|---------|
| `KeywordExtractor` | `cli/i2run/i2run_components.py` | Extracts parameters from .def.xml |
| `ArgumentBuilder` | `cli/i2run/i2run_components.py` | Creates argparse arguments |
| `PluginPopulator` | `cli/i2run/i2run_components.py` | Populates plugin from parsed args |
| `CCP4i2RunnerDjango` | `cli/i2run/CCP4i2RunnerDjango.py` | Executes jobs with Django backend |
| `i2run` context manager | `tests/i2run/utils.py` | Test helper for running jobs |

## Writing New Tests

### Basic Test Structure

```python
import gemmi
import xml.etree.ElementTree as ET
from .utils import i2run, demoData

def test_my_task():
    """Test description."""
    args = ["task_name"]
    args += ["--INPUT_FILE", demoData("gamma", "some_file.mtz")]
    args += ["--PARAMETER", "value"]

    with i2run(args) as job:
        # Validate outputs exist
        assert (job / "OUTPUT.mtz").exists()

        # Validate file contents
        mtz = gemmi.read_mtz_file(str(job / "OUTPUT.mtz"))
        assert mtz.spacegroup.hm == "P 21 21 21"

        # Check program.xml for results
        tree = ET.parse(job / "program.xml")
        result = tree.find(".//SomeResult")
        assert float(result.text) > 0.5
```

### Using Downloaded Fixtures

For tests requiring external data:

```python
from .utils import download, i2run

def test_with_pdb_redo(cif8xfm, mtz8xfm):
    """Test using PDB-REDO data (fixtures defined in conftest.py)."""
    args = ["refmac"]
    args += ["--XYZIN", cif8xfm]
    args += ["--HKLIN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]

    with i2run(args) as job:
        assert (job / "XYZOUT.cif").exists()
```

Available session fixtures (defined in `conftest.py`):
- `cif7beq`, `mtz7beq` - PDB entry 7BEQ
- `cif8xfm`, `mtz8xfm`, `seq8xfm` - PDB entry 8XFM

### Complex Parameter Syntax

For parameters with subvalues:

```python
def test_with_column_labels():
    args = ["aimless"]
    # File with column selection
    args += ["--HKLIN", f"fullPath={mtz_path}", "columnLabels=/*/*/[I,SIGI]"]

    with i2run(args) as job:
        pass

def test_with_list_items():
    args = ["split_mtz"]
    # Multiple COLUMNGROUPLIST items
    args += ["--COLUMNGROUPLIST", "columnGroupType=Phs", "contentFlag=1"]
    args += ["--COLUMNGROUPLIST", "columnGroupType=FandPhi", "contentFlag=2"]

    with i2run(args) as job:
        pass
```

### Expecting Errors

For tests that should produce errors:

```python
def test_expected_failure():
    args = ["phaser_simple"]
    args += ["--HKLIN", mtz_path]
    args += ["--ENSEMBLEFILE", pdb_path]
    # Missing required parameters...

    with i2run(args, allow_errors=True) as job:
        # Job ran but may have errors in diagnostic.xml
        xml_path = job / "diagnostic.xml"
        errors = ET.parse(xml_path).findall(".//errorReport")
        assert len(errors) > 0
```

### Skip Decorators

```python
import pytest

@pytest.mark.skip(reason="Task not yet in registry")
def test_future_task():
    pass

@pytest.mark.skipif(
    not Path("/path/to/optional/program").exists(),
    reason="Optional program not installed"
)
def test_optional_feature():
    pass
```

## Debugging Failed Tests

### Preserved Directories

Failed tests preserve their project directories:

```bash
# List test directories
ls ~/.cache/ccp4i2-tests/

# Example output:
# 20240201_143052_a1b2_acorn_test_acorn/
# 20240201_143105_c3d4_aimless_test_gamma/

# Inspect a failed test
cd ~/.cache/ccp4i2-tests/20240201_143052_a1b2_acorn_test_acorn/
ls CCP4_JOBS/job_1/
# diagnostic.xml  - Error reports
# input_params.xml - Input parameters
# program.xml     - Program output/results
# log.txt         - Program log (if any)
```

### Key Files to Check

| File | Purpose |
|------|---------|
| `diagnostic.xml` | Error and warning reports |
| `input_params.xml` | Parameters passed to task |
| `program.xml` | Task output and results |
| `*.log` | Program stdout/stderr |

### Cleanup

```bash
# Remove all test projects
rm -rf ~/.cache/ccp4i2-tests/*

# Remove projects from a specific date
rm -rf ~/.cache/ccp4i2-tests/20240201_*

# Keep only recent (e.g., last 10)
ls -t ~/.cache/ccp4i2-tests/ | tail -n +11 | xargs -I {} rm -rf ~/.cache/ccp4i2-tests/{}
```

## Environment Variables

| Variable | Purpose |
|----------|---------|
| `CCP4I2_ROOT` | Root of ccp4i2 installation (server/ directory) |
| `DJANGO_SETTINGS_MODULE` | Django settings (default: `ccp4i2.config.test_settings`) |
| `CCP4I2_PROJECTS_DIR` | Directory for test projects |
| `CCP4` | CCP4 installation root (set by ccp4.setup-sh) |
| `CLIBD_MON` | CCP4 monomer library path |

## Test Coverage

### Current Test Files

| File | Tasks Tested |
|------|--------------|
| `test_acorn.py` | ACORN phase improvement |
| `test_acedrg.py` | AceDRG ligand dictionary |
| `test_aimless.py` | Aimless scaling |
| `test_arcimboldo.py` | ARCIMBOLDO ab initio phasing |
| `test_arpwarp.py` | ARP/wARP model building |
| `test_asu_contents.py` | ASU contents analysis |
| `test_auspex.py` | Auspex data analysis |
| `test_buster.py` | BUSTER refinement |
| `test_crank2.py` | Crank2 experimental phasing |
| `test_csymmatch.py` | CSYMMATCH symmetry matching |
| `test_dimple.py` | DIMPLE quick refinement |
| `test_editbfac.py` | Edit B-factors (AlphaFold models) |
| `test_find_waters.py` | Water finding |
| `test_freerflag.py` | FreeR flag assignment |
| `test_import_merged.py` | MTZ import |
| `test_lorestr.py` | LORESTR low-resolution refinement |
| `test_metalcoord.py` | Metal coordination analysis |
| `test_model_asu_check.py` | Model ASU checking |
| `test_modelcraft.py` | ModelCraft model building |
| `test_molrep.py` | MolRep molecular replacement |
| `test_mrbump.py` | MrBUMP MR pipeline |
| `test_mrparse.py` | MrParse sequence analysis |
| `test_parrot.py` | Parrot density modification |
| `test_phaser_*.py` | Phaser molecular replacement |
| `test_refmac.py` | Refmac5 refinement |
| `test_servalcat.py` | Servalcat refinement |
| `test_sheetbend.py` | Sheetbend model morphing |
| `test_split_mtz.py` | MTZ column extraction |
| `test_substitute_ligand.py` | Ligand substitution |
| `test_validate.py` | Structure validation |
| `test_xia2.py` | xia2 data processing |

---

## Architecture Notes

The i2run test infrastructure demonstrates CCP4i2's ability to:

1. **Programmatic Job Definition** - Define complete jobs via command-line arguments
2. **Automatic Parameter Discovery** - Extract parameters from XML definitions
3. **Django Integration** - Execute jobs with full database tracking
4. **Validation** - Verify crystallographic outputs with gemmi

This enables both automated testing and scriptable job execution for pipelines.
