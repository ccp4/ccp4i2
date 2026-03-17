# Writing i2run Tests for CCP4i2

This guide covers how to write integration tests for CCP4i2 task wrappers and
pipelines using the `i2run` test harness.  The infrastructure lives under
`server/ccp4i2/tests/i2run/` and is documented in CLAUDE.md under **Testing**.

## Quick start

```python
# server/ccp4i2/tests/i2run/test_mytask.py
from .utils import i2run, demoData

def test_mytask_basic():
    args = ["mytask"]
    args += ["--HKLIN", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]

    with i2run(args) as job:
        assert (job / "XYZOUT.cif").exists()
```

Run it:

```bash
cd server
./run_test.sh ccp4i2/tests/i2run/test_mytask.py -v
```

## How the harness works

`i2run(args)` is a context manager that:

1. Creates an isolated project directory with a unique timestamped name
2. Creates a per-test SQLite database inside that directory
3. Invokes the Django `i2run` management command **in-process**
4. Parameters are auto-discovered from the task's `.def.xml`
5. After the task completes, `diagnostic.xml` is checked for errors (severity >= 4)
6. The job directory (`CCP4_JOBS/job_1/`) is yielded as a `pathlib.Path`
7. On success the directory is cleaned up; on failure it is preserved for debugging

## Test data sources

### Local demo data

Use `demoData()` to reference files shipped with CCP4i2:

```python
from .utils import demoData

mtz = demoData("gamma", "merged_intensities_Xe.mtz")
pdb = demoData("gamma", "gamma_model.pdb")
```

### Remote fixtures (PDB, PDB-REDO, RCSB)

For tests that need real crystallographic data, download from public
repositories using session-scoped fixtures.  Downloads happen once per test
session and are cached.

**Step 1** — Add URL builders to your test or use existing ones from `urls.py`:

```python
# urls.py provides these builders:
from .urls import (
    pdbe_pdb,       # PDBe PDB format        — pdbe_pdb("1h1s")
    pdbe_mmcif,     # PDBe mmCIF             — pdbe_mmcif("7beq")
    pdbe_sfcif,     # PDBe structure factors  — pdbe_sfcif("7beq")
    pdbe_fasta,     # PDBe FASTA sequence     — pdbe_fasta("8xfm")
    redo_mtz,       # PDB-REDO MTZ            — redo_mtz("1h1s")
    redo_cif,       # PDB-REDO refined mmCIF  — redo_cif("8xfm")
    rcsb_pdb,       # RCSB PDB format         — rcsb_pdb("1h1s")
    rcsb_mmcif,     # RCSB mmCIF              — rcsb_mmcif("7beq")
    rcsb_ligand_cif,# RCSB ligand definition  — rcsb_ligand_cif("ATP")
    rcsb_ligand_sdf,# RCSB ligand 3D coords   — rcsb_ligand_sdf("ATP")
)
```

**Step 2** — Add session-scoped fixtures in `conftest.py`:

```python
from .urls import redo_mtz, pdbe_pdb
from .utils import download

@fixture(scope="session")
def pdb1h1s():
    with download(pdbe_pdb("1h1s")) as path:
        yield path

@fixture(scope="session")
def mtz1h1s():
    with download(redo_mtz("1h1s")) as path:
        yield path
```

Session scope means the file is downloaded once and shared across all tests in
the run.  The `download()` context manager handles temp-file cleanup.

**Step 3** — Use the fixtures in your test:

```python
def test_mytask_with_remote_data(pdb1h1s, mtz1h1s):
    args = ["mytask"]
    args += ["--XYZIN", pdb1h1s]
    args += ["--F_SIGF", f"fullPath={mtz1h1s}", "columnLabels=/*/*/[FP,SIGFP]"]

    with i2run(args) as job:
        assert (job / "XYZOUT.cif").exists()
```

**Tip**: PDB-REDO MTZ files (`redo_mtz`) are preferred over PDBe structure
factor CIFs because they are already in MTZ format with standard column names
(FP, SIGFP, FREE).  This avoids conversion steps in your test.

## Parameter syntax

Parameters follow the `.def.xml` keyword names with `--` prefix.

### Simple parameters

```python
args += ["--RESOLUTION_HIGH", "2.0"]
args += ["--RUNREFMAC", "True"]
args += ["--NCYCLES", "10"]
```

### File parameters with column selection

When an MTZ contains many columns, use `columnLabels` to select specific ones:

```python
args += ["--F_SIGF", f"fullPath={mtz_path}", "columnLabels=/*/*/[FP,SIGFP]"]
args += ["--FREERFLAG", f"fullPath={mtz_path}", "columnLabels=/*/*/[FREE]"]
```

The `columnLabels` syntax is `/*/*/[COL1,COL2,...]` — this triggers MTZ
splitting to extract only the named columns.

### File parameters with sub-keywords

Some file types accept additional sub-keywords:

```python
args += ["--ASUIN", f"seqFile={seq_path}"]
```

### List parameters (CList items)

To add multiple items to a CList, repeat the keyword:

```python
args += ["--SELECTIONS", "text=A/"]
args += ["--SELECTIONS", "text=B/"]
args += ["--SELECTIONS", "text=C/"]
```

Each `--SELECTIONS` adds a new item to the list.  Sub-keywords like `text=`
set attributes on the list item.

### Compound list items

For list items with multiple sub-values:

```python
args += ["--ENSEMBLES", "pdbFile=/path/to/model.pdb", "identity=0.9"]
args += ["--ENSEMBLES", "pdbFile=/path/to/model2.pdb", "identity=0.8"]
```

## Common assertion patterns

### Check output files exist

```python
with i2run(args) as job:
    assert (job / "XYZOUT.cif").exists(), "No output coordinates"
    assert (job / "HKLOUT.mtz").exists(), "No output reflections"
```

### Validate crystallographic files with gemmi

```python
import gemmi

with i2run(args) as job:
    st = gemmi.read_structure(str(job / "XYZOUT.cif"))
    assert len(st) > 0, "Empty structure"

    mtz = gemmi.read_mtz_file(str(job / "HKLOUT.mtz"))
    assert mtz.nreflections > 0
```

### Parse program.xml for results

```python
import xml.etree.ElementTree as ET

with i2run(args) as job:
    tree = ET.parse(job / "program.xml")
    rfree = float(tree.find(".//r_free").text)
    assert rfree < 0.35, f"R-free too high: {rfree}"
```

### Check for expected errors

```python
with i2run(args, allow_errors=True) as job:
    xml = ET.parse(job / "diagnostic.xml")
    errors = xml.findall(".//errorReport")
    assert len(errors) > 0, "Expected an error report"
```

## Controlling test order

Some tests (e.g., Phaser) are sensitive to module-level state pollution from
packages like RDKit.  Use `pytest.mark.order` to run them first:

```python
import pytest

@pytest.mark.order("first")
def test_phaser_runs_first(mtz8xfm, cif8xfm):
    ...
```

## Running tests

```bash
cd server

# Single test file
./run_test.sh ccp4i2/tests/i2run/test_mytask.py -v

# Single test function
./run_test.sh ccp4i2/tests/i2run/test_mytask.py::test_mytask_basic -v

# All i2run tests
./run_test.sh ccp4i2/tests/i2run/ -v

# Parallel (recommended — most tasks are single-threaded)
./run_test.sh ccp4i2/tests/i2run/ -n 8

# With stdout visible
./run_test.sh ccp4i2/tests/i2run/ -v -s
```

## Debugging failed tests

Failed tests preserve their project directory under `~/.cache/ccp4i2-tests/`.
Key files to inspect:

| File | Contents |
|------|----------|
| `diagnostic.xml` | Error reports and warnings |
| `program.xml` | Task results and statistics |
| `input_params.xml` | Parameters as received by the task |
| `log.txt` | Full program log output |
| `job_N/` | Sub-job directories (for pipelines) |

```bash
# List preserved test directories
ls -la ~/.cache/ccp4i2-tests/

# Clean up all
rm -rf ~/.cache/ccp4i2-tests/*

# Clean up a specific day
rm -rf ~/.cache/ccp4i2-tests/20260317_*
```

## Checklist for a new test

1. Identify the task name (as it appears in `.def.xml` / plugin registry)
2. Find suitable test data — demo_data for simple cases, remote fixtures for
   real-world data
3. Add any new session-scoped download fixtures to `conftest.py`
4. Write the test using `i2run(args)` context manager
5. Assert on output files, parsed XML results, or gemmi validation
6. Run locally: `./run_test.sh ccp4i2/tests/i2run/test_mytask.py -v`
7. Run full suite to check for regressions: `./run_test.sh ccp4i2/tests/i2run/ -v`
