# i2run — Reference Guide

`i2run` is the command-line (and programmatic) entry point for running CCP4i2
tasks in Django mode.  It handles project creation, container population,
job execution and result persistence.

---

## Command-line invocation

```bash
source /path/to/ccp4-YYYYMMDD/bin/ccp4.setup-sh
export CCP4I2_BACKEND=django
export DJANGO_SETTINGS_MODULE=ccp4x.config.settings

cd server
ccp4-python manage.py i2run <task_name> [--project_name NAME] [--PARAM value ...]
```

### Positional argument

| Argument | Description |
|---|---|
| `task_name` | Plugin TASKNAME, e.g. `aimless_pipe`, `scaleit`, `refmac5` |

### Fixed optional arguments

| Argument | Description |
|---|---|
| `--project_name NAME` | Project directory name (auto-generated if omitted) |
| `--project_path PATH` | Explicit absolute path for the project directory |
| `--delay` | Queue the job for later execution (do not run immediately) |
| `--batch` | Batch mode (suppress interactive prompts) |

---

## Specifying task parameters

Every task-specific parameter is introduced with `--PARAM_NAME`, where the
name matches the `id` attribute in the task's `*.def.xml` container.

### Single-value parameters

```bash
--RESOLUTION_MAX 2.5
--SPACEGROUP "P 21 21 21"
```

### File parameters

File objects accept one or more `key=value` tokens after the flag.  All tokens
are collected together by the parser (they are not separate `--` flags).

| Key | Meaning |
|---|---|
| `file=<path>` | Shorthand path key — recognised by some file types (e.g. unmerged MTZ). **Use `fullPath=` when in doubt.** |
| `fullPath=<path>` | Absolute filesystem path — works for all `CDataFile`-derived types |
| `columnLabels=/*/*/[F,SIGF]` | MTZ column selection (triggers automatic column extraction) |
| `annotation=<text>` | Human-readable label shown in the GUI |

Examples:

```bash
# Minimal — just a path
--HKLIN file=/data/native.mtz

# With column selection
--HKLIN fullPath=/data/merged.mtz columnLabels=/*/*/[F,SIGF]

# With annotation
--FREERFLAG file=/data/freer.mtz annotation="Free-R from aimless"
```

### List parameters (CList / `CMiniMtzDataFileList` etc.)

A list parameter accepts **one item per flag invocation**.  Repeat `--PARAM`
once for each item you want to add to the list.

```bash
# Two files in MERGEDFILES (e.g. for scaleit)
--MERGEDFILES file=/data/native.mtz \
--MERGEDFILES file=/data/derivative.mtz

# Three unmerged sweeps for aimless_pipe
--UNMERGEDFILES file=/data/sweep1.mtz \
--UNMERGEDFILES file=/data/sweep2.mtz \
--UNMERGEDFILES file=/data/sweep3.mtz
```

Each invocation appends one element to the list; there is no indexed notation
for list-type top-level parameters (but see Nested access below).

### Nested / sub-object access

Use `/` to navigate into sub-objects:

```bash
--UNMERGEDFILES file=/data/sweep1.mtz crystalName=mycrystal dataset=DS1
```

Square-bracket indexing is used for sub-CLists *inside* a parameter (e.g.
`columnList` entries in `splitMtz`):

```bash
--COLUMNGROUPLIST \
  columnGroupType=Phs contentFlag=1 dataset=ds1 selected=True \
  "columnList[0]/columnLabel=HLA" "columnList[0]/columnType=A" \
  "columnList[1]/columnLabel=HLB" "columnList[1]/columnType=A"
```

---

## Using i2run in tests

The `i2run()` context manager in `tests/i2run/utils.py` wraps the management
command for in-process use from pytest:

```python
from ccp4i2.tests.i2run.utils import demoData, i2run

def test_my_task():
    mtz = demoData("gamma", "merged_intensities_native.mtz")
    args = ["my_task", "--HKLIN", f"file={mtz}"]
    with i2run(args) as job:
        # job is a pathlib.Path to the CCP4_JOBS/job_N directory
        mtz_out = job / "HKLOUT.mtz"
        assert mtz_out.exists()
```

### `demoData(*paths)`

Constructs a path rooted at `<ccp4i2_root>/demo_data/`:

```python
demoData("gamma", "gamma_native.mtz")        # unmerged native
demoData("gamma", "merged_intensities_native.mtz")  # merged I+/I- native
demoData("gamma", "merged_intensities_Xe.mtz")      # merged I+/I- Xe derivative
demoData("gamma", "freeR.mtz")               # free-R flags
demoData("gamma", "initial_phases.mtz")      # HL phase coefficients
demoData("mdm2",  "mdm2_unmerged.mtz")       # unmerged mdm2
```

### What `i2run()` does

1. Generates a unique project name (based on test function name + timestamp).
2. Sets `sys.argv` and calls the Django `i2run` management command in-process.
3. Finds the created job directory (`CCP4_JOBS/job_N/`).
4. Parses `diagnostic.xml` and fails the test if any severity-4 errors are
   present (unless `allow_errors=True`).
5. Yields the `Path` to the job directory for assertions.
6. Cleans up the project directory on success (preserves it on failure).

### Asserting on results

```python
import xml.etree.ElementTree as ET
from gemmi import read_mtz_file

with i2run(args) as job:
    # Verify an MTZ output is readable
    read_mtz_file(str(job / "HKLOUT.mtz"))

    # Parse program.xml for task-specific results
    tree = ET.parse(job / "program.xml")
    spacegroup = tree.find(".//Result/Dataset/SpacegroupName").text
    assert spacegroup == "P 21 21 21"
```

Common output files present in every job directory:

| File | Description |
|---|---|
| `program.xml` | Task-specific result XML written by `processOutputFiles()` |
| `diagnostic.xml` | Error/warning report from the CCP4i2 framework |
| `log.txt` (or `*.log`) | Program log file |
| `*.mtz` | Output reflection data files |
| `*.pdb` / `*.cif` | Output coordinate files |

---

## Common patterns by task type

### Scaling / merging (aimless_pipe)

```python
mtz = demoData("gamma", "gamma_native.mtz")
args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}"]
```

### Scaling with derivatives (scaleit)

```python
native = demoData("gamma", "merged_intensities_native.mtz")
xe     = demoData("gamma", "merged_intensities_Xe.mtz")
args = ["scaleit",
        "--MERGEDFILES", f"fullPath={native}",
        "--MERGEDFILES", f"fullPath={xe}"]
```

### Import merged MTZ with free-R

```python
args = ["import_merged"]
args += ["--HKLIN",    demoData("gamma", "merged_intensities_Xe.mtz")]
args += ["--FREERFLAG", demoData("gamma", "freeR.mtz")]
```

### Phasing / refinement (needs coordinates + data)

```python
args = ["refmac5"]
args += ["--XYZIN",  demoData("gamma", "gamma_model.pdb")]
args += ["--F_SIGF", f"file={demoData('gamma', 'mergedForMakingCif.mtz')}",
                     "columnLabels=/*/*/[F_SIGF_F,F_SIGF_SIGF]"]
args += ["--FREERFLAG", f"file={demoData('gamma', 'freeR.mtz')}"]
```

---

## How the argument parser works (internals)

The i2run management command passes `sys.argv[2:]` directly to
`CCP4i2RunnerDjango.parseArgs()`.  The runner:

1. Instantiates the plugin class to discover its container tree.
2. Registers one `argparse` argument per leaf parameter, computing the minimum
   unambiguous path for each.
3. All arguments use `nargs="+"` (each flag collects one or more tokens).
4. `CList`-type parameters additionally use `action="append"` so the flag can
   be repeated.
5. `PluginPopulator.populate()` walks the parsed namespace and sets values on
   the container using `key=value` token parsing and `/`-separated path
   navigation.

This means:
- A `CList` grows by repeating its flag — one item per invocation.
- Sub-attributes of a file/object are extra tokens on the *same* invocation,
  not separate flags.
- There is no `--PARAM[0]` index syntax for top-level list parameters.
