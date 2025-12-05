# i2run: CCP4i2 Task Runner Command

`i2run` is the command-line interface for configuring and running CCP4i2 tasks directly from the command line. It provides direct access to all CPluginScript parameters, supports complex data structures, and integrates with the Django project database.

## Table of Contents

- [Quick Start](#quick-start)
- [Basic Invocation](#basic-invocation)
- [Non-Task-Specific Arguments](#non-task-specific-arguments)
- [Parameter Setting Syntax](#parameter-setting-syntax)
  - [Simple Parameters](#simple-parameters)
  - [File Parameters](#file-parameters)
  - [Complex Data Objects](#complex-data-objects)
  - [List Parameters](#list-parameters)
- [fileUse Syntax](#fileuse-syntax)
- [Examples](#examples)
- [Debugging Tips](#debugging-tips)

---

## Quick Start

```bash
# Run a task with basic parameters
python manage.py i2run import_merged_mtz --project_name toxd --HKLIN /path/to/data.mtz

# Run with file sub-parameters
python manage.py i2run refmac5 --project_name toxd \
    --HKLIN fullPath=/path/to/reflections.mtz \
    --XYZIN fullPath=/path/to/model.pdb

# Configure only (don't execute)
python manage.py i2run refmac5 --project_name toxd --XYZIN fullPath=/model.pdb --i2run_configure
```

---

## Basic Invocation

```
python manage.py i2run <task_name> [--project_name <name>] [--param1 <value>] [--param2 <value>] ...
```

### Arguments

- `<task_name>` - (Required, positional) The name of the task/plugin to run (e.g., `refmac5`, `import_merged_mtz`, `servalcat_pipe`)
- All other arguments are dynamic and depend on the task's def.xml definition

---

## Non-Task-Specific Arguments

These arguments are common to all tasks and control i2run behavior:

| Argument | Description |
|----------|-------------|
| `--project_name <name>` | Project name to create job in (creates project if needed) |
| `--project_path <path>` | Custom project directory path |
| `--delay` | Delay execution (for GUI integration) |
| `--batch` | Run in batch mode (non-interactive) |
| `--i2run_configure` | Configure job parameters but don't execute |

### Example: Configure Only

The `--i2run_configure` flag is useful for setting up a job without running it:

```bash
python manage.py i2run refmac5 --project_name toxd \
    --XYZIN fullPath=/path/to/model.pdb \
    --HKLIN fullPath=/path/to/data.mtz \
    --i2run_configure
```

This creates the job in the database with all parameters set, ready to be run later via the GUI or API.

---

## Parameter Setting Syntax

### Simple Parameters

For simple parameters (CString, CInt, CFloat, CBoolean), use direct values:

```bash
# Integer parameter
--NCYCLES 10

# Float parameter
--RESOLUTION 2.5

# String parameter
--TITLE "Refinement round 1"

# Boolean parameter (use True/False or 1/0)
--USE_TWIN True
```

### File Parameters

File parameters (CDataFile types) support sub-attributes using the `key=value` syntax:

```bash
# Set file path
--XYZIN fullPath=/path/to/model.pdb

# Set annotation
--XYZIN annotation="Input model structure"

# Multiple sub-attributes for the same parameter
--XYZIN fullPath=/path/to/model.pdb --XYZIN annotation="My model"
```

#### Common File Sub-Attributes

| Attribute | Description |
|-----------|-------------|
| `fullPath=` | Full filesystem path to the file |
| `annotation=` | Human-readable description |
| `columnLabels=` | MTZ column label selections (for MTZ files) |
| `selection=` | Atom selection string (for PDB files) |
| `selection/text=` | Selection text within selection object |

### Complex Data Objects

For complex CData objects with nested structure, use slash-separated paths:

```bash
# Set nested attribute
--SELECTION selection/text="chain A and resname ALA"

# For MTZ column labels
--F_SIGF columnLabels/F=FP --F_SIGF columnLabels/SIGF=SIGFP
```

### List Parameters

For list parameters (CList types), you can add multiple items:

```bash
# Add multiple ensemble items with nested attributes
--ENSEMBLES sequence=MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH \
--ENSEMBLES nCopies=1 \
--ENSEMBLES pdbItemList/structure/fullPath=/path/to/search_model.pdb
```

For multi-item lists, specify the index explicitly or let i2run auto-append.

---

## Type-Based Parsing

i2run uses intelligent type-based parsing:

1. **Fundamental types** (CString, CInt, CFloat, CBoolean) - Values are taken literally, even if they contain `=`

   ```bash
   # SMILES strings work correctly - the = is literal, not key=value
   --SMILESIN "CN1CCC(=O)CC4"
   ```

2. **Composite types** (CDataFile, CData, CContainer) - Support `key=value` syntax for sub-attributes

   ```bash
   # File parameters support sub-attributes
   --HKLIN fullPath=/data.mtz columnLabels/F=FP
   ```

---

## fileUse Syntax

Reference files from previous jobs using the `fileUse:` syntax:

```
fileUse:<project_name>/<task_name>/<job_index>/<param_name>/<file_index>
```

### Components

| Component | Description |
|-----------|-------------|
| `project_name` | Project containing the source file |
| `task_name` | Task type (optional filter) |
| `job_index` | 0-based job index within matching jobs |
| `param_name` | Parameter name in the source job |
| `file_index` | 0-based index if multiple files (optional) |

### Examples

```bash
# Reference first HKLOUT from first import_merged_mtz job
--HKLIN "fileUse:toxd/import_merged_mtz/0/HKLOUT"

# Reference XYZOUT from second refmac5 job
--XYZIN "fileUse:toxd/refmac5/1/XYZOUT"
```

---

## Examples

### Import Merged MTZ Data

```bash
python manage.py i2run import_merged_mtz \
    --project_name toxd \
    --HKLIN fullPath=/data/toxd_free.mtz \
    --HKLIN annotation="Toxd native data"
```

### Run Refmac5 Refinement

```bash
python manage.py i2run refmac5 \
    --project_name toxd \
    --XYZIN fullPath=/models/toxd_initial.pdb \
    --HKLIN fullPath=/data/toxd_free.mtz \
    --NCYCLES 10 \
    --WEIGHT_TYPE AUTO
```

### Run Servalcat Pipeline

```bash
python manage.py i2run servalcat_pipe \
    --project_name toxd \
    --XYZIN fullPath=/models/toxd.pdb \
    --HKLIN fullPath=/data/toxd.mtz \
    --PROSMART_PROTEIN True \
    --NCYCLES 5
```

### Phaser Molecular Replacement

```bash
python manage.py i2run phaser_simple \
    --project_name mr_project \
    --HKLIN fullPath=/data/target.mtz \
    --XYZIN fullPath=/models/search_model.pdb \
    --NCOPIES 2 \
    --ID_RMS ID \
    --SEARCHSEQUENCEIDENTITY 0.5
```

### AceDRG Ligand

```bash
python manage.py i2run acedrgNew \
    --project_name ligand_project \
    --SMILESIN "CC(=O)Oc1ccccc1C(=O)O"
```

### Create Monomer Link

```bash
python manage.py i2run MakeLink \
    --project_name glyco_project \
    --MON_1_TYPE TLC \
    --RES_NAME_1_TLC NAG \
    --ATOM_NAME_1 C1 \
    --MON_2_TYPE TLC \
    --RES_NAME_2_TLC ASN \
    --ATOM_NAME_2 ND2
```

### Chain Jobs with fileUse

```bash
# First: import data
python manage.py i2run import_merged_mtz \
    --project_name workflow \
    --HKLIN fullPath=/data/native.mtz

# Second: import model
python manage.py i2run import_coords \
    --project_name workflow \
    --XYZIN fullPath=/models/initial.pdb

# Third: run refinement using outputs from previous jobs
python manage.py i2run refmac5 \
    --project_name workflow \
    --HKLIN "fileUse:workflow/import_merged_mtz/0/HKLOUT" \
    --XYZIN "fileUse:workflow/import_coords/0/XYZOUT"
```

---

## Parameter Discovery

To discover available parameters for a task, use the plugins CLI command:

```bash
# List all available parameters for refmac5
ccp4i2 plugins show refmac5

# Or create a job and examine the container
python manage.py i2run refmac5 --project_name test --i2run_configure
# This outputs the configured plugin XML showing all parameters
```

---

## Debugging Tips

### 1. Use --i2run_configure

Configure without running to verify parameters are set correctly:

```bash
python manage.py i2run task_name --project_name test --PARAM value --i2run_configure
```

### 2. Enable Debug Logging

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

### 3. Check Parameter Paths

The i2run system computes minimum unique paths for parameters. If a parameter name is ambiguous, use the full path:

```bash
# If XYZIN is ambiguous (appears in multiple containers)
--inputData.XYZIN fullPath=/path/to/file.pdb
```

### 4. Common Issues

| Issue | Solution |
|-------|----------|
| "Plugin not found" | Check task name spelling, regenerate plugin registry |
| "Parameter not found" | Use full path, check def.xml for correct parameter name |
| "fileUse not resolved" | Verify project/job exists, check job index |
| "Value parsing error" | Check quoting, ensure special characters are escaped |

---

## Architecture

The i2run system consists of:

1. **CCP4i2RunnerBase** - Common argument parsing and workflow
2. **CCP4i2RunnerDjango** - Django-specific implementation with ORM integration
3. **KeywordExtractor** - Extracts parameter metadata from plugin definitions
4. **ArgumentBuilder** - Builds argparse arguments with backward-compatible aliases
5. **PluginPopulator** - Populates plugin instances with parsed command-line arguments

---

## See Also

- [CLI Quick Reference](CLI_QUICK_REFERENCE.md) - Modern CLI commands
- [CLI Documentation](CLI.md) - Full CLI reference
- [QUICK_REFERENCE.md](../QUICK_REFERENCE.md) - Plugin development patterns
