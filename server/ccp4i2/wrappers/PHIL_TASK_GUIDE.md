# Implementing a PHIL-Based Task in CCP4i2

This guide covers how to wrap a PHIL-based tool (Phenix, PhaserTNG, DIALS, etc.) as a CCP4i2 task using the `PhilPluginScript` infrastructure. It uses `phasertng_picard` as the reference implementation.

## Overview

Traditional CCP4i2 tasks define **all** their parameters in `.def.xml` files — hundreds of lines of hand-written XML specifying every input, output, and control parameter. For PHIL-based tools that already carry a complete `master_phil` with types, defaults, help text, and validation, this is redundant and error-prone.

`PhilPluginScript` solves this by:

- **Importing** the tool's `master_phil` at runtime and converting it to CData objects automatically
- Keeping CCP4i2's **rich file types** (CMtzDataFile, CPdbDataFile, etc.) for file I/O in `inputData` / `outputData`
- Using **shims** to bridge rich file types → PHIL path=value pairs at execution time
- Assembling a validated `working.phil` via `master_phil.fetch()` instead of concatenating strings

The result: a few hundred lines of Python replaces thousands of lines of XML, and the parameter UI stays in sync with upstream tool updates automatically.

## Architecture

```
┌──────────────────────────────────────────────────────────┐
│                   PhilPluginScript                       │
│                                                          │
│  ┌────────────┐  ┌────────────────┐  ┌───────────────┐  │
│  │ inputData  │  │ controlParams  │  │  outputData   │  │
│  │ (def.xml)  │  │ (from PHIL)    │  │  (def.xml)    │  │
│  │            │  │                │  │               │  │
│  │ HKLIN      │  │ picard.control │  │ XYZOUT        │  │
│  │ XYZIN      │  │ phasertng.*    │  │ HKLOUT        │  │
│  │ ASUIN      │  │ (~400 params)  │  │               │  │
│  │ DICT       │  │                │  │               │  │
│  └─────┬──────┘  └───────┬────────┘  └───────┬───────┘  │
│        │                 │                    │          │
│        ▼                 ▼                    │          │
│   ┌─────────┐    ┌──────────────┐             │          │
│   │  Shims  │───>│ working.phil │             │          │
│   └─────────┘    └──────┬───────┘             │          │
│                         │                     │          │
│                         ▼                     │          │
│                  ┌──────────────┐              │          │
│                  │   Execute    │──────────────┘          │
│                  │  phasertng   │  processOutputFiles()   │
│                  └──────────────┘                         │
└──────────────────────────────────────────────────────────┘
```

### The Three Layers

| Layer | Source | Purpose |
|-------|--------|---------|
| `inputData` / `outputData` | `.def.xml` | Rich file types with CCP4i2 features (file browser, validation, DB registration, column selection) |
| `controlParameters` | `master_phil` at runtime | All non-file parameters — auto-generated from the tool's PHIL tree |
| Shims | Python code | Bridge layer: convert rich inputData files → PHIL `name=value` pairs at execution time |

## Step-by-Step Implementation

### 1. Create the Directory Structure

```
wrappers/
  my_tool/
    script/
      my_tool.py          # Plugin script (PhilPluginScript subclass)
      my_tool.def.xml     # inputData + outputData only
```

### 2. Write the `.def.xml`

Define **only** `inputData` and `outputData`. Leave `controlParameters` empty — it gets populated at runtime from the PHIL tree.

```xml
<?xml version='1.0' encoding='UTF-8'?>
<ccp4i2>
  <ccp4i2_header>
    <function>DEF</function>
    <pluginName>my_tool</pluginName>
  </ccp4i2_header>
  <ccp4i2_body id="my_tool">

    <container id="inputData">
      <!-- Rich file types with CCP4i2 features -->
      <content id="HKLIN">
        <className>CMtzDataFile</className>
        <qualifiers>
          <guiLabel>Reflection data</guiLabel>
          <allowUndefined>False</allowUndefined>
          <mustExist>True</mustExist>
          <toolTip>Input MTZ file</toolTip>
        </qualifiers>
      </content>
      <content id="XYZIN">
        <className>CPdbDataFile</className>
        <qualifiers>
          <guiLabel>Model coordinates</guiLabel>
          <allowUndefined>True</allowUndefined>
          <mustExist>True</mustExist>
        </qualifiers>
      </content>
    </container>

    <container id="outputData">
      <content id="XYZOUT">
        <className>CPdbDataFile</className>
        <qualifiers>
          <guiLabel>Output coordinates</guiLabel>
          <default><subType>1</subType></default>
        </qualifiers>
      </content>
      <content id="HKLOUT">
        <className>CMtzDataFile</className>
        <qualifiers>
          <guiLabel>Output reflections</guiLabel>
        </qualifiers>
      </content>
    </container>

    <!-- Populated at runtime from master_phil -->
    <container id="controlParameters">
    </container>

  </ccp4i2_body>
</ccp4i2>
```

### 3. Write the Plugin Script

Subclass `PhilPluginScript` and implement four methods:

```python
import os
import shutil
import logging

from ccp4i2.core.PhilPluginScript import PhilPluginScript
from ccp4i2.utils.phil_shims import MtzFileShim, PdbFileShim

logger = logging.getLogger(__name__)


class my_tool(PhilPluginScript):
    TASKNAME = "my_tool"
    TASKCOMMAND = "my-tool"

    # PHIL scopes to hide from the GUI (handled by shims or CCP4i2)
    PHIL_EXCLUDE_SCOPES = [
        "input.hklin",          # Handled by MtzFileShim
        "input.xyzin",          # Handled by PdbFileShim
        "output",               # Managed by CCP4i2
    ]

    def get_master_phil(self):
        """Return the tool's master PHIL scope."""
        from my_tool import master_phil
        return master_phil

    def get_shim_definitions(self):
        """Map CCP4i2 rich file types to PHIL parameters."""
        return [
            MtzFileShim("HKLIN", "input.hklin"),
            PdbFileShim("XYZIN", "input.xyzin"),
        ]

    def get_command_target(self):
        """Return the executable name or path."""
        cmd = shutil.which("my-tool")
        return cmd or "my-tool"

    def _elevate_to_job_dir(self, source_path):
        """Copy a file into the job root if it isn't there already."""
        work_dir = str(self.getWorkDirectory())
        dest = os.path.join(work_dir, os.path.basename(source_path))
        if os.path.normpath(source_path) != os.path.normpath(dest):
            shutil.copy2(source_path, dest)
        return dest

    def processOutputFiles(self):
        """Find output files and register them as outputData.

        If the tool writes to a subdirectory, copy files into the job
        root first — glean_job_files() only registers files there.
        """
        import glob

        work_dir = str(self.getWorkDirectory())
        sub_dir = os.path.join(work_dir, "my_tool_output")
        out = self.container.outputData

        # Search subdirectory first, fall back to job root
        search_dirs = [sub_dir, work_dir] if os.path.isdir(sub_dir) else [work_dir]

        for d in search_dirs:
            pdbs = glob.glob(os.path.join(d, "*.pdb"))
            if pdbs:
                dest = self._elevate_to_job_dir(pdbs[0])
                out.XYZOUT.setFullPath(dest)
                break

        for d in search_dirs:
            mtzs = glob.glob(os.path.join(d, "*.mtz"))
            if mtzs:
                dest = self._elevate_to_job_dir(mtzs[0])
                out.HKLOUT.setFullPath(dest)
                break

        return PhilPluginScript.SUCCEEDED
```

### 4. Register the Task

Add an entry in `core/tasks.py`:

```python
"my_tool": Task(
    title="My Tool - Description",
    description="What the tool does.",
    shortTitle="My Tool",
    pluginPath="ccp4i2.wrappers.my_tool.script.my_tool:my_tool",
    defXmlPath="wrappers/my_tool/script/my_tool.def.xml",
),
```

And add `"my_tool"` to the appropriate category list in the `TASK_CATEGORIES` dict.

## The Subclass Contract

`PhilPluginScript` requires four methods. Here's what each does and when it's called:

### `get_master_phil() → libtbx.phil.scope`

**When called:** During `__init__` (via `_merge_phil_parameters()`) and during `build_working_phil()`.

**What to return:** The tool's complete master PHIL scope object.

**Notes:**
- Some tools use custom PHIL types (e.g., phasertng has `filesystem`, `mtzcol`, `scatterer`). These custom converters must be registered before parsing. If your tool uses `iotbx.cli_parser.CCTBXParser`, instantiate it to register converters automatically:

```python
def get_master_phil(self):
    from phasertng.programs import picard
    from iotbx.cli_parser import CCTBXParser
    parser = CCTBXParser(
        program_class=picard.Program,
        logger=None,
        parse_phil=False,
    )
    return parser.master_phil
```

- If the tool isn't installed, return `None` gracefully — the plugin will still instantiate (without PHIL parameters).

### `get_shim_definitions() → list[PhilShim]`

**When called:** During `build_working_phil()`, at execution time.

**What to return:** A list of shim instances that convert CCP4i2 rich file types into PHIL name=value pairs.

**Available shims:**

| Shim | Input Type | What it does |
|------|-----------|--------------|
| `MtzFileShim` | `CMtzDataFile` | Extracts the file path |
| `PdbFileShim` | `CPdbDataFile` | Extracts the file path; supports multiple PHIL targets |
| `PdbFileListShim` | `CList` of `CPdbDataFile` | Iterates the list, emitting one `(phil_path, value)` per item — for PHIL params with `.multiple = True` |
| `AsuContentShim` | `CAsuDataFile` | Loads ASU XML → writes FASTA file → returns path |
| `DictFileShim` | `CDictDataFile` | Extracts the CIF dictionary file path |
| `DagFileShim` | `CPhaserTngDagFile` | Extracts the DAG file path |

**Writing a custom shim:** Subclass `PhilShim` and implement `convert(container, work_directory) → list[(phil_path, value)]`.

### `get_command_target() → str`

**When called:** By the default `makeCommandAndScript()`.

**What to return:** The executable name or path (e.g., `"phasertng.picard"` or the result of `shutil.which()`).

### `processOutputFiles() → int`

**When called:** After the external program completes successfully.

**What to do:** Find output files on disk, copy them into the job root directory if necessary, and call `setFullPath()` on the appropriate `outputData` entries. The automatic gleaning system (`glean_job_files()`) then registers any set+existing files in the database.

**Critical:** `glean_job_files()` only registers files that live directly in the job directory (`getWorkDirectory()`). If the wrapped program writes output to a subdirectory (common for tools like PhaserTNG, ModelCraft, SIMBAD, etc.), you **must** copy those files into the job root before calling `setFullPath()`. Use `shutil.copy2()` to preserve metadata:

```python
def _elevate_to_job_dir(self, source_path):
    """Copy a file into the job root if it isn't there already."""
    work_dir = str(self.getWorkDirectory())
    dest = os.path.join(work_dir, os.path.basename(source_path))
    if os.path.normpath(source_path) != os.path.normpath(dest):
        shutil.copy2(source_path, dest)
    return dest
```

**Important:** `setFullPath()` overrides whatever `set_output_file_names()` configured pre-execution — no bespoke logic needed there.

## How Phil2CData Works

The `Phil2CData` converter (in `utils/phil_to_cdata.py`) walks the PHIL tree and creates a mirroring CData hierarchy:

### Type Mapping

| PHIL type | CData type | Notes |
|-----------|-----------|-------|
| `str` | `CString` | |
| `int` | `CInt` | `value_min`/`value_max` → `min`/`max` qualifiers |
| `float` | `CFloat` | `value_min`/`value_max` → `min`/`max` qualifiers |
| `bool` | `CBoolean` | Unless the default is `Auto` (see ternary below) |
| `choice` | `CString` | With `enumerators`, `onlyEnumerators`, `menuText` qualifiers |
| `path` | `CString` | Rich file types go in inputData, not here |
| `bool` (ternary) | `CString` | Default is `Auto`/`None` → enumerators `[True, False, Auto]` |
| Custom types | `CString` | `floats`, `ints`, `filesystem`, `mtzcol`, `scatterer`, etc. |

### Qualifier Mapping

| PHIL attribute | CData qualifier | Used for |
|---------------|----------------|----------|
| `short_caption` | `guiLabel` | Display label in the UI |
| `help` | `toolTip` | Tooltip text |
| `expert_level` | `expertLevel` | Parameter visibility level |
| `caption` | `menuText` | Dropdown labels for choice types |
| `value_min` | `min` | Minimum value constraint |
| `value_max` | `max` | Maximum value constraint |

### Name Mangling

PHIL uses dotted paths (`picard.control.use_shortcuts`). CData uses Python attribute names. The converter replaces dots with double underscores:

```
picard.control.use_shortcuts  →  picard__control__use_shortcuts
```

Each leaf also stores its original PHIL path in a `philPath` qualifier for reverse mapping at execution time.

### Scope → CContainer

Each PHIL scope becomes a nested `CContainer`:

```
picard {           →  CContainer("picard")
  control {        →    CContainer("picard__control")
    use_jog = F    →      CBoolean("picard__control__use_jog")
  }
}
```

### Default State Tracking

All values from `master_phil` are set with `ValueState.DEFAULT`. This means `obj.isSet(allowDefault=False)` returns `False` for defaults. When the user changes a value in the UI, the state becomes `EXPLICITLY_SET`. At execution time, `extract_phil_parameters()` only collects explicitly set values — defaults are already in the master_phil.

## Execution Flow

Here's the complete sequence from job launch to DB registration:

```
1. Plugin.__init__()
   ├─ super().__init__() loads .def.xml → creates inputData/outputData
   ├─ _merge_phil_parameters() imports master_phil → Phil2CData → controlParameters
   └─ loadDataFromXml() overlays any saved user params

2. processInputFiles()
   └─ Standard CPluginScript input validation

3. makeCommandAndScript()   ← can override or use default
   ├─ build_working_phil()
   │   ├─ extract_phil_parameters() → user-changed params only
   │   ├─ shims.convert() → rich file types → PHIL values
   │   ├─ master_phil.fetch(user_phil) → validated, complete PHIL
   │   └─ write working.phil to job directory
   └─ Construct command line: [executable, working.phil]

4. startProcess()
   └─ Execute external command

5. processOutputFiles()     ← you implement this
   ├─ Find output files on disk (may be in subdirectories)
   ├─ Copy files into job root if needed (glean_job_files only looks there)
   └─ Call out.XYZOUT.setFullPath(), out.HKLOUT.setFullPath(), etc.

6. set_output_file_names()
   └─ Configures project-relative paths (automatic, no custom code needed)

7. glean_job_files()        ← automatic
   ├─ Iterates outputData.find_all_files()
   ├─ Checks isSet() and exists() for each
   └─ Calls register_output_file() → creates File + FileUse DB records
```

## PHIL Scope Exclusion

The `PHIL_EXCLUDE_SCOPES` list tells Phil2CData to skip certain PHIL scopes entirely. This is critical for scopes that:

1. **Are handled by shims** — file I/O paths that come from CCP4i2's rich file types. Without exclusion, the user would see both a CCP4i2 file browser (inputData) AND a raw text field (controlParameters) for the same parameter.

2. **Are managed by CCP4i2** — output directory, file naming, etc. that CCP4i2 controls via its own mechanisms.

```python
PHIL_EXCLUDE_SCOPES = [
    "picard.hklin",      # → MtzFileShim handles this
    "picard.seqin",      # → AsuContentShim handles this
    "picard.xyzin",      # → PdbFileShim handles this
    "output",            # → CCP4i2 manages output paths
]
```

Exclusion works hierarchically: excluding `"picard"` would exclude everything under it.

## Testing

Three levels of testing are appropriate:

### Unit Tests (no external tool dependency)

Use synthetic PHIL via `libtbx.phil.parse()` to test the infrastructure:

```python
# test_phil_to_cdata.py — tests Phil2CData converter
# test_phil_shims.py — tests shim conversion
# test_phil_plugin_script.py — tests PhilPluginScript with mock subclass

class MockPhilPlugin(PhilPluginScript):
    TASKNAME = "mock_phil_plugin"
    TASKCOMMAND = "echo"

    def get_master_phil(self):
        return parse("""
            refinement {
              resolution = 2.0
                .type = float
              cycles = 10
                .type = int
            }
        """)

    def get_shim_definitions(self):
        return []

    def get_command_target(self):
        return "echo"
```

### Integration Tests (requires the external tool)

Gate with `pytest.mark.skipif`:

```python
try:
    from my_tool import master_phil
    _has_my_tool = True
except ImportError:
    _has_my_tool = False

pytestmark = pytest.mark.skipif(
    not _has_my_tool,
    reason="my_tool not installed"
)

def test_instantiation():
    cls = get_plugin_class("my_tool")
    plugin = cls()
    # Verify PHIL parameters were merged
    ...

def test_i2run():
    args = ["my_tool", "--HKLIN", demoData("gamma", "data.mtz"), ...]
    with i2run(args) as job:
        assert (job / "working.phil").exists()
        ...
```

### Running Tests

```bash
# Unit tests (always available)
./run_test.sh tests/test_phil_to_cdata.py -v
./run_test.sh tests/test_phil_shims.py -v

# Integration tests (requires the tool)
./run_test.sh tests/i2run/test_my_tool.py -v

# All PHIL tests
./run_test.sh -k "phil" -v
```

## Inspecting a Tool's master_phil

Before writing a wrapper, inspect the tool's PHIL to understand parameter names, types, and whether inputs accept multiple values. Use `ccp4-python` with the CCP4 environment sourced.

### Quick dump of the full PHIL tree

```bash
source ../ccp4-20251105/bin/ccp4.setup-sh

# For tools with a simple master_phil attribute:
ccp4-python -c "
from my_tool import master_phil
master_phil.show()
"

# For tools that need custom PHIL types registered (e.g., PhaserTNG):
ccp4-python -c "
from phasertng.programs import picard
from iotbx.cli_parser import CCTBXParser
parser = CCTBXParser(program_class=picard.Program, logger=None, parse_phil=False)
parser.master_phil.show()
"
```

### Checking whether a parameter is multiple (repeatable)

PHIL parameters with `.multiple = True` accept repeated values. This is critical for deciding whether the CCP4i2 input should be a single file or a `CList`:

```bash
ccp4-python -c "
from phasertng.programs import picard
from iotbx.cli_parser import CCTBXParser
parser = CCTBXParser(program_class=picard.Program, logger=None, parse_phil=False)

# Walk the PHIL tree looking for a specific parameter
for obj in parser.master_phil.all_definitions():
    if 'xyzin' in obj.path:
        print(f'{obj.path}  type={obj.object.type}  multiple={obj.object.multiple}')
"
```

Output like `picard.xyzin  type=path  multiple=True` tells you to use `CList` + `PdbFileListShim` rather than a single `CPdbDataFile` + `PdbFileShim`.

### Inspecting scopes and their attributes

```bash
ccp4-python -c "
from phasertng.programs import picard
from iotbx.cli_parser import CCTBXParser
parser = CCTBXParser(program_class=picard.Program, logger=None, parse_phil=False)

# Show a specific scope with full attributes
scope = parser.master_phil.get('picard')
if scope:
    scope.show(attributes_level=2)
"
```

`attributes_level=2` reveals `.multiple`, `.expert_level`, `.short_caption`, `.type`, value constraints, etc.

### Finding PHIL source files in the CCP4 installation

PHIL definitions typically live alongside the tool's Python package:

```bash
# Find where a package is installed
ccp4-python -c "import phasertng; print(phasertng.__file__)"

# Common locations for PHIL definitions:
#   <package>/phil/master_phil_file.py
#   <package>/phil/master_auto_file.py
#   <package>/programs/<tool>.py  (often has master_phil as a class attribute)
```

## Reference Implementation

The `phasertng_picard` wrapper serves as the canonical example:

| File | Purpose |
|------|---------|
| `wrappers/phasertng_picard/script/phasertng_picard.py` | Plugin script |
| `wrappers/phasertng_picard/script/phasertng_picard.def.xml` | inputData/outputData definition |
| `core/PhilPluginScript.py` | Base class |
| `utils/phil_to_cdata.py` | PHIL → CData converter |
| `utils/phil_shims.py` | Rich file type → PHIL value bridges |
| `tests/test_phil_to_cdata.py` | Phil2CData unit tests |
| `tests/test_phil_shims.py` | Shim unit tests |
| `tests/test_phil_plugin_script.py` | PhilPluginScript unit tests |
| `tests/i2run/test_phasertng_picard.py` | Integration tests |

## Design Decisions

### Why not generate `.def.xml` from PHIL?

Considered and rejected. PHIL types don't map 1:1 to CCP4i2's XML schema (e.g., PHIL has no concept of `CMtzDataFile` with column selection, `CAsuDataFile` with sequence parsing, or `CDictDataFile` with dictionary preview). Generating XML would either lose CCP4i2 features or require complex translation logic that's harder to maintain than the shim approach.

### Why runtime conversion instead of a build step?

The PHIL tree is available at runtime via a simple import. A build step would add complexity, create stale artifacts, and break when the upstream tool updates its PHIL definitions. Runtime conversion means the parameter UI always matches the installed tool version.

### Why keep `.def.xml` for files?

CCP4i2's rich file types provide features that plain PHIL paths don't: file browser integration, column selection UI for MTZ files, ASU content editing, database registration, file validation, and autopopulation from previous jobs. These features justify the small amount of XML.

### Why `master_phil.fetch()` instead of string concatenation?

PHIL's `fetch()` method handles type validation, scope merging, proper quoting, and default propagation correctly. String concatenation (the old approach) was fragile and could produce invalid PHIL that would fail at runtime with cryptic errors.
