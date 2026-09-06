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

    #: Where the tool's PHIL comes from (see "Declaring the PHIL source").
    PHIL_SCOPE = "my_tool:master_phil"

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

### Declaring the PHIL source

**When resolved:** During `__init__` (via `_merge_phil_parameters()`) and again
during `build_working_phil()`.

Set exactly one of these class attributes. The base class does the import and
reports a missing tool consistently, so a wrapper does not write that code.

| Attribute | Form | Use when |
|-----------|------|----------|
| `PHIL_SCOPE` | `"module.path:attribute"` | The module already exposes a scope object |
| `PHIL_PARAMS_FILE` | `"package:relative/file.params"` | The tool ships its parameters as a file inside the package |
| `PHIL_PROGRAM` | `"module.path:ClassName"` | The tool is a CCTBX program template |

```python
class my_task(PhilPluginScript):
    PHIL_SCOPE = "xia2.cli.multiplex:phil_scope"
```

**`PHIL_PROGRAM` also registers custom PHIL converters.** Tools built on the
CCTBX program template often define their own types — phasertng has
`filesystem`, `mtzcol`, `scatterer` — and those converters are registered as a
side effect of constructing `CCTBXParser`. Reaching for the scope directly
would fail to parse them, so use `PHIL_PROGRAM` for such tools even where an
attribute looks available.

**A missing tool resolves to `None`, not an exception.** The task still opens
in the interface and can say what it needs; `_merge_phil_parameters()` treats
`None` as "no parameters to merge". A raised `ImportError` would stop the task
loading at all.

**Override `get_master_phil()` directly** for anything these three do not
cover — assembling a scope from several sources, say, or one that needs
arguments. The declarations are a shortcut for the common shapes, not a
replacement for the method.

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

Output like `picard.xyzin  type=path  multiple=True` tells you to use `CList` + `PdbFileListShim` rather than a single `CPdbDataFile` + `PdbFileShim` for an *input file* you want typed on the CCP4i2 side.

### Style tokens

`.style` is free text the Phenix GUI interprets. The tokens a CCP4i2 GUI can
act on become qualifiers (`parse_phil_style` in `phil_to_cdata.py`); the rest
stay in the raw `style` qualifier and are ignored:

| token | qualifier |
|---|---|
| `spinner min=N max=M` | `min`/`max`, only where the `.type` declared none |
| `hidden` | `hidden: True` — the container element does not draw it |
| `height:N` | `guiMode: multiLine` (str only) |
| `input_file`, `file_type:X` | `philInputFile`, `philFileType` — a tag, not a file object: files belong in `inputData` via a shim, and the tag is the list of shims to write |
| `directory` | `isDirectory` |
| `phaser:mode:A,B*` or `tng:input:+a+b` | `philModes: [...]` |
| `phaser:ignore` | `philIgnored: True` (recorded, not hidden — these are the GUI's own control-flow choices) |

### One task per mode

A tool that runs in one of several modes and tags its parameters by mode
(`phaser:mode:`, `tng:input:`) becomes one CCP4i2 task per mode: set
`PHIL_MODE = "EP_AUTO"` and `PHIL_MODE_PATH = "phaser.mode"` on the
wrapper. Only the parameters whose tags match are offered — an untagged
parameter takes the tag of its nearest tagged scope, `*` and `MR*` are
wildcards (`match_modes` in `phil_to_cdata.py`, Phaser's own rule) — the
mode parameter leaves the tree, and the working phil opens with
`phaser.mode = EP_AUTO`.

### A worked example: `phaser_mr_auto_phil`

`wrappers/phaser_mr_auto_phil` is Phaser's MR_AUTO as a mode task: `PHIL_MODE`
fixes the mode, the def.xml declares only typed inputs (reflections, search
models, a composition source), and the shims in
`wrappers/phaser_phil/script/phaser_shims.py` write `phaser.hklin`/`labin`,
`phaser.ensemble`/`search` blocks and `phaser.composition` from them. Phaser
runs in-process (`phaser_run.run_mode`): Phaser's own driver builds the
`phaser.Input` object from the working phil, and a `PhaserRecorder` is the
callback that writes `program.xml` as the run proceeds -- module timeline,
progress, warnings, graphs, the verdict. After the run the solutions come
from the `mr_solution` object (typed fields; only the annotation is
tokenised, against the grammar Phaser documents), and the search-strategy
narrative from the fixed control sentences in Phaser's summary blocks, with
a count of any block that matched none of them. Nothing in the report is
inferred from prose. `phaser_ep_auto_phil` is EP_AUTO the same way: an
`EpCrystalShim` writes the one crystal block (anomalous pairs, labels read
back from the file with iotbx, substructure as `crystal.pdb_file`), the
hands, sites and figures of merit come from `ResultEP`, and the
substructure-completion cycles from the SAD summary block. The sections the
two reports share live in `phaser_report_base.py`.

The other MR modes are thin subclasses of `phaser_mr_auto_phil`, each setting
`PHIL_MODE` and composing the base class's harvest steps (`harvestCoordinates`,
`writeSolutions`, `recordRun`) as the mode's Result allows:

| task | mode | takes | writes |
|---|---|---|---|
| `phaser_mr_frf_phil` | MR_FRF | search models | a rotation list (`RFILEOUT`) |
| `phaser_mr_ftf_phil` | MR_FTF | a rotation list (`RFILEIN`) | solutions (`SOLOUT`) |
| `phaser_mr_pak_phil` | MR_PAK | solutions (`SOLIN`) | the ones that pack |
| `phaser_mr_rnp_phil` | MR_RNP | solutions, or ensembles placed at origin | refined solutions, models, maps |

A mode that works on placed solutions sets `SEARCHES_ENSEMBLES = False` so
the ensembles need not ask for copies, and `SOLUTION_INPUT` names the typed
input the `SolutionHook` hands to `setSOLU`. Phaser's rule on the kind of
file -- the translation function wants a rotation list, every other mode
solutions with translations -- is checked before the run (`solutions_kind_check`,
code 118), as is that the solutions name this job's ensembles (code 117).

### A pipeline that hosts the tool's PHIL

`pipelines/phaser_pipeline_phil` runs `phaser_mr_auto_phil` as a sub-job and
then sheetbend and refmac. Its parameters *are* the task's: the same
`PHIL_PARAMS_FILE` and `PHIL_MODE`, and `get_phil_exclude_scopes()` returns
the task's exclusions plus `phaser_mr_auto_phil.phil_shim_targets()`, so the
two trees have the same shape and `hand_phil_to(sub_job)` copies every set
value across. The pipeline keeps no keyword snapshot of its own; its typed
inputs are the task's plus what it adds (a Free-R set, a reference structure,
two switches), copied by name. The sub-job's `xml_responders` let the pipeline
embed the live record in its own `program.xml`. `phaser_simple_phil` is the
one-model case, building the ensemble list from `XYZIN` before validation and
before the run. `phaser_rnp_pipeline_phil` hosts MR_RNP: a parent model cut
into rigid bodies by atom selections, each an ensemble placed at the origin of
its own coordinates (`FIXENSEMBLES`, Phaser's `solution_at_origin`), so no
solution file changes hands.

One trap: the base constructor asks a task for its shims -- to keep their
targets out of the parameter tree -- before the subclass `__init__` runs, so
a shim that needs the plugin must be created lazily (a property), not in
`__init__`.

### Replacing a classic task

A PHIL task that replaces a def.xml one is named as its `successor` in
`tasks.py`. The classic task stays registered -- its jobs still open and
their reports render -- but the chooser no longer offers it (the task lookup
reports `supersededBy`), and cloning one of its jobs makes a job of the
successor, which adopts the old job's front page through
`PhilPluginScript.adopt_legacy_container(old_container)`: typed inputs of
the same name, inputs renamed in `LEGACY_INPUT_RENAMES` (new name -> old),
and values that were parameters there and are PHIL here, listed in
`LEGACY_PHIL_VALUES` (old name -> PHIL path; the Phaser tasks map the
resolution limits). Everything else takes the PHIL defaults: the old job's
keyword snapshot is not carried over, by design.

### Shims write blocks too, and own their targets

A shim's `convert()` may return `(path, [entries])` for one instance of a
repeated scope, inner paths relative to the scope, alongside plain
`(path, value)` pairs. And every path a shim writes (its `phil_*`
attributes, see `PhilShim.phil_targets()`) is excluded from the generic tree
automatically, so a file chosen as a typed input is not also offered as a
bare path string under the parameters.

### Repeated scopes and definitions in controlParameters

Nothing needs doing for `.multiple` parameters that stay in `controlParameters`:
`Phil2CData` turns a `.multiple = True` **scope** into a `CList` whose items are
containers shaped like the scope (Phaser's `composition.chain`, `ensemble`,
`crystal.dataset` — nested repeats included), and a `.multiple = True`
**definition** into a `CList` of the leaf type. Both start empty, which is
what libtbx's `fetch()` gives the tool unless the working phil supplies
instances; each `makeItem()` carries the scope's defaults. `build_working_phil()`
writes one `path { ... }` block per item, and one `path = value` line per item
of a repeated definition.

Two libtbx facts shape this. A block with no assignments in it is dropped by
`fetch()`, and so is an instance whose values all equal the master template —
PHIL cannot express "a repeat that is all defaults", so an item the user adds
and leaves untouched is not an instance as far as the tool is concerned. Only
user-set values are written inside a block, exactly as at the top level.

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

Five more tasks were ported to this pattern from generated `.def.xml` files,
and between them cover the awkward cases:

| Task | What it shows |
|------|---------------|
| `xia2_dials` | The plain shape: exclusions, arguments after the phil file, a `validity()` that reads PHIL-built paths |
| `xia2_xds` | Subclassing a ported task; narrowing a choice whose upstream default is not valid here |
| `xia2_multiplex` | Arguments built from a `CList` of input files |
| `xia2_ssx_reduce` | CCP4i2's own parameters surviving the merge, and one delivered as a nested phil file |
| `phaser_phil` | The minimal case — two declarations and `get_command_target()`, nothing else |

## Porting a Task That Already Has a Generated `.def.xml`

Some tasks predate `PhilPluginScript` and carry a `.def.xml` **generated
offline** from an upstream PHIL scope and checked in, usually by a
`create_def_xml.py` beside it. Those are exactly the tasks worth porting: the
checked-in XML is a snapshot that stops tracking upstream the moment it is
written.

**Find them:** a `create_def_xml.py` in the task's `script/` directory, or
`content id` values joined with `__` (`dials__index__method`), which is what
`PhilTaskCreator` emits and nobody writes by hand.

### Check first: is upstream well enough annotated?

Generated def.xml files sometimes exist because upstream lacked
`short_caption` and `expert_level`, and someone vendored a patched copy of the
scope to add them. Before porting, compare:

```python
for label, scope in (("upstream", upstream_phil), ("vendored", vendored_phil)):
    defs = list(scope.all_definitions())
    print(label, len(defs),
          sum(1 for o in defs if o.object.expert_level is not None),
          sum(1 for o in defs if o.object.short_caption is not None))
```

When the xia2 family was ported, upstream had *overtaken* its 2016 fork —
338 definitions against 167, 124 expert levels against 84, 241 short captions
against 104 — so the port lost one expert level and no captions at all. It
also *gained* labels, because `Phil2CData` falls back to a scope's name where
upstream gives no caption: 67 of 67 folders labelled, against 28 of 66 in the
generated XML.

If upstream is genuinely poorer, porting trades freshness for annotation. Say
so and decide deliberately.

### The port

1. Subclass `PhilPluginScript`; declare the PHIL source (`PHIL_SCOPE` and
   friends), add `get_command_target()` and `PHIL_EXCLUDE_SCOPES`. The old `_elts_to_remove` in `create_def_xml.py`
   translates directly — swap `__` for `.`.
2. Delete `extract_parameters()` and the hand-rolled phil writing. Override
   `makeCommandAndScript()` only if arguments follow the phil file.
3. Cut the `.def.xml` down to `inputData` and `outputData`, **keeping an empty
   `controlParameters` container** — `_merge_phil_parameters()` merges into it
   and cannot create it.
4. Delete `create_def_xml.py` and any vendored phil module.
5. Register the interface (see below) and run the tool end to end.

Parameter paths do not change: `Phil2CData` builds the same nested containers
with the same `__`-joined names as `PhilTaskCreator` did, so saved
`input_params.xml` still loads and any `validity()` accessors still resolve.

### Traps

**Command-line arguments must follow the phil file.** `build_working_phil()`
writes the *whole fetched scope, defaults included* — unlike the hand-rolled
approach, which wrote only what the user had set. An argument placed before the
file is overridden by the default inside it. `xia2 pipeline=3dii working.phil`
silently ran DIALS.

**A tool's default may not be legal for your task.** `xia2_xds` offers only the
XDS pipelines, but xia2's scope marks `dials` as default. Narrow the
enumerators *and* assert the choice on the command line; narrowing alone leaves
an out-of-range value that fails validation.

**Parameters that are CCP4i2's, not the tool's**, stay declared in the
`.def.xml`: `_merge_phil_parameters()` skips names already present. They must
then be kept *out* of the phil by overriding `extract_phil_parameters()`, or
they are written under names the tool does not know. `xia2_ssx_reduce` has two
— one shown in the interface and never sent, one delivered indirectly as a
phil file of its own because xia2 passes `symmetry.phil` through to
`dials.cosym` rather than exposing its parameters.

**Keep the `.def.xml` encoding declaration as UTF-8.** At least one guiLabel in
the tree contains an en-dash; rewriting the header as ASCII makes the file
unparseable.

**Guard empty list slots.** An unset entry in a `CList` of files becomes
`reflections=` on the command line, which the tool reads as a bad argument.

## The Interface

A PHIL task's parameters are far too numerous to enumerate by hand — xia2/DIALS
has 286. Render the container and let expert level do the filtering.

```tsx
import { ExpertLevelContext } from "../task-elements/expert-level-context";
import {
  EXCLUDE_EXPERT_LEVEL,
  PhilExpertLevelSelector,
  usePhilExpertLevel,
} from "../task-elements/phil-expert-level";

const { expertLevel, changeExpertLevel } = usePhilExpertLevel(props.job);

<PhilExpertLevelSelector expertLevel={expertLevel} onChange={changeExpertLevel} />
<ExpertLevelContext.Provider value={expertLevel}>
  <CCP4i2ContainerElement
    {...props}
    itemName="controlParameters"
    qualifiers={{ guiLabel: "Parameters" }}
    containerHint="FolderLevel"
    excludeItems={EXCLUDE_EXPERT_LEVEL}
  />
</ExpertLevelContext.Provider>
```

Lay out by hand only what needs prose around it — the file inputs, and any
pair of parameters that are alternatives.

**`PHIL_EXPERT_LEVEL` is not a display setting.** `extract_phil_parameters()`
also uses it to decide which parameters reach `working.phil`, so the interface
must bind to the container parameter rather than keep its own copy in React
state; two copies can disagree and silently drop a value the user set. That is
what `usePhilExpertLevel` is for. Exclude it from the rendered tree — it is a
control, not a parameter.

**Where expert level lives differs by route.** PHIL-generated containers carry
it as a plain qualifier; def.xml tasks nest it inside `guiDefinition`.
`CCP4i2ContainerElement` reads either, so an interface works before and after a
port.

## Design Decisions

### Why not generate `.def.xml` from PHIL?

Considered and rejected — and the tasks that did it bear the argument out. The xia2 family's generated files had fallen 34 to 171 definitions behind upstream, and one carried a vendored fork of the scope that upstream had long since overtaken. PHIL types don't map 1:1 to CCP4i2's XML schema (e.g., PHIL has no concept of `CMtzDataFile` with column selection, `CAsuDataFile` with sequence parsing, or `CDictDataFile` with dictionary preview). Generating XML would either lose CCP4i2 features or require complex translation logic that's harder to maintain than the shim approach.

### Why runtime conversion instead of a build step?

The PHIL tree is available at runtime via a simple import. A build step would add complexity, create stale artifacts, and break when the upstream tool updates its PHIL definitions. Runtime conversion means the parameter UI always matches the installed tool version.

### Why keep `.def.xml` for files?

CCP4i2's rich file types provide features that plain PHIL paths don't: file browser integration, column selection UI for MTZ files, ASU content editing, database registration, file validation, and autopopulation from previous jobs. These features justify the small amount of XML.

### Why `master_phil.fetch()` instead of string concatenation?

PHIL's `fetch()` method handles type validation, scope merging, proper quoting, and default propagation correctly. String concatenation (the old approach) was fragile and could produce invalid PHIL that would fail at runtime with cryptic errors.
