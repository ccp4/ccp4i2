# Task Manager and Plugin/Report Registry

This directory contains the plugin and report registry systems for cdata-codegen.

## Directory Structure

### Plugin Registry Files
- **plugin_registry.py** - Auto-generated lazy-loading plugin registry
- **plugin_lookup.json** - Plugin metadata index
- **plugin_lookup.py** - Script to regenerate plugin registry

### Report Registry Files
- **report_registry.py** - Auto-generated lazy-loading report registry
- **report_lookup.json** - Report metadata index
- **report_lookup.py** - Script to regenerate report registry

### DEF.XML and Parameter Files
- **defxml_lookup.json** - Index of .def.xml files for plugin definitions
- **defxml_lookup.py** - Script to regenerate defxml index
- **def_xml_handler.py** - Parser for .def.xml plugin definition files
- **params_xml_handler.py** - Parser for .params.xml parameter files

### UI Metadata
- **task_module_map.json** - Maps task names to UI module categories
- **task_metadata.json** - UI display metadata (titles, descriptions)

---

## Regenerating All Lookup Files

Use the CTaskManager `--rebuild` option to regenerate all lookup files at once:

```bash
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
.venv/bin/python core/CCP4TaskManager.py --rebuild
```

This will regenerate:
1. `defxml_lookup.json` - Index of .def.xml files
2. `plugin_lookup.json` and `plugin_registry.py` - Plugin metadata and lazy loader
3. `report_lookup.json` and `report_registry.py` - Report metadata and lazy loader

---

## Plugin Registry

### Regenerating Plugin Registry Only

```bash
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
.venv/bin/python core/task_manager/plugin_lookup.py
```

**Important Notes:**
1. Use the virtual environment's Python (`.venv/bin/python`) to ensure all dependencies are available
2. The script scans `wrappers/`, `wrappers2/`, and `pipelines/` directories
3. Plugins that fail to import (due to missing dependencies like qtgui, mmut, etc.) will be skipped with warnings
4. The generated files contain metadata for ~144 plugins that can be successfully imported
5. Module paths are automatically corrected to be relative to CCP4I2_ROOT

### Plugin Registry Usage

```python
from core.CCP4TaskManager import TASKMANAGER

task_mgr = TASKMANAGER()

# Get plugin class (lazy loaded on first access)
plugin_class = task_mgr.get_plugin_class('ctruncate')

# Get plugin metadata without importing
metadata = task_mgr.get_plugin_metadata('ctruncate')

# List all available plugins
plugins = task_mgr.list_plugins()  # Returns ~144 plugin names

# Locate .def.xml file for a plugin
def_xml_path = task_mgr.locate_def_xml('refmac')
```

### Available Plugins (~144)

- **Data Processing**: ctruncate, aimless, pointless, dials_*, xia2_*
- **Phasing**: phaser_*, shelx_*, crank2_*
- **Refinement**: refmac, buster, prosmart
- **Model Building**: buccaneer, parrot, coot_*
- **Analysis**: molprobity, edstats, moorhen_*
- **Utilities**: pdbset, rwcontents, superpose

---

## Report Registry

### Overview

The report registry discovers report classes (subclasses of `Report` from `report.CCP4ReportParser`) in the plugin directories. Reports are used to generate job status/result displays in the UI.

### Regenerating Report Registry Only

```bash
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
.venv/bin/python core/task_manager/report_lookup.py
```

**Discovery Process:**
1. Scans `wrappers/`, `wrappers2/`, and `pipelines/` for `*_report.py` files
2. Identifies classes inheriting from `Report`
3. Extracts attributes: `TASKNAME`, `TASKTITLE`, `TASKVERSION`, `RUNNING`, `WATCHED_FILE`, `FAILED`, `USEPROGRAMXML`, `SEPARATEDATA`, `ERROR_CODES`
4. Generates both `report_lookup.json` (for reference) and `report_registry.py` (lazy loader)

### Report Registry Usage

```python
from core.CCP4TaskManager import TASKMANAGER

task_mgr = TASKMANAGER()

# Check if a report exists for a task (fast, no import)
if task_mgr.has_report('refmac'):
    print("refmac has a report class")

# Get report metadata without importing
metadata = task_mgr.get_report_metadata('refmac')
print(metadata['RUNNING'])        # True - supports live updates
print(metadata['WATCHED_FILE'])   # File to watch for changes

# Get report class (lazy loaded)
report_class = task_mgr.getReportClass('refmac')

# Get specific attribute (tries metadata first, imports if needed)
running = task_mgr.getReportAttribute('refmac', 'RUNNING')

# List all available reports
reports = task_mgr.list_reports()  # Returns ~136 report names
```

### Report Metadata Fields

| Field | Type | Description |
|-------|------|-------------|
| `TASKNAME` | str | Task name this report is for |
| `TASKTITLE` | str | Human-readable title |
| `RUNNING` | bool | If True, report supports live updates while job runs |
| `WATCHED_FILE` | str/None | File to monitor for changes during job execution |
| `FAILED` | bool | Whether this is a failure report |
| `USEPROGRAMXML` | bool | Whether to use program.xml for data |
| `SEPARATEDATA` | bool | Whether data is stored separately |
| `ERROR_CODES` | dict | Error code definitions |

### Available Reports (~136)

Reports generally follow the pattern `{taskname}_report.py` with class `{taskname}_report(Report)`.

---

## DEF.XML Lookup

### Regenerating DEF.XML Lookup Only

```bash
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
python core/task_manager/defxml_lookup.py
```

This scans for `.def.xml` files and updates `defxml_lookup.json`.

### Usage

```python
from core.CCP4TaskManager import TASKMANAGER

task_mgr = TASKMANAGER()

# Locate .def.xml file for a task
def_xml_path = task_mgr.locate_def_xml('refmac')
# Returns: Path('/Users/.../wrappers/refmac_i2/script/refmac.def.xml')
```

---

## CCP4I2_ROOT Environment Variable

**CRITICAL**: All tests and the plugin/report system require `CCP4I2_ROOT` to be set:

```bash
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
```

The project root contains:
- `wrappers/` - Legacy ccp4i2 plugin wrappers (~115 plugins)
- `wrappers2/` - Additional legacy plugins
- `pipelines/` - Multi-step pipeline plugins
- `demo_data/` - Test data files
- `core/` - New Python implementation
- `server/` - Django backend
- `report/` - Report parser and base classes

---

## Import Paths

Plugins and reports are imported using their module paths relative to CCP4I2_ROOT:

- **Wrappers**: `wrappers.{plugin_name}.script.{module}`
- **Wrappers2**: `wrappers2.{plugin_name}.script.{module}`
- **Pipelines**: `pipelines.{plugin_name}.script.{module}`

The registries handle all import logic automatically.

---

## Task Tree and UI Integration

CTaskManager provides methods for building the task browser UI:

```python
task_mgr = TASKMANAGER()

# Get task tree structure for UI navigation
tree = task_mgr.task_tree(show_wrappers=False)
# Returns: [(module_name, module_title, [task_names...]), ...]

# Get task lookup for displaying task cards
lookup = task_mgr.task_lookup
# Returns: {taskName: {version: {TASKTITLE, DESCRIPTION, MAINTAINER, shortTitle}}}

# Get icon paths for module folders
icons = task_mgr.task_icon_lookup
# Returns: {module_name: 'qticons/icon.png'}
```

---

## Troubleshooting

### Plugin/Report Not Found

1. Ensure `CCP4I2_ROOT` is set correctly
2. Check if the plugin/report file exists in the expected location
3. Regenerate the registries with `--rebuild`
4. Check for import errors in the console output

### Import Errors During Regeneration

Some plugins require dependencies not available in the venv:
- `qtgui` - Qt GUI components (skip these, backend-only)
- `mmut` - mmCIF utilities (use stub in `stubs/mmut.py.stub`)
- `ccp4mg`, `mmdb2`, `ccp4srs` - Use stubs in `stubs/` directory

### Stale Registry Data

If plugins or reports have been added/removed/modified:
```bash
.venv/bin/python core/CCP4TaskManager.py --rebuild
```
