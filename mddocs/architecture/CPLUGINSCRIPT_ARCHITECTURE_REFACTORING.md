# CPluginScript Architecture Refactoring

**Date**: 2025-10-31
**Issue**: Current utilities create standalone CContainers, bypassing CPluginScript + dbHandler architecture
**Impact**: File parameters don't persist properly, no Django ↔ CData synchronization
**Solution**: Refactor to use CPluginScript with dbHandler as the starting point for ALL operations

---

## Current Problem

### What We're Doing Now (WRONG)

```python
# server/ccp4x/lib/utils/parameters/set_param.py
def set_parameter(job: Job, path: str, value: Any) -> Result[Dict]:
    # ❌ Creates standalone container with no database connection
    container = get_job_container(job.uuid)

    # ❌ Sets parameter on disconnected container
    set_object_value(container, path, value)

    # ❌ Manually saves to XML
    save_params_xml(container, job.directory / "params.xml")

    # ❌ File parameters don't persist properly!
    return Result.ok({"path": path, "value": value})
```

### What We SHOULD Be Doing (RIGHT)

```python
# Using CPluginScript + dbHandler
def set_parameter(job: Job, path: str, value: Any) -> Result[Dict]:
    # ✅ Get plugin instance with dbHandler attached
    plugin = get_job_plugin_with_db(job)

    # ✅ Set parameter through plugin's container
    # This triggers proper file handling through CDataFile.setFullPath()
    set_object_value(plugin.container, path, value)

    # ✅ Save through plugin (handles DB sync)
    plugin.container.saveDataToXml(str(job.directory / "params.xml"))

    # ✅ Update database via dbHandler
    if plugin._dbHandler:
        plugin._dbHandler.updateJobStatus(
            jobId=str(job.uuid),
            container=plugin.container
        )

    # ✅ File parameters persist correctly!
    return Result.ok({"path": path, "value": value})
```

---

## Why CPluginScript + dbHandler Matters

### 1. Database Integration
- **dbHandler** provides `updateJobStatus(container=...)` for syncing CData → Django
- Plugin has `_dbJobId`, `_dbJobNumber`, `_dbProjectId` for context
- File operations can query/update database records

### 2. File Path Resolution
- **CDataFile.setFullPath()** checks for `plugin._dbHandler` via `_find_plugin_parent()`
- When DB-aware, calls `_update_from_database(path, plugin)` to get dbFileId
- Creates/updates Django File records automatically

### 3. Proper Hierarchy
```
CPluginScript (has dbHandler, knows job context)
  └── container (CContainer)
       ├── inputData (CContainer)
       │    └── XYZIN (CDataFile) ← knows parent is plugin, can access dbHandler
       ├── outputData (CContainer)
       └── controlParameters (CContainer)
```

### 4. Signal System
- Plugin emits `finished` signal when done
- Django handlers can connect to update UI/status
- Async workflows work properly

---

## Architecture Pattern

### The Golden Rule

**Every operation on job parameters MUST:**
1. Start by getting `CPluginScript` instance via `get_job_plugin()`
2. Attach `CCP4i2DjangoDbHandler` instance
3. Set job context: `plugin.set_db_job_id(job.uuid)`
4. Operate on `plugin.container`, NOT standalone containers
5. Save via `plugin.container.saveDataToXml()`
6. Update DB via `plugin._dbHandler.updateJobStatus()` if needed

---

## Implementation Plan

### Step 1: Create Unified Plugin Context Utility

**File**: `server/ccp4x/lib/utils/plugins/plugin_context.py`

```python
"""
Unified utility for getting CPluginScript instances with database context.
This is the ONE place we create plugins - all other utilities use this.
"""
import logging
from typing import Optional
from pathlib import Path

from ccp4x.db.models import Job, Project
from ccp4x.db.ccp4i2_django_db_handler import CCP4i2DjangoDbHandler
from ccp4x.lib.job_utils.get_job_plugin import get_job_plugin
from core.CCP4PluginScript import CPluginScript
from ccp4x.lib.response import Result

logger = logging.getLogger(__name__)


def get_plugin_with_context(
    job: Job,
    params_file: Optional[Path] = None,
    create_db_handler: bool = True
) -> Result[CPluginScript]:
    """
    Get CPluginScript instance with full database context.

    This is the CANONICAL way to get a plugin for any operation.

    Args:
        job: Django Job model instance
        params_file: Optional specific params file to load
        create_db_handler: Whether to attach dbHandler (default True)

    Returns:
        Result[CPluginScript] with plugin instance or error

    Example:
        result = get_plugin_with_context(job)
        if result.success:
            plugin = result.data
            # Now use plugin.container for all operations
            set_object_value(plugin.container, "inputData.XYZIN", "/path/to/file.pdb")
            plugin.container.saveDataToXml(str(job.directory / "params.xml"))
    """
    try:
        # Create dbHandler
        dbHandler = None
        if create_db_handler:
            dbHandler = CCP4i2DjangoDbHandler(project_id=str(job.project.uuid))

        # Get plugin instance (loads params.xml or input_params.xml)
        plugin = get_job_plugin(
            the_job=job,
            parent=None,
            dbHandler=dbHandler
        )

        if plugin is None:
            return Result.fail(
                f"Failed to load plugin for task '{job.task_name}'",
                details={"job_id": str(job.uuid), "task_name": job.task_name}
            )

        # Set database context on plugin
        plugin.set_db_job_id(str(job.uuid))
        plugin.set_db_job_number(job.number)
        plugin._dbProjectId = str(job.project.uuid)
        plugin._dbProjectName = job.project.name

        logger.debug(
            "Loaded plugin %s with DB context: job=%s, project=%s",
            job.task_name, job.uuid, job.project.uuid
        )

        return Result.ok(plugin)

    except Exception as e:
        logger.exception("Failed to get plugin with context for job %s", job.uuid)
        return Result.fail(
            f"Error loading plugin: {str(e)}",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "error_type": type(e).__name__
            }
        )
```

### Step 2: Refactor set_parameter

**File**: `server/ccp4x/lib/utils/parameters/set_param.py`

```python
"""
Set job parameter using CPluginScript architecture.
"""
import logging
from typing import Any
from pathlib import Path

from ccp4x.db.models import Job
from ccp4x.lib.response import Result
from ccp4x.lib.utils.plugins.plugin_context import get_plugin_with_context
from ccp4x.lib.job_utils.set_parameter import set_object_value

logger = logging.getLogger(__name__)


def set_parameter(job: Job, path: str, value: Any) -> Result[dict]:
    """
    Set a parameter on a job using proper CPluginScript + dbHandler architecture.

    Args:
        job: Django Job instance
        path: Object path (e.g., "inputData.XYZIN" or "inputData.XYZIN.baseName")
        value: Value to set (can be file path, number, string, etc.)

    Returns:
        Result[dict] with parameter info or error

    Example:
        # Set file parameter
        result = set_parameter(job, "inputData.XYZIN", "/path/to/file.pdb")

        # Set nested attribute
        result = set_parameter(job, "inputData.XYZIN.annotation", "Input structure")

        # Set control parameter
        result = set_parameter(job, "container.NCYCLES", 10)
    """
    # Get plugin with database context
    plugin_result = get_plugin_with_context(job)
    if not plugin_result.success:
        return Result.fail(f"Failed to load plugin: {plugin_result.error}")

    plugin = plugin_result.data

    try:
        # Set parameter through plugin's container
        # This ensures proper file handling, validation, and hierarchy
        logger.debug("Setting parameter %s = %s on job %s", path, value, job.uuid)
        set_object_value(plugin.container, path, value)

        # Save parameters to XML
        params_file = job.directory / "params.xml"
        logger.debug("Saving parameters to %s", params_file)
        plugin.container.saveDataToXml(str(params_file), check=False)

        # Update database via dbHandler (if available)
        if plugin._dbHandler:
            logger.debug("Updating database via dbHandler for job %s", job.uuid)
            plugin._dbHandler.updateJobStatus(
                jobId=str(job.uuid),
                container=plugin.container
            )

        # Get the actual object for return info
        obj = plugin.container
        parts = path.split('.')
        for part in parts:
            if hasattr(obj, part):
                obj = getattr(obj, part)
            else:
                obj = None
                break

        # Return detailed info about what was set
        result_data = {
            "path": path,
            "value": value,
            "object_type": type(obj).__name__ if obj else "Unknown",
        }

        # Add file-specific info if it's a CDataFile
        if obj and hasattr(obj, 'getFullPath'):
            result_data["file_path"] = obj.getFullPath()
            if hasattr(obj, 'dbFileId') and hasattr(obj.dbFileId, 'value'):
                result_data["db_file_id"] = str(obj.dbFileId.value)

        logger.info("Successfully set parameter %s on job %s", path, job.uuid)
        return Result.ok(result_data)

    except AttributeError as e:
        return Result.fail(
            f"Parameter path '{path}' not found",
            details={"path": path, "error": str(e)}
        )
    except Exception as e:
        logger.exception("Failed to set parameter %s on job %s", path, job.uuid)
        return Result.fail(
            f"Error setting parameter: {str(e)}",
            details={"path": path, "value": value, "error_type": type(e).__name__}
        )
```

### Step 3: Refactor validate_job

**File**: `server/ccp4x/lib/utils/jobs/validate.py`

```python
"""
Validate job parameters using CPluginScript architecture.
"""
import logging
from xml.etree import ElementTree as ET

from ccp4x.db.models import Job
from ccp4x.lib.response import Result
from ccp4x.lib.utils.plugins.plugin_context import get_plugin_with_context

logger = logging.getLogger(__name__)


def validate_job(job: Job) -> Result[ET.Element]:
    """
    Validate job parameters using plugin's validation system.

    Args:
        job: Django Job instance

    Returns:
        Result[ET.Element] with error report XML or error
    """
    # Get plugin with database context
    plugin_result = get_plugin_with_context(job)
    if not plugin_result.success:
        return Result.fail(f"Failed to load plugin: {plugin_result.error}")

    plugin = plugin_result.data

    try:
        # Use plugin's built-in validation (validates container)
        logger.debug("Validating job %s (task: %s)", job.uuid, job.task_name)

        # CPluginScript has container.validateInput() method
        error_report = plugin.container.validateInput()

        # Convert error report to XML
        error_tree = error_report.toEtree()

        errors = error_tree.findall('.//error')
        warnings = error_tree.findall('.//warning')

        logger.info(
            "Validation complete for job %s: %d errors, %d warnings",
            job.uuid, len(errors), len(warnings)
        )

        return Result.ok(error_tree)

    except Exception as e:
        logger.exception("Failed to validate job %s", job.uuid)
        return Result.fail(
            f"Validation failed: {str(e)}",
            details={"job_id": str(job.uuid), "error_type": type(e).__name__}
        )
```

### Step 4: Refactor get_job_reports

**File**: `server/ccp4x/lib/utils/jobs/reports.py`

```python
"""
Get job reports (params, report, diagnostic) using CPluginScript architecture.
"""
import logging
from xml.etree import ElementTree as ET
from pathlib import Path

from ccp4x.db.models import Job
from ccp4x.lib.response import Result
from ccp4x.lib.utils.plugins.plugin_context import get_plugin_with_context

logger = logging.getLogger(__name__)


def get_job_params_xml(job: Job, regenerate: bool = False) -> Result[ET.Element]:
    """
    Get params XML for job using plugin's container.

    Args:
        job: Django Job instance
        regenerate: If True, regenerate from container; if False, read existing file

    Returns:
        Result[ET.Element] with params XML or error
    """
    if not regenerate:
        # Try to read existing params.xml
        params_file = job.directory / "params.xml"
        if not params_file.exists():
            params_file = job.directory / "input_params.xml"

        if params_file.exists():
            try:
                tree = ET.parse(str(params_file))
                return Result.ok(tree.getroot())
            except ET.ParseError as e:
                logger.warning("Failed to parse existing params.xml: %s", e)
                # Fall through to regenerate

    # Regenerate from plugin container
    plugin_result = get_plugin_with_context(job)
    if not plugin_result.success:
        return Result.fail(f"Failed to load plugin: {plugin_result.error}")

    plugin = plugin_result.data

    try:
        # Save container to XML and parse
        params_file = job.directory / "params.xml"
        plugin.container.saveDataToXml(str(params_file), check=False)

        tree = ET.parse(str(params_file))
        return Result.ok(tree.getroot())

    except Exception as e:
        logger.exception("Failed to get params XML for job %s", job.uuid)
        return Result.fail(
            f"Error generating params XML: {str(e)}",
            details={"job_id": str(job.uuid), "error_type": type(e).__name__}
        )


def get_job_report_xml(job: Job) -> Result[ET.Element]:
    """
    Get report XML for job (if job has finished).

    Args:
        job: Django Job instance

    Returns:
        Result[ET.Element] with report XML or error
    """
    report_file = job.directory / "report.xml"

    if not report_file.exists():
        return Result.fail(
            "Report file not found (job may not have finished)",
            details={"job_id": str(job.uuid), "report_file": str(report_file)}
        )

    try:
        tree = ET.parse(str(report_file))
        return Result.ok(tree.getroot())
    except ET.ParseError as e:
        return Result.fail(
            f"Failed to parse report.xml: {str(e)}",
            details={"job_id": str(job.uuid), "report_file": str(report_file)}
        )


def get_job_diagnostic_xml(job: Job) -> Result[ET.Element]:
    """
    Get diagnostic XML by validating job parameters.

    Args:
        job: Django Job instance

    Returns:
        Result[ET.Element] with diagnostic XML or error
    """
    from .validate import validate_job
    return validate_job(job)
```

---

## Migration Strategy

### Phase 1: Create Infrastructure (30 min)
1. ✅ Create `server/ccp4x/lib/utils/plugins/plugin_context.py`
2. ✅ Test `get_plugin_with_context()` with simple job

### Phase 2: Refactor Core Utilities (1 hour)
1. ✅ Refactor `set_parameter()` to use plugin context
2. ✅ Refactor `validate_job()` to use plugin context
3. ✅ Refactor `get_job_reports()` to use plugin context
4. ✅ Test each utility as it's refactored

### Phase 3: Test End-to-End (30 min)
1. ✅ Create job
2. ✅ Set file parameter via `set_parameter()`
3. ✅ Verify parameter persists in params.xml
4. ✅ Verify dbFileId is set correctly
5. ✅ Validate job works
6. ✅ Get reports work

### Phase 4: Update Other Utilities (Optional)
- `clone_job` - Already uses `get_job_plugin()` ✅
- `run_job` - Already uses plugin architecture ✅
- `export_job` - Could benefit from plugin context

---

## Benefits

### Immediate
1. ✅ **File parameters persist correctly** - CDataFile.setFullPath() with DB awareness
2. ✅ **Database sync works** - dbHandler.updateJobStatus() called properly
3. ✅ **Validation works** - Uses plugin's built-in validation
4. ✅ **Proper hierarchy** - All objects know their parent plugin

### Long-term
1. ✅ **Consistent architecture** - ONE way to work with jobs
2. ✅ **Maintainability** - Single point of truth for plugin creation
3. ✅ **Extensibility** - Easy to add new operations
4. ✅ **Type safety** - Plugin has all attributes we need
5. ✅ **Signal system** - Can connect to plugin.finished for async workflows

---

## Key Insight

**Your original observation was spot-on**: We were creating orphaned CContainers and manually juggling XML files. The CPluginScript + dbHandler architecture was DESIGNED to handle all of this, but we bypassed it.

By using `get_plugin_with_context()` as the entry point for ALL operations, we get:
- ✅ Proper file handling
- ✅ Database synchronization
- ✅ Validation
- ✅ Signal system
- ✅ Error reporting
- ✅ Async support

All of this **already exists** in CPluginScript - we just need to use it!

---

**Next Step**: Implement Phase 1 (create `plugin_context.py`) and test it?
