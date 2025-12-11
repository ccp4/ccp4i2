# Database Integration Summary

## What I've Learned

After examining the Django database models and the existing job execution workflow, I now understand:

### 1. Current Database Schema (models.py)

- **Project**: Tracks project directory, last job number, creation metadata
- **Job**: Hierarchical job tracking with parent/child relationships
  - Job numbers like "1", "1.1", "1.2" for nested execution
  - Status tracking (PENDING → RUNNING → FINISHED/FAILED)
  - Directory convention: `project_dir/CCP4_JOBS/job_1/job_2/...`
- **File**: Tracks input and output files
  - Directory types: JOB_DIR (1) or IMPORT_DIR (2)
  - Linked to FileType via mimetype
  - Stores content_flag and sub_type for crystallographic context
- **FileUse**: Links files to jobs with roles (IN/OUT)
- **FileImport**: Tracks imported external files with checksums
- **JobValueKey / JobFloatValue / JobCharValue**: Store performance indicators (KPIs)

### 2. Current Workflow (run_job.py → import_files.py → glean_job_files.py)

**Job Execution Stages**:
1. **Setup** - Create database handler, Qt event loop
2. **Import Files** - Copy external files to `CCP4_IMPORTED_FILES/`, register in database
3. **Set Output Paths** - Pre-assign file paths based on job directory
4. **Execute** - Run plugin in Qt event loop
5. **Glean Files** - Extract output files from container, register in database
6. **Glean KPIs** - Extract performance indicators from output data

**Key Observation**: All container interaction uses **legacy CContainer API**:
- `.qualifiers("mimeTypeName")` - String-based access, no type safety
- `getattr(item, "subType", None)` - Fragile attribute access with fallbacks
- Manual type conversions with try/except blocks
- Recursive search with `find_objects(container, lambda, True)`
- String parsing of `.objectPath()` to extract parameter names

### 3. Opportunities with New CData

The new CData system provides **much more robust** introspection:

✅ **Metadata System**: `get_merged_metadata('qualifiers')` provides type-safe access
✅ **Hierarchical Traversal**: `.childNames()` + `HierarchicalObject` for clean iteration
✅ **Type Safety**: `.value` property returns correctly-typed values
✅ **Value State Tracking**: `.isSet()`, `.getValueState()` for precise state checking
✅ **Object Paths**: `.object_path()` for debugging and parameter name extraction
✅ **Async Compatible**: Can be made fully async for non-blocking I/O

## What I've Created

### 1. Modern Async Database Handler (`async_db_handler.py`)

Provides a clean async API for database-backed job tracking:

```python
handler = AsyncDatabaseHandler(project_uuid=project.uuid)

# Automatic tracking with context manager
async with handler.track_job(plugin):
    # Job created, status updated automatically via signals
    await plugin.execute()
    # Files gleaned automatically on completion
```

**Key Features**:
- Async/await throughout - non-blocking database operations
- Automatic job lifecycle tracking via signal connections
- Context manager for clean resource management
- Type hints and modern Python patterns
- Replaces legacy `CCP4i2DjangoDbHandler` string-based API

### 2. Modern CData Utilities (`cdata_utils.py`)

Type-safe utilities for container introspection:

```python
# Find all files using hierarchical traversal
output_files = find_all_files(plugin.outputData)

# Extract metadata using CData system
for file_obj in output_files:
    metadata = extract_file_metadata(file_obj)
    print(f"File: {metadata['name']}")
    print(f"Type: {metadata['file_type']}")
    print(f"Content: {metadata.get('content_flag')}")

# Extract KPIs with type safety
kpis = find_objects_by_type(plugin.outputData, CPerformanceIndicator)
for kpi in kpis:
    values = extract_kpi_values(kpi)
    print(f"R-factor: {values.get('Rfactor')}")
```

**Functions Provided**:
- `find_all_files()` - Hierarchical file discovery
- `find_objects_by_type()` - Type-safe object search
- `extract_file_metadata()` - Type-safe metadata extraction
- `extract_parameter_name()` - Clean parameter name extraction
- `extract_kpi_values()` - Type-safe KPI extraction
- `validate_file_metadata_completeness()` - Metadata validation
- `debug_print_container_structure()` - Debugging tool

### 3. Comprehensive Documentation

- **`USAGE_EXAMPLES.md`** - Side-by-side comparison of legacy vs. modern API
- **`WORKFLOW_ANALYSIS.md`** - Detailed analysis of current workflow and improvement opportunities
- **`INTEGRATION_SUMMARY.md`** - This document

## Comparison: Legacy vs. Modern

### File Gleaning Example

**Legacy Approach** (current `glean_job_files.py`):
```python
# String-based, fragile, lots of error handling
file_type = item.qualifiers("mimeTypeName")  # ❌ String access
sub_type = getattr(item, "subType", None)     # ❌ getattr fallback
try:
    sub_type = int(sub_type)  # ❌ Manual conversion
except AttributeError:
    sub_type = None

content = getattr(item, "contentFlag", None)
try:
    content = int(content)
except AttributeError:
    content = None

# Regex parsing for parameter name
job_param_name = extract_from_first_bracketed(item.objectPath())
```

**Modern Approach** (using new utilities):
```python
# Type-safe, clean, minimal error handling needed
metadata = extract_file_metadata(file_obj)

# Everything is type-safe and validated
file_type = metadata['file_type']
sub_type = metadata.get('sub_type')  # None if not present
content_flag = metadata.get('content_flag')  # None if not present
param_name = metadata['name']  # Direct from object
```

### Container Traversal Example

**Legacy Approach**:
```python
# Lambda-based recursive search
outputs = find_objects(
    container.outputData,
    lambda a: isinstance(a, CDataFile),
    recursive=True
)
```

**Modern Approach**:
```python
# Clean hierarchical traversal using CData system
outputs = find_all_files(container.outputData)

# Or more specific searches
mtz_files = find_objects_matching(
    container.outputData,
    lambda obj: isinstance(obj, CDataFile) and
                obj.get_merged_metadata('qualifiers').get('mimeTypeName') == 'application/CCP4-mtz'
)
```

### KPI Extraction Example

**Legacy Approach**:
```python
# Manual iteration with lots of type checking
kpis = find_objects(container, lambda a: isinstance(a, CPerformanceIndicator), True)
for kpi in kpis:
    for kpi_param_name in kpi.dataOrder():  # ❌ Legacy method
        value = getattr(kpi, kpi_param_name)
        if value is None:
            continue
        try:
            if isinstance(value, CFloat):
                value = float(value)  # ❌ Manual conversion
            elif isinstance(value, CString):
                value = str(value)
            # ... more type checking
        except TypeError:
            continue
```

**Modern Approach**:
```python
# Clean, type-safe extraction
kpis = find_objects_by_type(container.outputData, CPerformanceIndicator)
for kpi in kpis:
    values = extract_kpi_values(kpi)  # Returns Dict[str, Union[float, int, str, bool]]

    # Type-safe access
    r_factor = values.get('Rfactor')  # Already converted to float
    if r_factor is not None:
        print(f"R-factor: {r_factor:.4f}")
```

## Key Improvements

| Feature | Legacy | Modern |
|---------|--------|--------|
| **Metadata Access** | `.qualifiers("key")` | `.get_merged_metadata('qualifiers')['key']` |
| **Type Safety** | `getattr()` + try/except | Direct access with type checks |
| **Traversal** | `find_objects(lambda)` | `.childNames()` + hierarchy |
| **Parameter Names** | Regex parsing | `.name` attribute |
| **Value Conversion** | Manual `int()`/`float()` | `.value` property |
| **Error Handling** | try/except everywhere | Minimal - types known upfront |
| **Execution** | Qt event loop (sync) | async/await |
| **Status Tracking** | Manual updates | Signal-based automatic |
| **API Style** | camelCase, strings | snake_case, type hints |

## Next Steps / Recommendations

### Immediate Actions

1. **Extend `AsyncDatabaseHandler`** with additional methods:
   ```python
   async def register_imported_file(self, ...):
       """For tracking input files copied to CCP4_IMPORTED_FILES"""

   async def register_job_float_value(self, job_uuid, key, value):
       """For KPI tracking"""

   async def register_job_char_value(self, job_uuid, key, value):
       """For KPI tracking"""
   ```

2. **Create async file operations module**:
   ```python
   async def import_input_files_async(job, plugin, db_handler):
       """Replace import_files.py with async version"""

   async def glean_output_files_async(job_uuid, output_container, db_handler):
       """Replace glean_job_files.py with async version using cdata_utils"""
   ```

3. **Update CPluginScript integration**:
   - Add `_dbHandler` attribute to store `AsyncDatabaseHandler` instance
   - Connect `statusChanged` signal to automatic database updates
   - Add helper methods for file gleaning

4. **Create comprehensive tests**:
   - Test modern utilities with real CData objects
   - Test async database handler with Django test database
   - Integration tests for full job lifecycle
   - Comparison tests (legacy vs. modern output should match)

### Gradual Migration Strategy

**Phase 1: Utilities (Completed ✅)**
- Created `cdata_utils.py` with modern introspection functions
- Created `async_db_handler.py` with async API
- Documented differences and improvements

**Phase 2: Async File Operations (Next)**
- Create `async_import_files.py` using `cdata_utils`
- Create `async_glean_files.py` using `cdata_utils`
- Add KPI extraction using `extract_kpi_values()`

**Phase 3: Async Job Runner**
- Create `async_run_job.py` replacing Qt event loop
- Use `AsyncDatabaseHandler.track_job()` context manager
- Support nested job execution with proper hierarchy

**Phase 4: Testing & Validation**
- Ensure output matches legacy implementation
- Performance testing (async should be faster)
- Integration tests with real plugins

**Phase 5: Deprecation**
- Mark legacy handlers as deprecated
- Update documentation
- Provide migration guide

### Open Questions for Discussion

1. **Container API**: Should we update `CContainer.dataOrder()` and `.find()` methods to be more Pythonic? Or just use the new utilities?

2. **Signal Integration**: Should file gleaning happen automatically on `statusChanged` signal, or remain manual?

3. **Backward Compatibility**: Do we need to maintain the old `CCP4i2DjangoDbHandler` API for existing code?

4. **Metadata Validation**: Should we add schema validation for file objects to ensure they have required qualifiers?

5. **Performance**: Should we batch database operations for better performance with many files?

## Usage Examples

### Example 1: Simple Job Execution with Database Tracking

```python
import uuid
from server.ccp4x.db.async_db_handler import AsyncDatabaseHandler
from ccp4i2.core.CCP4PluginScript import CPluginScript

async def run_ctruncate_job(project_uuid: uuid.UUID):
    """Run a ctruncate job with automatic database tracking"""

    # Create plugin
    plugin = CPluginScript(taskName="ctruncate", name="my_ctruncate_job")
    plugin.inputData.HKLIN.set("/path/to/input.mtz")
    plugin.inputData.COLANO.set("/*/*/[IMEAN,SIGIMEAN]")

    # Create database handler
    handler = AsyncDatabaseHandler(project_uuid=project_uuid)

    # Execute with automatic tracking
    async with handler.track_job(plugin):
        # Job created in database
        # Status automatically updated as plugin runs
        # Files automatically gleaned on completion
        result = await plugin.execute()

    return result
```

### Example 2: Manual File Gleaning

```python
from server.ccp4x.lib.cdata_utils import find_all_files, extract_file_metadata

async def glean_output_files(job_uuid, output_container, db_handler):
    """Manually glean output files from a container"""

    # Find all files using modern hierarchical traversal
    output_files = find_all_files(output_container)

    for file_obj in output_files:
        # Check if file exists on disk
        if not file_obj.exists():
            continue

        # Extract metadata using CData system
        metadata = extract_file_metadata(file_obj)

        # Register in database
        file_record = await db_handler.register_output_file(
            job_uuid=job_uuid,
            file_path=Path(str(file_obj)),
            file_type=metadata['file_type'],
            param_name=metadata['name'],
            sub_type=metadata.get('sub_type'),
            content_flag=metadata.get('content_flag'),
            annotation=metadata.get('annotation', ''),
        )

        # Link back to container
        file_obj.dbFileId.set(str(file_record.uuid))
```

### Example 3: Nested Job Execution

```python
async def run_pipeline_with_db(project_uuid):
    """Run a pipeline with nested jobs, all tracked in database"""

    handler = AsyncDatabaseHandler(project_uuid=project_uuid)

    # Parent job
    parent_plugin = CPluginScript(taskName="copycell", name="pipeline_parent")

    async with handler.track_job(parent_plugin):
        # Child job 1
        ctruncate = parent_plugin.makePluginObject(
            taskName="ctruncate",
            name="truncate_step"
        )

        async with handler.track_job(ctruncate):
            await ctruncate.execute()
            # Job number: "1.1"
            # Automatically linked to parent in database

        # Child job 2
        refmac = parent_plugin.makePluginObject(
            taskName="refmac",
            name="refine_step"
        )

        async with handler.track_job(refmac):
            await refmac.execute()
            # Job number: "1.2"

        await parent_plugin.execute()
        # Parent job number: "1"
```

## Conclusion

The new CData system provides **significantly more robust** introspection capabilities than the legacy CContainer API. By leveraging:

- **Metadata registry** for type-safe access
- **Hierarchical object system** for clean traversal
- **Value state tracking** for precise state management
- **Async/await** for non-blocking I/O
- **Signal system** for automatic status tracking

We can create a **cleaner, safer, and more maintainable** database integration that eliminates:

- ❌ String-based attribute access
- ❌ Fragile `getattr()` with fallbacks
- ❌ Manual type conversions with try/except
- ❌ Regex parsing of object paths
- ❌ Qt event loop dependency

The new `AsyncDatabaseHandler` + `cdata_utils` provides everything needed for modern database-backed execution. The next step is to create async versions of `import_files` and `glean_files` using these utilities.
