# ViewSet Modernization Map

This document maps ViewSet endpoints to their modern utility functions, following the pattern established in management commands.

## Pattern: Management Commands → ViewSet Endpoints

Management commands in `server/ccp4x/db/management/commands/` demonstrate the **correct modern approach**:
- Use utilities from `ccp4x.lib.utils.*`
- Work with CPluginScript architecture
- Return `Result[T]` objects with `.success`, `.data`, `.error`
- Handle async database operations properly

## JobViewSet - Current vs Modern

| Endpoint | Current Import | Should Use | Notes |
|----------|---------------|------------|-------|
| `clone/` | ✅ `lib.utils.jobs.clone` | Already modern | Fixed! |
| `set_parameter/` | ✅ `lib.utils.parameters.set_param` | Already modern | Lines 405-430 |
| `params_xml/` | ✅ `lib.utils.jobs.reports.get_job_params_xml` | Already modern | |
| `report_xml/` | ✅ `lib.utils.jobs.reports.get_job_report_xml` | Already modern | |
| `diagnostic_xml/` | ✅ `lib.utils.jobs.reports.get_job_diagnostic_xml` | Already modern | |
| `validation/` | ✅ `lib.utils.jobs.validate.validate_job` | Already modern | |
| `run/` | ✅ `lib.utils.jobs.context_run` | Already modern | |
| `run_local/` | ✅ `lib.utils.jobs.run` | Already modern | |
| `i2run_command/` | ✅ `lib.utils.jobs.i2run.i2run_for_job` | Already modern | |
| **`container/`** | ❌ `lib.job_utils.json_for_job_container` | **`lib.utils.containers.json_for_container`** | **HIGH PRIORITY** |
| **`def_xml/`** | ❌ `lib.job_utils.load_nested_xml` | **`lib.utils.parameters.load_xml`** | **MEDIUM** |
| **`object_method/`** | ❌ `lib.job_utils.object_method` | **`lib.utils.helpers.object_method`** | **MEDIUM** |
| **`set_context_job/`** | ❌ `lib.job_utils.set_input_by_context_job` | **`lib.utils.parameters.set_input_by_context`** | **MEDIUM** |
| **`upload_file_param/`** | ❌ `lib.job_utils.upload_file_param` | **`lib.utils.files.upload_param`** | **MEDIUM** |
| **`digest/`** | ❌ `lib.job_utils.digest_param_file` | **`lib.utils.files.digest`** | **HIGH PRIORITY** |
| **`digest_param_file/`** | ❌ `lib.job_utils.digest_param_file` | **`lib.utils.files.digest`** | **HIGH PRIORITY** |
| **`what_next/`** | ❌ `lib.job_utils.get_what_next` | **`lib.utils.navigation.what_next`** | **HIGH PRIORITY** |
| **`preview/`** | ❌ `lib.job_utils.preview_job` | **`lib.utils.jobs.preview`** | **LOW** |
| `dependent_jobs/` | `lib.job_utils.find_dependent_jobs` | **`lib.utils.navigation.dependencies`** | Could modernize |
| `delete` | `lib.job_utils.delete_job_and_dependents` | `lib.utils.navigation.dependencies` | Works OK |
| `export_job/` | Direct file operations | `lib.utils.jobs.export` | Check if available |
| `export_job_file_menu/` | TaskManager direct | Keep as-is | Plugin metadata |

## FileViewSet - Current vs Modern

| Endpoint | Current Import | Should Use | Notes |
|----------|---------------|------------|-------|
| **`digest/`** | ❌ `lib.job_utils.digest_file` | **`lib.utils.files.digest.digest_file`** | **HIGH PRIORITY** |
| **`digest_by_uuid/`** | ❌ `lib.job_utils.digest_file` | **`lib.utils.files.digest.digest_file`** | **HIGH PRIORITY** |
| **`preview/`** | ❌ `lib.job_utils.preview_file` | **`lib.utils.files.preview.preview_file`** | **LOW** |
| All others | Direct Django ORM | Keep as-is | Standard CRUD |

## ProjectViewSet - Current vs Modern

| Endpoint | Current Import | Should Use | Notes |
|----------|---------------|------------|-------|
| **`directory/`** | ❌ `lib.job_utils.list_project` | **`lib.utils.navigation.list_project`** | **MEDIUM** |
| **`create_task/`** | ❌ `lib.job_utils.create_task` | **`lib.utils.tasks.create.create_task`** | **HIGH PRIORITY** |
| **`preview_file/`** | ❌ `lib.job_utils.preview_file` | **`lib.utils.files.preview.preview_file`** | **LOW** |
| `destroy` | `lib.job_utils.delete_job_and_dependents` | `lib.utils.navigation.dependencies` | Works OK |
| All others | Direct Django ORM | Keep as-is | Standard operations |

## Modern Utility Structure

```
lib/utils/
├── containers/          # Container operations
│   ├── get_container.py       # Get job/plugin container
│   ├── json_for_container.py  # Serialize container to JSON
│   ├── find_objects.py        # Navigate object paths
│   ├── validate.py            # Validate container
│   └── remove_defaults.py     # Clean up defaults
│
├── jobs/               # Job lifecycle operations
│   ├── clone.py              # Clone job (Result pattern)
│   ├── create.py             # Create new job
│   ├── execute.py            # Execute job
│   ├── run.py                # Run job locally
│   ├── context_run.py        # Context-aware run
│   ├── reports.py            # Get params/report/diagnostic XML
│   ├── validate.py           # Validate job
│   ├── i2run.py              # Generate i2run command
│   ├── preview.py            # Preview job
│   ├── export.py             # Export job
│   └── get_container.py      # Get job container
│
├── files/              # File operations
│   ├── digest.py             # Digest file contents
│   ├── preview.py            # Preview with viewer
│   ├── upload_param.py       # Upload file parameter
│   ├── import_files.py       # Import files
│   ├── glean_files.py        # Glean output files
│   ├── export.py             # Export files
│   ├── detect_type.py        # Detect file type
│   └── patch_paths.py        # Patch file paths
│
├── parameters/         # Parameter management
│   ├── set_param.py          # Set parameter (Result pattern)
│   ├── set_parameter.py      # Alternative setter
│   ├── load_xml.py           # Load task def XML
│   ├── save_params.py        # Save params to XML
│   ├── set_input_by_context.py  # Set input from context
│   └── unset_output_data.py  # Clear output data
│
├── navigation/         # Project/job navigation
│   ├── what_next.py          # Get workflow suggestions
│   ├── dependencies.py       # Find dependent jobs
│   ├── list_project.py       # List project contents
│   └── task_tree.py          # Get task hierarchy
│
├── helpers/            # Helper utilities
│   ├── object_method.py      # Call methods on objects
│   └── terminal.py           # Terminal operations
│
└── plugins/            # Plugin management
    ├── get_plugin.py         # Get plugin instance
    └── plugin_context.py     # Plugin context manager
```

## Migration Priority

### HIGH PRIORITY (Core data access - breaks frequently)
1. **`container/` endpoint**
   - Current: `job_utils.json_for_job_container`
   - Modern: `lib.utils.containers.json_for_container`
   - Impact: Core data serialization

2. **`digest/` and `digest_param_file/` endpoints**
   - Current: `job_utils.digest_param_file`, `job_utils.digest_file`
   - Modern: `lib.utils.files.digest`
   - Impact: File introspection

3. **`create_task/` endpoint**
   - Current: `job_utils.create_task`
   - Modern: `lib.utils.tasks.create.create_task`
   - Impact: Job creation

4. **`what_next/` endpoint**
   - Current: `job_utils.get_what_next`
   - Modern: `lib.utils.navigation.what_next`
   - Impact: Workflow suggestions

### MEDIUM PRIORITY (Parameter management)
5. **`upload_file_param/` endpoint**
   - Current: `job_utils.upload_file_param`
   - Modern: `lib.utils.files.upload_param`

6. **`object_method/` endpoint**
   - Current: `job_utils.object_method`
   - Modern: `lib.utils.helpers.object_method`

7. **`set_context_job/` endpoint**
   - Current: `job_utils.set_input_by_context_job`
   - Modern: `lib.utils.parameters.set_input_by_context`

8. **`def_xml/` endpoint**
   - Current: `job_utils.load_nested_xml`
   - Modern: `lib.utils.parameters.load_xml`

9. **`directory/` endpoint**
   - Current: `job_utils.list_project`
   - Modern: `lib.utils.navigation.list_project`

### LOW PRIORITY (External tools - work OK)
10. **Preview endpoints** (`preview/`, `preview_file/`, `preview_job/`)
    - Work with external applications
    - Low risk

## Update Pattern

For each endpoint, follow this pattern:

### Before (Legacy):
```python
from ..lib.job_utils.json_for_job_container import json_for_job_container

@action(detail=True, methods=["get"])
def container(self, request, pk=None):
    job = models.Job.objects.get(id=pk)
    result_dict = json_for_job_container(job)  # Returns dict directly
    return JsonResponse({"status": "Success", "result": result_dict})
```

### After (Modern):
```python
from ..lib.utils.containers.json_for_container import json_for_container
from ..lib.utils.jobs.get_container import get_job_container

@action(detail=True, methods=["get"])
def container(self, request, pk=None):
    job = models.Job.objects.get(id=pk)
    result = get_job_container(job)  # Returns Result[CContainer]

    if result.success:
        container = result.data
        json_result = json_for_container(container)  # Serialize
        return Response({"status": "Success", "result": json_result})
    else:
        return Response(
            {"status": "Failed", "reason": result.error},
            status=400
        )
```

## Key Differences

### Result Pattern
Modern utilities return `Result[T]` objects:
```python
class Result:
    success: bool
    data: T | None      # Available if success=True
    error: str | None   # Available if success=False

    def to_dict(self) -> dict:
        """Convert to JSON-serializable dict"""
```

### Error Handling
- **Legacy**: Raises exceptions, returns None, or returns success/failure dict inconsistently
- **Modern**: Always returns `Result[T]`, consistent error handling

### Database Sync
- **Legacy**: May not handle CData → DB synchronization
- **Modern**: Properly syncs CData changes to database

### Type Safety
- **Legacy**: Untyped, implicit contracts
- **Modern**: Type hints, explicit Result types

## Testing After Migration

For each endpoint updated:

1. **Run specific test**:
   ```bash
   pytest server/ccp4x/tests/api/test_viewsets_comprehensive.py::JobViewSetTests::test_job_container -xvs
   ```

2. **Verify Result pattern**:
   - Check `result.success` is True for valid inputs
   - Check `result.data` contains expected data
   - Check `result.error` is None on success

3. **Test error cases**:
   - Invalid job ID
   - Missing parameters
   - Invalid parameter paths

4. **Check database sync**:
   - Verify changes persist
   - Check async_glean_files works
   - Verify file uses created

## Next Steps

1. Update JobViewSet imports for HIGH PRIORITY endpoints
2. Update endpoint implementations to use Result pattern
3. Run tests to verify behavior unchanged
4. Repeat for MEDIUM and LOW priority
5. Remove old job_utils modules once all ViewSets updated

## Benefits

- ✅ Consistent Result pattern across all endpoints
- ✅ Proper async database synchronization
- ✅ Better error handling and reporting
- ✅ Type safety with Result[T]
- ✅ Follows CPluginScript architecture
- ✅ Management commands and API endpoints use same code
- ✅ Easier to test and maintain
