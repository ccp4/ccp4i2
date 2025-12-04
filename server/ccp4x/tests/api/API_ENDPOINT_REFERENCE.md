# API Endpoint Reference

This document catalogs all API endpoints, their purposes, and their dependencies on `job_utils` modules.

## Purpose

This reference enables systematic refactoring of `job_utils` by:
1. Documenting which endpoints depend on which legacy utilities
2. Providing test coverage to ensure behavior remains unchanged
3. Enabling safe replacement of legacy code with modern CData introspection

## FileViewSet

Base URL: `/files/`

| Endpoint | Method | Purpose | job_utils Dependencies |
|----------|--------|---------|----------------------|
| `/` | GET | List all files | None |
| `/{id}/` | GET | Retrieve file by ID | None |
| `/{uuid}/by_uuid/` | GET | Retrieve file by UUID | None |
| `/{id}/download/` | GET | Download file by ID | None |
| `/{uuid}/download_by_uuid/` | GET | Download file by UUID | None |
| `/{id}/digest/` | GET | Get file digest/summary by ID | `digest_file` |
| `/{uuid}/digest_by_uuid/` | GET | Get file digest/summary by UUID | `digest_file` |
| `/{id}/preview/` | POST | Preview file with external viewer | `preview_file` |

**Refactoring Targets:**
- `digest_file`: Should use modern CData file introspection
- `preview_file`: Launches external applications, may keep as-is

## JobViewSet

Base URL: `/jobs/`

| Endpoint | Method | Purpose | job_utils Dependencies |
|----------|--------|---------|----------------------|
| `/` | GET | List all jobs | None |
| `/{id}/` | GET | Retrieve job by ID | None |
| `/{id}/` | PUT | Update job | None |
| `/{id}/` | DELETE | Delete job and dependents | `delete_job_and_dependents` |
| `/{id}/what_next/` | GET | Get suggested next steps | `get_what_next` |
| `/{id}/set_context_job/` | POST | Set context job for inputs | `set_input_by_context_job`, `getEtree` |
| `/{id}/object_method/` | POST | Execute method on job object | `object_method`, `getEtree` |
| `/{id}/params_xml/` | GET | Get job parameters as XML | Uses unified utility (no job_utils) |
| `/{id}/report_xml/` | GET | Get job report as XML | Uses unified utility (no job_utils) |
| `/{id}/dependent_jobs/` | GET | Get dependent jobs | `find_dependent_jobs` |
| `/{id}/clone/` | POST | Clone job | `clone_job` |
| `/{id}/run/` | POST | Execute job (context-aware) | Uses unified utility (no job_utils) |
| `/{id}/run_local/` | POST | Execute job locally | Uses unified utility (no job_utils) |
| `/{id}/container/` | GET | Get job container JSON | `json_for_job_container` |
| `/{id}/diagnostic_xml/` | GET | Get diagnostic XML | Uses unified utility (no job_utils) |
| `/{id}/digest/` | GET | Digest object at path | `digest_param_file` |
| `/{id}/i2run_command/` | GET | Get i2run command line | `i2run_for_job` |
| `/{id}/digest_param_file/` | GET | Digest parameter file | `digest_param_file` |
| `/{id}/def_xml/` | GET | Get task definition XML | `load_nested_xml` |
| `/{id}/validation/` | GET | Validate job parameters | Uses unified utility (no job_utils) |
| `/{id}/set_parameter/` | POST | Set job parameter | Uses unified utility (no job_utils) |
| `/{id}/upload_file_param/` | POST | Upload file parameter | `upload_file_param`, `getEtree` |
| `/{id}/preview/` | POST | Preview job with viewer | `preview_job` |
| `/{id}/files/` | GET | Get job files | None |
| `/{id}/export_job/` | GET | Export job as ZIP | None (uses export_project_to_zip) |
| `/{id}/export_job_file_menu/` | GET | Get export file menu | None (uses TASKMANAGER) |
| `/{id}/export_job_file/` | GET | Export specific job file | `export_job_file` |

**Refactoring Targets (Priority Order):**

### High Priority (Core data access)
1. **`json_for_job_container`**: Should use modern CData serialization
   - Currently navigates container hierarchy manually
   - Should use HierarchicalObject.find_children() and CData metadata

2. **`digest_param_file`**: Should use CData file introspection
   - Currently uses legacy container traversal
   - Should use find_all_files() from cdata_utils

3. **`get_what_next`**: Should use modern workflow analysis
   - Currently uses legacy database queries
   - Should use CData metadata and FileUse relationships

### Medium Priority (Parameter management)
4. **`set_input_by_context_job`**: Already being replaced
   - Uses unified utility pattern
   - Some legacy code may remain

5. **`object_method`**: Should use HierarchicalObject.find_by_path()
   - Currently navigates container manually
   - Should use modern object_path() and find_by_path()

6. **`upload_file_param`**: Should use modern file handling
   - Currently uses legacy container methods
   - Should use AsyncDatabaseHandler + CData

### Low Priority (Utility functions)
7. **`clone_job`**: Job lifecycle operation
   - May use modern utilities internally
   - Less critical for CData refactoring

8. **`find_dependent_jobs`**: Database query utility
   - Already uses Django ORM
   - Minimal CData dependency

9. **`delete_job_and_dependents`**: Cleanup operation
   - Already uses Django ORM
   - Minimal CData dependency

10. **`i2run_for_job`**: Command generation
    - Utility function, minimal CData access
    - Low priority

11. **`load_nested_xml`**: XML processing
    - Not CData-related
    - Can remain as-is

12. **`getEtree`**: Error formatting
    - Not CData-related
    - Can remain as-is

### Keep As-Is
- **`preview_job`**: Launches external applications
- **`preview_file`**: Launches external applications
- **`export_job_file`**: File export utility

## ProjectViewSet

Base URL: `/projects/`

| Endpoint | Method | Purpose | job_utils Dependencies |
|----------|--------|---------|----------------------|
| `/` | GET | List all projects | None |
| `/{id}/` | GET | Retrieve project by ID | None |
| `/{id}/` | DELETE | Delete project | `delete_job_and_dependents` |
| `/import_project/` | POST | Import project from ZIP | None (uses management command) |
| `/{id}/files/` | GET | Get project files | None |
| `/{id}/file_uses/` | GET | Get project file uses | None |
| `/{id}/jobs/` | GET | Get project jobs | None |
| `/{id}/job_float_values/` | GET | Get KPI float values | None |
| `/{id}/job_char_values/` | GET | Get KPI char values | None |
| `/{id}/tags/` | GET | Get project tags | None |
| `/{id}/tags/` | POST | Add tag to project | None |
| `/{id}/tags/{tag_id}/` | DELETE | Remove tag from project | None |
| `/{id}/directory/` | GET | Get directory listing | `list_project` |
| `/{id}/project_file/` | GET | Get project file | None |
| `/{id}/preview_file/` | POST | Preview file with viewer | `preview_file` |
| `/{id}/create_task/` | POST | Create new task/job | `create_task` |
| `/{id}/export/` | POST | Export project | None (uses subprocess) |
| `/{id}/exports/` | GET | Get export history | None |

**Refactoring Targets:**
1. **`list_project`**: Should use modern directory traversal
2. **`create_task`**: Should use modern CData initialization
3. **`preview_file`**: Keep as-is (external app)

## Simple CRUD ViewSets

These ViewSets provide standard REST operations with no job_utils dependencies:

- **FileImportViewSet** (`/fileimports/`): CRUD for file imports
- **FileUseViewSet** (`/fileuses/`): CRUD for file uses
- **FileTypeViewSet** (`/filetypes/`): CRUD for file types
- **ProjectTagViewSet** (`/projecttags/`): CRUD for project tags
- **ProjectExportViewSet** (`/projectexports/`): CRUD + download for exports

## Testing Strategy

### Phase 1: Baseline Tests (Complete)
- ✅ Comprehensive test suite created
- ✅ All endpoints documented
- ✅ Dependencies cataloged

### Phase 2: Run Tests
```bash
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
pytest server/ccp4x/tests/api/test_viewsets_comprehensive.py -v
```

### Phase 3: Refactor job_utils
For each job_utils module:
1. Identify all API endpoints that depend on it
2. Run baseline tests to capture current behavior
3. Replace with modern CData introspection
4. Re-run tests to ensure behavior unchanged
5. Remove old job_utils code

### Phase 4: Integration Tests
After refactoring, run full test suite:
```bash
pytest server/ccp4x/tests/api/ -v
pytest tests/ -v
```

## Example Refactoring: json_for_job_container

**Current approach (job_utils):**
```python
# In job_utils/json_for_job_container.py
def json_for_job_container(job):
    container = get_job_container(job)  # Legacy container access
    # Manual traversal of container hierarchy
    return manually_serialize(container)
```

**Modern approach (CData introspection):**
```python
# In lib/utils/jobs/container.py
from core.cdata_utils import serialize_cdata_to_json

def json_for_job_container(job):
    plugin = get_job_plugin(job)
    # Use modern CData serialization
    return serialize_cdata_to_json(plugin.container)
```

## Modernization Utilities

Create these utilities in `server/ccp4x/lib/cdata_utils/`:

1. **`serialize_cdata_to_json(obj)`**: Convert CData to JSON
   - Use get_merged_metadata() for type info
   - Use find_children() for hierarchy
   - Use get_qualifier() for constraints

2. **`find_all_files_in_container(container)`**: Find all files
   - Use HierarchicalObject.find_children()
   - Filter by isinstance(obj, CDataFile)
   - Extract metadata with extract_file_metadata()

3. **`digest_cdata_file(file_obj)`**: Digest file using metadata
   - Use file_obj.mimeTypeName for type
   - Use appropriate parser (gemmi, BioPython, etc.)
   - Return structured digest

4. **`navigate_to_path(container, object_path)`**: Navigate hierarchy
   - Use HierarchicalObject.find_by_path() if available
   - Or implement dot-notation navigation
   - Return object at path

## Migration Checklist

- [ ] Run baseline tests (capture current behavior)
- [ ] Create modern CData utilities
- [ ] Refactor `json_for_job_container` ← Start here
- [ ] Refactor `digest_param_file`
- [ ] Refactor `get_what_next`
- [ ] Refactor `set_input_by_context_job`
- [ ] Refactor `object_method`
- [ ] Refactor `upload_file_param`
- [ ] Refactor `list_project`
- [ ] Refactor `create_task`
- [ ] Run full test suite
- [ ] Remove old job_utils code
- [ ] Update documentation

## Notes

- All endpoints using "unified utilities" have already been migrated
- Focus refactoring on job_utils modules that navigate CData hierarchies
- Keep preview/export utilities as-is (they launch external processes)
- Prioritize high-traffic endpoints first (clone, container, set_parameter)
