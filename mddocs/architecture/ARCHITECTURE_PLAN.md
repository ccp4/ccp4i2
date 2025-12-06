# CCP4i2 REST API & CLI Integration Architecture Plan

## Executive Summary

This document outlines a plan to create a **unified library architecture** that serves both:
1. **Django REST Framework API endpoints** (ViewSets)
2. **Django management commands** (CLI tools)

The goal is to eliminate code duplication and create reusable business logic that can be invoked from either interface.

---

## Current State Analysis

### What Exists

#### 1. **Django REST API** (`server/ccp4x/api/`)
Currently implemented ViewSets:
- `JobViewSet` - 30+ custom action endpoints for job management
- `ProjectViewSet` - Project CRUD and related operations
- `FileViewSet` - File operations (download, digest, preview)
- `FileImportViewSet` - File import functionality
- `FileUseViewSet` - File usage tracking
- `ProjectExportViewSet` - Project export operations
- `ProjectTagViewSet` - Tag management

#### 2. **Job Utilities Library** (`server/ccp4x/lib/job_utils/`)
45 utility modules providing core functionality:

**Job Lifecycle:**
- `create_job.py` - Create new jobs
- `clone_job.py` - Clone existing jobs
- `run_job.py` - Execute jobs
- `context_dependent_run.py` - Environment-aware job execution
- `find_dependent_jobs.py` - Dependency tracking and deletion

**Parameter Management:**
- `set_parameter.py` - Set job parameters
- `upload_file_param.py` - Upload and set file parameters
- `set_input_by_context_job.py` - Auto-configure from context
- `get_job_container.py` - Get job CContainer
- `save_params_for_job.py` - Persist parameters

**Validation & Analysis:**
- `validate_container.py` - Validate job parameters
- `digest_file.py` - Analyze file contents
- `ccp4i2_report.py` - Generate job reports
- `get_what_next.py` - Suggest next workflow steps

**File Operations:**
- `import_files.py` - Import files into jobs
- `export_job_file.py` - Export job files
- `glean_job_files.py` - Discover job files
- `patch_output_file_paths.py` - Update file paths

**Container/Object Operations:**
- `json_for_job_container.py` - Serialize containers to JSON
- `value_dict_for_object.py` - Extract object values
- `find_objects.py` - Navigate object hierarchies
- `object_method.py` - Invoke methods on objects

**Project Operations:**
- `list_project.py` - List project contents
- `create_task.py` - Create tasks in projects
- `get_task_tree.py` - Get available task types

**Utilities:**
- `i2run_for_job.py` - Generate CLI commands
- `preview_file.py` / `preview_job.py` - Launch external viewers
- `detect_file_type.py` - Identify file formats
- `mtz_as_dict.py` - Parse MTZ files
- `json_encoder.py` - Custom JSON encoding

#### 3. **Django Management Commands** (`server/ccp4x/db/management/commands/`)
Currently implemented CLI commands:
- `create_project.py` - Create new projects
- `list_projects.py` - List all projects
- `list_jobs.py` - List jobs for a project (NEW)
- `create_job.py` - Create new jobs
- `clone_job.py` - Clone jobs
- `run_job.py` - Run jobs
- `i2run.py` - Run i2 commands
- `set_job_status.py` - Update job status
- `import_ccp4_project_zip.py` - Import projects
- `export_project.py` - Export projects
- `preview_file.py` - Preview files

#### 4. **CData System** (Pure Python, No Django)
- **Core Classes**: `CContainer`, `CDataFile`, `CData` hierarchy
- **Plugin System**: 148 plugins via lazy-loading registry
- **Task Manager**: `CTaskManager` for plugin orchestration
- **Signal System**: Qt-free event handling

---

## Gap Analysis

### What's Missing

#### A. Library-Level Operations

Many API endpoints contain business logic directly in ViewSet methods. We need to extract this into reusable library functions:

**Missing from `job_utils`:**
1. **Job Report Operations**
   - `get_job_params_xml(job)` - Extract params.xml
   - `get_job_report_xml(job)` - Generate/retrieve report
   - `get_job_diagnostic_xml(job)` - Get diagnostics
   - `get_job_def_xml(job)` - Get task definition

2. **Job File Operations**
   - `get_job_files(job)` - List all job files
   - `export_job_archive(job, output_path)` - Export job as ZIP
   - `get_export_job_file_menu(job)` - Get exportable files menu

3. **Job Validation & Status**
   - `validate_job_params(job)` - Full validation with error report
   - `get_job_container_json(job)` - Get container as JSON (exists but needs normalization)

4. **Project Operations**
   - Many project operations exist in ViewSet but not in reusable form

#### B. Management Command Coverage

We have CLI commands for basic operations but many advanced operations lack CLI equivalents:

**Missing CLI Commands:**
1. `validate_job` - Validate job parameters
2. `get_job_report` - Get job report
3. `set_job_parameter` - Set individual parameters
4. `upload_job_file` - Upload file to job
5. `export_job` - Export single job
6. `digest_file` - Analyze file contents
7. `get_what_next` - Get workflow suggestions
8. `set_context_job` - Set context for auto-configuration

#### C. Consistent Error Handling

Currently, error handling varies:
- Some functions return dicts with `{"status": "Success/Failed", "reason": ...}`
- Some raise exceptions
- Some return None or empty values

We need a **consistent error handling strategy**.

---

## Proposed Architecture

### Design Principles

1. **Separation of Concerns**:
   - **Library Layer**: Pure business logic (no Django REST knowledge)
   - **API Layer**: Thin wrappers that handle HTTP concerns
   - **CLI Layer**: Thin wrappers that handle command-line concerns

2. **Single Source of Truth**:
   - All business logic lives in `server/ccp4x/lib/`
   - API endpoints call library functions
   - Management commands call library functions

3. **Consistent Interfaces**:
   - Library functions accept Django model instances or primitives
   - Library functions return structured dictionaries or raise well-defined exceptions
   - Both API and CLI layers handle formatting and presentation

4. **Type Safety**:
   - Use type hints throughout
   - Define clear input/output contracts

### Directory Structure

```
server/ccp4x/lib/
├── job_utils/           # Job-related operations (EXISTING, enhance)
│   ├── lifecycle/       # NEW: Create, clone, run, delete
│   ├── parameters/      # NEW: Set, get, validate parameters
│   ├── files/           # NEW: Import, export, digest files
│   ├── reporting/       # NEW: Reports, diagnostics, XMLs
│   └── workflow/        # NEW: What-next, context, dependencies
│
├── project_utils/       # NEW: Project operations
│   ├── lifecycle.py     # Create, delete, export projects
│   ├── listing.py       # List contents, directory operations
│   └── import_export.py # Import/export operations
│
├── file_utils/          # NEW: File operations
│   ├── detection.py     # File type detection
│   ├── digest.py        # File content analysis
│   └── preview.py       # External viewer launching
│
├── container_utils/     # NEW: CContainer/CData operations
│   ├── serialization.py # JSON/XML serialization
│   ├── navigation.py    # Object path navigation
│   └── validation.py    # Container validation
│
├── response/            # NEW: Standardized response handling
│   ├── result.py        # Result/Error types
│   └── exceptions.py    # Custom exceptions
│
└── utils/               # NEW: Cross-cutting concerns
    ├── cdata_utils.py   # CData manipulation helpers (EXISTING, enhance)
    └── json_encoder.py  # Custom encoders (EXISTING)
```

### Proposed Library Organization

#### 1. **Refactor Existing `job_utils/`**

Reorganize the 45 existing files into logical subdirectories:

```python
# job_utils/lifecycle/
from .create_job import create_job
from .clone_job import clone_job
from .run_job import run_job, run_job_local, run_job_async
from .delete_job import delete_job_and_dependents

# job_utils/parameters/
from .get_container import get_job_container
from .set_parameter import set_parameter
from .upload_file import upload_file_param
from .set_context import set_input_by_context_job
from .save_params import save_params_for_job

# job_utils/files/
from .import_files import import_files
from .export_file import export_job_file
from .glean_files import glean_job_files
from .digest import digest_param_file

# job_utils/reporting/
from .params_xml import get_job_params_xml      # NEW
from .report_xml import get_job_report_xml      # NEW
from .diagnostic_xml import get_job_diagnostic_xml  # NEW
from .def_xml import get_job_def_xml            # NEW
from .ccp4i2_report import generate_job_report      # EXISTING

# job_utils/workflow/
from .what_next import get_what_next
from .dependencies import find_dependent_jobs, find_dependency_tree  # NEW
from .validation import validate_job_container  # NEW wrapper

# job_utils/serialization/
from .json_container import json_for_job_container
from .value_dict import value_dict_for_object
```

#### 2. **Create New `project_utils/`**

Extract project operations from `ProjectViewSet`:

```python
# project_utils/lifecycle.py
def create_project(name, description, directory) -> models.Project:
    """Create a new CCP4i2 project."""

def delete_project(project: models.Project) -> None:
    """Delete project and all jobs/files."""

def export_project(project: models.Project, output_path: Path,
                   job_selection: set = None) -> Path:
    """Export project to ZIP archive."""

# project_utils/listing.py
def list_project_contents(project: models.Project) -> dict:
    """List project directory contents."""

def get_project_jobs(project: models.Project, filters: dict = None) -> list:
    """Get jobs for project with optional filters."""

def get_project_files(project: models.Project) -> list:
    """Get all files in project."""

# project_utils/import_export.py
def import_project_zip(zip_path: Path) -> models.Project:
    """Import project from ZIP file."""
```

#### 3. **Create Standardized Response Types**

```python
# lib/response/result.py
from typing import TypeVar, Generic, Optional
from dataclasses import dataclass

T = TypeVar('T')

@dataclass
class Result(Generic[T]):
    """Standardized result type for library operations."""
    success: bool
    data: Optional[T] = None
    error: Optional[str] = None
    error_details: Optional[dict] = None

    @classmethod
    def ok(cls, data: T) -> 'Result[T]':
        return cls(success=True, data=data)

    @classmethod
    def fail(cls, error: str, details: dict = None) -> 'Result[T]':
        return cls(success=False, error=error, error_details=details)

    def to_dict(self) -> dict:
        """Convert to dictionary for API/CLI responses."""
        if self.success:
            return {"status": "Success", "data": self.data}
        else:
            result = {"status": "Failed", "reason": self.error}
            if self.error_details:
                result["details"] = self.error_details
            return result


# lib/response/exceptions.py
class CCP4OperationError(Exception):
    """Base exception for CCP4 operations."""
    def __init__(self, message: str, details: dict = None):
        self.message = message
        self.details = details or {}
        super().__init__(message)

class JobNotFoundError(CCP4OperationError):
    """Job not found."""

class ValidationError(CCP4OperationError):
    """Parameter validation failed."""

class FileOperationError(CCP4OperationError):
    """File operation failed."""
```

---

## Implementation Strategy

### Phase 1: Foundation (Week 1)

**Goal**: Establish architecture and refactor core utilities

1. **Create new directory structure**
   ```bash
   mkdir -p server/ccp4x/lib/{project_utils,file_utils,container_utils,response}
   ```

2. **Implement response types**
   - Create `lib/response/result.py`
   - Create `lib/response/exceptions.py`
   - Write unit tests

3. **Refactor 5 critical job_utils functions**
   - Choose: `create_job`, `run_job`, `clone_job`, `set_parameter`, `get_job_container`
   - Update to use `Result` type
   - Add comprehensive type hints
   - Write unit tests
   - Update callers (ViewSets)

**Deliverables**:
- Response type system implemented
- 5 refactored utilities with tests
- Updated ViewSet endpoints using new utilities
- Documentation for new patterns

### Phase 2: Job Operations Library (Week 2)

**Goal**: Complete job operations library

1. **Reorganize existing job_utils**
   - Move files into subdirectories (lifecycle/, parameters/, files/, etc.)
   - Update imports across codebase
   - Ensure all tests pass

2. **Create missing job operation functions**
   - `get_job_params_xml()`
   - `get_job_report_xml()`
   - `get_job_diagnostic_xml()`
   - `get_job_def_xml()`
   - `validate_job_container()` (wrapper around existing)
   - `export_job_archive()`

3. **Refactor remaining ViewSet methods**
   - Extract business logic to library
   - ViewSet methods become thin wrappers
   - Update error handling to use Result/exceptions

**Deliverables**:
- Organized job_utils with subdirectories
- All missing job operations implemented
- JobViewSet fully refactored to use library
- Comprehensive test coverage

### Phase 3: Project & File Operations (Week 3)

**Goal**: Implement project and file utilities

1. **Create project_utils library**
   - `lifecycle.py` - create, delete, export
   - `listing.py` - list contents, jobs, files
   - `import_export.py` - import/export operations

2. **Create file_utils library**
   - `detection.py` - file type detection
   - `digest.py` - consolidate digest operations
   - `preview.py` - external viewer launching

3. **Refactor ProjectViewSet and FileViewSet**
   - Extract all business logic
   - Use new libraries
   - Standardize error handling

**Deliverables**:
- project_utils library complete
- file_utils library complete
- ProjectViewSet and FileViewSet refactored
- Tests for all new libraries

### Phase 4: Management Commands (Week 4)

**Goal**: Create CLI equivalents for all major operations

1. **Implement missing management commands**
   ```
   - validate_job - Validate job parameters
   - get_job_report - Extract job report
   - set_job_parameter - Set parameter value
   - upload_job_file - Upload file parameter
   - export_job - Export single job
   - digest_file - Analyze file
   - get_what_next - Get workflow suggestions
   - set_context_job - Set job context
   - run_job_async - Run job asynchronously
   ```

2. **Each command follows template**:
   ```python
   from django.core.management.base import BaseCommand
   from ccp4x.lib.job_utils.lifecycle import run_job
   from ccp4x.lib.response.result import Result

   class Command(BaseCommand):
       help = "Run a CCP4i2 job"

       def add_arguments(self, parser):
           parser.add_argument('job_id', type=str)

       def handle(self, *args, **options):
           job_id = options['job_id']

           # Call library function
           result: Result = run_job(job_id)

           # Handle result
           if result.success:
               self.stdout.write(self.style.SUCCESS(
                   f"Job {job_id} started successfully"
               ))
           else:
               raise CommandError(result.error)
   ```

3. **Standardize command interface**
   - Consistent argument naming
   - JSON output option (--json)
   - Verbose mode (--verbose)
   - Error handling pattern

**Deliverables**:
- 9+ new management commands
- Consistent CLI interface
- Documentation for all commands
- Integration tests

### Phase 5: Documentation & Testing (Week 5)

**Goal**: Comprehensive documentation and test coverage

1. **API Documentation**
   - Document all library functions
   - Add usage examples
   - Create architecture diagrams

2. **Testing**
   - Unit tests for all libraries (target: 90%+ coverage)
   - Integration tests for API endpoints
   - Integration tests for management commands
   - End-to-end workflow tests

3. **User Documentation**
   - Update CLAUDE.md with new architecture
   - Create CLI command reference
   - Create API endpoint reference
   - Add migration guide for existing code

**Deliverables**:
- Comprehensive test suite
- Full API documentation
- CLI reference guide
- Architecture documentation

---

## Example: Before & After

### Before: Direct ViewSet Implementation

```python
# JobViewSet.py
@action(detail=True, methods=["get"])
def params_xml(self, request, pk=None):
    try:
        the_job = models.Job.objects.get(id=pk)
        params_path = the_job.directory / "params.xml"
        fallback_params_path = the_job.directory / "input_params.xml"
        if the_job.status in [models.Job.Status.UNKNOWN, models.Job.Status.PENDING]:
            params_path = the_job.directory / "input_params.xml"
            fallback_params_path = the_job.directory / "params.xml"
        with open(params_path, "r", encoding="UTF-8") as params_xml_file:
            params_xml = params_xml_file.read()
        return Response({"status": "Success", "xml": params_xml})
    except (ValueError, models.Job.DoesNotExist) as err:
        return Response({"status": "Failed", "reason": str(err)})
    except FileNotFoundError as err:
        try:
            with open(fallback_params_path, "r") as params_xml_file:
                params_xml = params_xml_file.read()
            return Response({"status": "Success", "xml": params_xml})
        except FileNotFoundError as err1:
            return Response({"status": "Failed", "reason": str(err1)})
```

### After: Library-Based Implementation

```python
# lib/job_utils/reporting/params_xml.py
from typing import Optional
from pathlib import Path
from ccp4x.db import models
from ccp4x.lib.response.result import Result
from ccp4x.lib.response.exceptions import JobNotFoundError, FileOperationError

def get_job_params_xml(job: models.Job) -> Result[str]:
    """
    Get job parameters XML.

    Args:
        job: Job model instance

    Returns:
        Result containing XML string or error

    Raises:
        JobNotFoundError: If job doesn't exist
        FileOperationError: If no params file found
    """
    # Determine which file to check first based on status
    if job.status in [models.Job.Status.UNKNOWN, models.Job.Status.PENDING]:
        primary_path = job.directory / "input_params.xml"
        fallback_path = job.directory / "params.xml"
    else:
        primary_path = job.directory / "params.xml"
        fallback_path = job.directory / "input_params.xml"

    # Try primary path
    try:
        with open(primary_path, "r", encoding="UTF-8") as f:
            return Result.ok(f.read())
    except FileNotFoundError:
        pass

    # Try fallback path
    try:
        with open(fallback_path, "r", encoding="UTF-8") as f:
            return Result.ok(f.read())
    except FileNotFoundError:
        return Result.fail(
            f"No params file found for job {job.uuid}",
            details={
                "tried_paths": [str(primary_path), str(fallback_path)],
                "job_status": job.status
            }
        )


# api/JobViewSet.py (simplified)
from ccp4x.lib.job_utils.reporting import get_job_params_xml

@action(detail=True, methods=["get"])
def params_xml(self, request, pk=None):
    """Get job parameters XML."""
    try:
        job = models.Job.objects.get(id=pk)
        result = get_job_params_xml(job)

        if result.success:
            return Response({"status": "Success", "xml": result.data})
        else:
            return Response(result.to_dict(), status=404)

    except models.Job.DoesNotExist:
        return Response(
            {"status": "Failed", "reason": f"Job {pk} not found"},
            status=404
        )


# db/management/commands/get_job_params.py (NEW)
from django.core.management.base import BaseCommand, CommandError
from ccp4x.db import models
from ccp4x.lib.job_utils.reporting import get_job_params_xml

class Command(BaseCommand):
    help = "Get job parameters XML"

    def add_arguments(self, parser):
        parser.add_argument('job_id', help='Job UUID or database ID')
        parser.add_argument('--output', '-o', help='Output file path')

    def handle(self, *args, **options):
        job_id = options['job_id']

        # Find job
        try:
            job = models.Job.objects.get(uuid=job_id)
        except models.Job.DoesNotExist:
            try:
                job = models.Job.objects.get(id=int(job_id))
            except (models.Job.DoesNotExist, ValueError):
                raise CommandError(f"Job not found: {job_id}")

        # Get params XML
        result = get_job_params_xml(job)

        if result.success:
            output_path = options.get('output')
            if output_path:
                with open(output_path, 'w') as f:
                    f.write(result.data)
                self.stdout.write(self.style.SUCCESS(
                    f"Params XML written to {output_path}"
                ))
            else:
                self.stdout.write(result.data)
        else:
            raise CommandError(result.error)
```

---

## Benefits of This Approach

### 1. **Code Reuse**
- Business logic written once
- Used by both API and CLI
- Easier to maintain

### 2. **Testability**
- Library functions are pure business logic
- Easy to unit test without Django REST
- Can test CLI commands independently

### 3. **Consistency**
- Same error handling everywhere
- Same validation logic
- Same business rules

### 4. **Flexibility**
- Easy to add new interfaces (GraphQL, gRPC, etc.)
- Library functions can be used in scripts
- Can build composite operations easily

### 5. **Maintainability**
- Clear separation of concerns
- Easy to locate functionality
- Refactoring is safer

### 6. **Documentation**
- Library functions are self-documenting
- Type hints provide clear contracts
- Easy to generate API docs

---

## Testing Strategy

### Unit Tests (Library Layer)
```python
# tests/lib/job_utils/test_params_xml.py
def test_get_job_params_xml_finished_job(mock_job):
    mock_job.status = models.Job.Status.FINISHED
    mock_job.directory = Path("/tmp/test_job")

    # Create params.xml
    params_path = mock_job.directory / "params.xml"
    params_path.write_text("<params>test</params>")

    result = get_job_params_xml(mock_job)

    assert result.success
    assert "<params>test</params>" in result.data
```

### Integration Tests (API Layer)
```python
# tests/api/test_job_viewset.py
def test_params_xml_endpoint(client, test_job):
    response = client.get(f'/api/jobs/{test_job.id}/params_xml/')

    assert response.status_code == 200
    assert response.json()['status'] == 'Success'
    assert 'xml' in response.json()
```

### Integration Tests (CLI Layer)
```python
# tests/management/test_get_job_params.py
def test_get_job_params_command(test_job):
    out = StringIO()
    call_command('get_job_params', str(test_job.uuid), stdout=out)

    output = out.getvalue()
    assert '<params>' in output
```

---

## Migration Strategy

### For Existing Code

1. **Gradual Migration**
   - Start with new functionality
   - Refactor existing endpoints one at a time
   - Keep old code until new is proven

2. **Backward Compatibility**
   - Existing API contracts unchanged
   - Old utility functions deprecated, not removed
   - Migration period of 2-3 versions

3. **Testing During Migration**
   - Run both old and new implementations
   - Compare results
   - Gradually increase confidence

---

## Success Metrics

### Quantitative
- **Code Reuse**: 80%+ of business logic in libraries
- **Test Coverage**: 90%+ for library layer
- **API/CLI Parity**: 95%+ of operations available in both
- **Lines of Code**: Reduce ViewSet code by 60%+

### Qualitative
- Clear separation between layers
- Easy to add new operations
- Simple to create new interfaces
- Well-documented architecture

---

## Risks & Mitigation

### Risk 1: Breaking Existing API Clients
**Mitigation**:
- Maintain exact API response formats
- Version APIs if needed
- Comprehensive integration tests

### Risk 2: Performance Overhead
**Mitigation**:
- Profile critical paths
- Optimize library functions
- Use caching where appropriate

### Risk 3: Increased Complexity
**Mitigation**:
- Clear documentation
- Consistent patterns
- Code review process

### Risk 4: Migration Takes Longer Than Expected
**Mitigation**:
- Phased approach
- Prioritize high-value operations
- Can ship incrementally

---

## Next Steps

1. **Review this plan** with team
2. **Prioritize operations** for Phase 1
3. **Set up project tracking** (GitHub issues, milestones)
4. **Begin Phase 1 implementation**

---

## Appendix A: Complete Endpoint Inventory

### JobViewSet Endpoints (30)
| Endpoint | Method | Current Location | Target Library |
|----------|--------|------------------|----------------|
| `/what_next/` | GET | `get_what_next()` | ✅ `job_utils/workflow/` |
| `/set_context_job/` | POST | Inline | `job_utils/parameters/` |
| `/object_method/` | POST | `object_method()` | ✅ `job_utils/container/` |
| `/params_xml/` | GET | Inline | `job_utils/reporting/` |
| `/report_xml/` | GET | Inline + `generate_job_report()` | `job_utils/reporting/` |
| `/dependent_jobs/` | GET | `find_dependent_jobs()` | ✅ `job_utils/workflow/` |
| `/clone/` | POST | `clone_job()` | ✅ `job_utils/lifecycle/` |
| `/run/` | POST | `context_dependent_run()` | ✅ `job_utils/lifecycle/` |
| `/run_local/` | POST | `context_dependent_run()` | ✅ `job_utils/lifecycle/` |
| `/container/` | GET | `json_for_job_container()` | ✅ `job_utils/serialization/` |
| `/diagnostic_xml/` | GET | Inline | `job_utils/reporting/` |
| `/digest/` | GET | `digest_param_file()` | ✅ `job_utils/files/` |
| `/i2run_command/` | GET | `i2run_for_job()` | ✅ `job_utils/utils/` |
| `/digest_param_file/` | GET | `digest_param_file()` | ✅ `job_utils/files/` |
| `/def_xml/` | GET | Inline | `job_utils/reporting/` |
| `/validation/` | GET | `validate_container()` | ✅ `job_utils/validation/` |
| `/set_parameter/` | POST | `set_parameter()` | ✅ `job_utils/parameters/` |
| `/upload_file_param/` | POST | `upload_file_param()` | ✅ `job_utils/parameters/` |
| `/preview/` | POST | `preview_job()` | ✅ `job_utils/utils/` |
| `/files/` | GET | Inline | `job_utils/files/` |
| `/export_job/` | GET | `export_project_to_zip()` | `job_utils/files/` |
| `/export_job_file_menu/` | GET | Inline | `job_utils/files/` |
| `/export_job_file/` | GET | `export_job_file()` | ✅ `job_utils/files/` |

### ProjectViewSet Endpoints (15)
| Endpoint | Method | Current Location | Target Library |
|----------|--------|------------------|----------------|
| `/import_project/` | POST | Inline | `project_utils/import_export/` |
| `/files/` | GET | Inline | `project_utils/listing/` |
| `/file_uses/` | GET | Inline | `project_utils/listing/` |
| `/jobs/` | GET | Inline | `project_utils/listing/` |
| `/job_float_values/` | GET | Inline | `project_utils/listing/` |
| `/job_char_values/` | GET | Inline | `project_utils/listing/` |
| `/tags/` | GET/POST | Inline | `project_utils/tags/` |
| `/directory/` | GET | `list_project()` | ✅ `project_utils/listing/` |
| `/project_file/` | GET | Inline | `project_utils/files/` |
| `/preview_file/` | POST | `preview_file()` | ✅ `project_utils/files/` |
| `/create_task/` | POST | `create_task()` | ✅ `project_utils/tasks/` |
| `/export/` | POST | Inline | `project_utils/export/` |
| `/exports/` | GET | Inline | `project_utils/export/` |
| `/remove_tag/` | DELETE | Inline | `project_utils/tags/` |

### FileViewSet Endpoints (6)
| Endpoint | Method | Current Location | Target Library |
|----------|--------|------------------|----------------|
| `/by_uuid/` | GET | Inline | `file_utils/lookup/` |
| `/download/` | GET | Inline | `file_utils/download/` |
| `/download_by_uuid/` | GET | Inline | `file_utils/download/` |
| `/digest/` | GET | `digest_file()` | ✅ `file_utils/digest/` |
| `/digest_by_uuid/` | GET | `digest_file()` | ✅ `file_utils/digest/` |

**Legend:**
- ✅ = Function already exists in `job_utils/`
- Blank = Needs to be created

---

## Appendix B: Priority Matrix

### High Priority (Phase 1-2)
Operations that are frequently used or have complex logic:
1. Job creation, running, cloning
2. Parameter setting and validation
3. File operations (import, export, digest)
4. Report generation

### Medium Priority (Phase 3)
Important but less complex operations:
1. Project listing and directory operations
2. Workflow suggestions (what_next)
3. Job dependency tracking
4. Container serialization

### Low Priority (Phase 4-5)
Nice-to-have or specialized operations:
1. External viewer launching
2. Tag management
3. Export file menus
4. Diagnostic XML

---

## Appendix C: Sample Management Command Template

```python
"""
Django management command template for CCP4i2 operations.

This template provides a consistent structure for all management commands,
ensuring uniform error handling, logging, and user interaction.
"""

from django.core.management.base import BaseCommand, CommandError
from ccp4x.db import models
from ccp4x.lib.response.result import Result
from ccp4x.lib.response.exceptions import CCP4OperationError
import json
import logging

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    """
    [Command description]

    Usage:
        python manage.py <command_name> [args] [options]

    Examples:
        python manage.py <command_name> --help
        python manage.py <command_name> arg1 arg2 --json
    """

    help = "[One-line description]"
    requires_system_checks = []  # Set to [] to skip Django checks

    def add_arguments(self, parser):
        """Add command-line arguments."""
        # Required arguments
        parser.add_argument(
            'arg_name',
            type=str,
            help='Argument description'
        )

        # Optional arguments
        parser.add_argument(
            '--option',
            type=str,
            default=None,
            help='Option description'
        )

        # Common flags
        parser.add_argument(
            '--json',
            action='store_true',
            help='Output as JSON'
        )

        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Verbose output'
        )

    def handle(self, *args, **options):
        """
        Execute the command.

        Args:
            *args: Positional arguments
            **options: Command options

        Raises:
            CommandError: If operation fails
        """
        # Extract arguments
        arg_value = options['arg_name']
        json_output = options['json']
        verbose = options['verbose']

        try:
            # Find/validate models
            try:
                model_instance = models.Model.objects.get(id=arg_value)
            except models.Model.DoesNotExist:
                raise CommandError(f"Not found: {arg_value}")

            # Call library function
            from ccp4x.lib.module_name import function_name
            result: Result = function_name(model_instance)

            # Handle result
            if result.success:
                # JSON output
                if json_output:
                    self.stdout.write(json.dumps(result.to_dict(), indent=2))
                else:
                    # Human-readable output
                    self.stdout.write(self.style.SUCCESS("Operation successful"))
                    if verbose:
                        self.stdout.write(str(result.data))
            else:
                # Error output
                if json_output:
                    self.stdout.write(json.dumps(result.to_dict(), indent=2))
                else:
                    raise CommandError(result.error)

        except CCP4OperationError as e:
            # Handle known errors
            logger.exception("Operation failed", exc_info=e)
            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Failed",
                    "reason": e.message,
                    "details": e.details
                }, indent=2))
            else:
                raise CommandError(e.message)

        except Exception as e:
            # Handle unexpected errors
            logger.exception("Unexpected error", exc_info=e)
            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Failed",
                    "reason": f"Unexpected error: {str(e)}"
                }, indent=2))
            else:
                raise CommandError(f"Unexpected error: {str(e)}")
```

---

**Document Version**: 1.0
**Created**: 2025-10-31
**Author**: Architecture Planning Session
**Status**: Proposal - Awaiting Review
