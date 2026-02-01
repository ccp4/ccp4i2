# CCP4i2 API Testing Guide

This document describes the API testing infrastructure for CCP4i2, covering test execution, prerequisites, coverage, and the underlying mechanism.

## Quick Start

```bash
# From the server/ directory
./run_test.sh ccp4i2/tests/api/ -n 8
```

## Prerequisites

### 1. CCP4 Environment

The tests require the CCP4 suite with `ccp4-python`. Before running tests:

```bash
# Option A: Source CCP4 setup directly
source /path/to/ccp4-20251105/bin/ccp4.setup-sh

# Option B: Place CCP4 in a sibling directory (auto-detected)
# ../ccp4-20251105/bin/ccp4.setup-sh

# Option C: Set CCP4_ROOT in .env file
echo "CCP4_ROOT=/path/to/ccp4-20251105" > .env
```

### 2. Python Packages

Install ccp4i2 with test dependencies:

```bash
ccp4-python -m pip install -e ".[full]"
```

Key packages required:
- `pytest` - Test framework
- `pytest-django` - Django integration for pytest
- `pytest-xdist` - Parallel test execution
- `djangorestframework` - REST API framework
- `django` - Web framework

### 3. Test Data

Many tests require pre-built project zips located at:
```
../../test101/ProjectZips/
```

These contain crystallographic test data (refmac_gamma_test_0, aimless_gamma_native_test_1, etc.). Tests are automatically skipped if this directory doesn't exist.

## Running Tests

### Basic Usage

```bash
# Run all API tests
./run_test.sh ccp4i2/tests/api/

# Run specific test file
./run_test.sh ccp4i2/tests/api/test_api.py

# Run specific test class or method
./run_test.sh ccp4i2/tests/api/test_api.py::TestCCP4i2API::test_clone
```

### Parallel Execution

Tests are designed for parallel execution using pytest-xdist:

```bash
# Run with 8 parallel workers
./run_test.sh ccp4i2/tests/api/ -n 8

# Auto-detect CPU count
./run_test.sh ccp4i2/tests/api/ -n auto

# Exclude pipeline tests (long-running crystallographic jobs)
./run_test.sh ccp4i2/tests/api/ -m "not pipeline" -n 8
```

### Verbose Output

```bash
# Show test names as they run
./run_test.sh ccp4i2/tests/api/ -v

# Show print statements and logs
./run_test.sh ccp4i2/tests/api/ -v -s
```

## Test Coverage

### Test Files

| File | Purpose |
|------|---------|
| `test_api.py` | Core API endpoints: projects, jobs, cloning, parameter setting, file upload |
| `test_viewsets_comprehensive.py` | Complete coverage of all ViewSet endpoints |
| `test_job_utils.py` | Job utility functions: clone, what_next, digest, export |
| `test_job_execution_via_api.py` | End-to-end job configuration and execution |
| `test_parameter_setting_api.py` | Parameter setting with type coercion and validation |
| `test_project_tag_api.py` | Project tag CRUD operations |

### Endpoint Coverage

The test suite covers all major API endpoints:

**ProjectViewSet** (`/api/ccp4i2/projects/`)
- List, retrieve, create, delete projects
- Get project files, jobs, tags
- Create tasks, export projects, preview files

**JobViewSet** (`/api/ccp4i2/jobs/`)
- List, retrieve, clone, delete jobs
- Get/set parameters, upload files
- Get container, params XML, report XML
- Run jobs, get what_next recommendations

**FileViewSet** (`/api/ccp4i2/files/`)
- List, retrieve, download files
- Get file digest, preview files

**Simple CRUD ViewSets**
- ProjectTagViewSet (`/projecttags/`)
- FileTypeViewSet (`/filetypes/`)
- FileImportViewSet (`/fileimports/`)
- FileUseViewSet (`/fileuses/`)
- ProjectExportViewSet (`/projectexports/`)

## Test Mechanism

### Database Isolation

Each test runs with its own isolated SQLite database to enable parallel execution:

1. **Per-test database**: The `isolated_test_db` fixture creates a unique SQLite database for each test at `~/.cache/ccp4i2-tests/<timestamp>_<test_name>/project.sqlite`

2. **WAL mode**: SQLite Write-Ahead Logging enables concurrent reads during subprocess job execution

3. **Automatic cleanup**: Passed tests have their directories removed; failed tests are preserved for debugging

### Permission Bypass

API tests bypass authentication to focus on endpoint functionality:

```python
@pytest.fixture(autouse=True)
def bypass_api_permissions(monkeypatch):
    """Patch all viewsets to allow unauthenticated access during tests."""
    from rest_framework.permissions import AllowAny
    monkeypatch.setattr(ProjectViewSet, 'permission_classes', [AllowAny])
    # ... other viewsets
```

### Test Project Setup

Tests typically import pre-built project zips:

```python
@pytest.fixture(autouse=True)
def setup(self, bypass_api_permissions, test_project_path):
    self.client = APIClient()
    test_project_path.mkdir(exist_ok=True)
    import_ccp4_project_zip(
        TEST_DATA_DIR / "refmac_gamma_test_0.ccp4_project.zip",
        relocate_path=test_project_path,
    )
```

### Key Fixtures

| Fixture | Scope | Purpose |
|---------|-------|---------|
| `bypass_api_permissions` | function | Allows unauthenticated API access |
| `test_project_path` | function | Unique directory for each test |
| `isolated_test_db` | function | Isolated SQLite database per test |
| `demo_data_dir` | session | Path to demo crystallographic data |
| `gamma_mtz`, `gamma_model_pdb`, etc. | session | Paths to specific demo files |

## Writing New Tests

### Basic Test Structure

```python
import pytest
from rest_framework.test import APIClient
from ccp4i2.db.import_i2xml import import_ccp4_project_zip
from ccp4i2.db import models

API_PREFIX = "/api/ccp4i2"
TEST_DATA_DIR = Path(__file__).parent.parent.parent.parent.parent.parent / "test101" / "ProjectZips"
SKIP_REASON = f"Test data not found: {TEST_DATA_DIR}"

@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason=SKIP_REASON)
class TestMyFeature:

    @pytest.fixture(autouse=True)
    def setup(self, bypass_api_permissions, test_project_path):
        self.client = APIClient()
        test_project_path.mkdir(exist_ok=True)
        import_ccp4_project_zip(
            TEST_DATA_DIR / "refmac_gamma_test_0.ccp4_project.zip",
            relocate_path=test_project_path,
        )

    def test_my_endpoint(self):
        response = self.client.get(f"{API_PREFIX}/projects/")
        assert response.status_code == 200
```

### Best Practices for Parallel Execution

1. **Avoid hardcoded IDs**: Use dynamic lookups instead of assuming project/job IDs
   ```python
   # Bad
   response = self.client.get(f"{API_PREFIX}/projects/1/files/")

   # Good
   project = models.Project.objects.first()
   response = self.client.get(f"{API_PREFIX}/projects/{project.id}/files/")
   ```

2. **Use flexible assertions**: API response formats may vary
   ```python
   # Bad - too brittle
   assert result == {"success": True, "data": {"value": 5}}

   # Good - check key indicators
   assert result.get("success") is True
   assert "data" in result
   ```

3. **Generate unique data**: Avoid conflicts between parallel tests
   ```python
   import uuid
   unique_name = f"Test Item {uuid.uuid4().hex[:8]}"
   ```

4. **Clone before modifying**: Can't modify completed jobs
   ```python
   job = models.Job.objects.first()
   clone_response = self.client.post(f"{API_PREFIX}/jobs/{job.id}/clone/")
   clone = clone_response.json()
   # Now modify the clone, not the original
   ```

## Debugging Failed Tests

### Preserved Project Directories

Failed tests preserve their project directories for inspection:

```bash
# List failed test directories
ls ~/.cache/ccp4i2-tests/

# Example: inspect a failed test
cd ~/.cache/ccp4i2-tests/20240101_120000_1234_api_test_my_feature/
cat project.sqlite  # Database
ls job_1/           # Job directory with input_params.xml, etc.
```

### Cleanup

```bash
# Remove all test projects
rm -rf ~/.cache/ccp4i2-tests/*

# Remove projects from a specific date
rm -rf ~/.cache/ccp4i2-tests/20240101_*
```

## Environment Variables

| Variable | Purpose |
|----------|---------|
| `CCP4I2_ROOT` | Root directory of ccp4i2 installation |
| `DJANGO_SETTINGS_MODULE` | Django settings (default: `ccp4i2.config.test_settings`) |
| `CCP4I2_DB_FILE` | Path to SQLite database (set per-test) |
| `CCP4I2_PROJECTS_DIR` | Directory for test projects |

## Markers

```bash
# Skip pipeline tests (crystallographic job execution)
./run_test.sh ccp4i2/tests/api/ -m "not pipeline"

# Run only pipeline tests
./run_test.sh ccp4i2/tests/api/ -m pipeline
```

---

## Refactoring Context

These tests enable systematic refactoring of legacy `job_utils` code by:

1. **Documenting Current Behavior**: Tests capture how endpoints currently work
2. **Ensuring Compatibility**: Tests verify behavior doesn't change during refactoring
3. **Enabling Confidence**: Comprehensive coverage means safe refactoring
4. **Tracking Progress**: Can run tests after each module replacement

### job_utils Migration Targets

For detailed refactoring guidance, see `API_ENDPOINT_REFERENCE.md`. Key modules:

| Module | Priority | Endpoints |
|--------|----------|-----------|
| `json_for_job_container` | HIGH | `/jobs/{id}/container/` |
| `digest_param_file` | HIGH | `/jobs/{id}/digest/`, `/jobs/{id}/digest_param_file/` |
| `get_what_next` | HIGH | `/jobs/{id}/what_next/` |
| `object_method` | MEDIUM | `/jobs/{id}/object_method/` |
| `upload_file_param` | MEDIUM | `/jobs/{id}/upload_file_param/` |
