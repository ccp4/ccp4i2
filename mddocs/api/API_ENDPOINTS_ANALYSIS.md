# API Endpoints Analysis & Refactoring Plan

**Date**: 2025-10-31
**Purpose**: Map API endpoints to utilities and management commands, identify refactoring opportunities

---

## Overview

The JobViewSet provides **24 custom action endpoints** plus standard REST operations. This analysis identifies:
1. Which endpoints use CPluginScript architecture (✅)
2. Which use legacy approaches (⚠️)
3. Which have corresponding management commands
4. Which should be refactored for consistency

---

## Endpoint Inventory

### Category 1: Parameter Management

| Endpoint | Method | Architecture | Utility Function | Management Command | Status |
|----------|--------|--------------|------------------|-------------------|--------|
| `set_parameter/` | POST | ⚠️ **OLD** - Uses `set_parameter()` from `job_utils` | `set_parameter()` | `set_job_parameter` | **NEEDS REFACTOR** |
| `upload_file_param/` | POST | ⚠️ **OLD** - Uses `upload_file_param()` | `upload_file_param()` | ❌ No command | **NEEDS REFACTOR** |
| `set_context_job/` | POST | ⚠️ **OLD** - Uses `set_input_by_context_job()` | `set_input_by_context_job()` | ❌ No command | **NEEDS REFACTOR** |

#### Analysis

**`set_parameter/` (Line 1111-1164)**:
```python
# CURRENT (OLD):
def set_parameter(self, request, pk=None):
    result = set_parameter(job, object_path, value)  # ← Uses legacy job_utils function
    return JsonResponse({"status": "Success", "updated_item": result})
```

**RECOMMENDATION**: Refactor to use new unified utility:
```python
# SHOULD BE:
from ccp4x.lib.utils.parameters.set_param import set_parameter as set_job_param

def set_parameter(self, request, pk=None):
    result = set_job_param(job, object_path, value)  # ← Uses CPluginScript
    if result.success:
        return JsonResponse({"status": "Success", "updated_item": result.data})
    else:
        return JsonResponse({"status": "Failed", "reason": result.error}, status=400)
```

**Benefits**:
- ✅ Proper DB synchronization via dbHandler
- ✅ Consistent Result[T] pattern
- ✅ Same logic as management command
- ✅ Better error handling

---

### Category 2: Validation

| Endpoint | Method | Architecture | Utility Function | Management Command | Status |
|----------|--------|--------------|------------------|-------------------|--------|
| `validation/` | GET | ⚠️ **PARTIAL** - Uses `get_job_container()` + `validate_container()` | `validate_container()` | `validate_job` | **NEEDS REFACTOR** |

#### Analysis

**`validation/` (Line 1049-1103)**:
```python
# CURRENT (PARTIAL - uses old get_job_container):
def validation(self, request, pk=None):
    container: CContainer = get_job_container(the_job)  # ← OLD: orphaned container
    error_etree: ET.Element = validate_container(container)
    return Response({"status": "Success", "xml": ET.tostring(error_etree)})
```

**RECOMMENDATION**: Refactor to use unified utility:
```python
# SHOULD BE:
from ccp4x.lib.utils.jobs.validate import validate_job

def validation(self, request, pk=None):
    result = validate_job(the_job)  # ← Uses CPluginScript
    if result.success:
        return Response({"status": "Success", "xml": ET.tostring(result.data)})
    else:
        return Response({"status": "Failed", "reason": result.error}, status=400)
```

**Benefits**:
- ✅ Uses CPluginScript architecture
- ✅ Consistent with management command
- ✅ Proper error reporting with new CErrorReport API

---

### Category 3: Reports & XML

| Endpoint | Method | Architecture | Utility Function | Management Command | Status |
|----------|--------|--------------|------------------|-------------------|--------|
| `params_xml/` | GET | ⚠️ **DIRECT** - Reads filesystem directly | ❌ No utility | `get_job_report --type params` | **COULD REFACTOR** |
| `report_xml/` | GET | ⚠️ **PARTIAL** - Uses `make_old_report()` | `make_old_report()` | `get_job_report --type report` | **COULD REFACTOR** |
| `diagnostic_xml/` | GET | ⚠️ **DIRECT** - Reads filesystem directly | ❌ No utility | `get_job_report --type diagnostic` | **COULD REFACTOR** |
| `def_xml/` | GET | ⚠️ **DIRECT** - Reads filesystem directly | ❌ No utility | ❌ No command | **OK AS IS** |

#### Analysis

**`params_xml/` (Line 367-418)**:
```python
# CURRENT (DIRECT READ):
def params_xml(self, request, pk=None):
    params_path = the_job.directory / "params.xml"
    with open(params_path, "r", encoding="UTF-8") as params_xml_file:
        params_xml = params_xml_file.read()
    return Response({"status": "Success", "xml": params_xml})
```

**RECOMMENDATION**: Optionally refactor to use unified utility:
```python
# COULD BE (for consistency):
from ccp4x.lib.utils.jobs.reports import get_job_params_xml

def params_xml(self, request, pk=None):
    result = get_job_params_xml(the_job)
    if result.success:
        return Response({"status": "Success", "xml": result.data})
    else:
        return Response({"status": "Failed", "reason": result.error}, status=404)
```

**Benefits**:
- ✅ Consistent error handling
- ✅ Same logic as management command
- ✅ Centralized file path logic

**Priority**: LOW - Current approach works fine, but refactoring would improve consistency

---

### Category 4: Job Lifecycle

| Endpoint | Method | Architecture | Utility Function | Management Command | Status |
|----------|--------|--------------|------------------|-------------------|--------|
| `clone/` | POST | ✅ **GOOD** - Uses `clone_job()` utility | `clone_job()` | `clone_job` | ✅ **CONSISTENT** |
| `run/` | POST | ✅ **GOOD** - Uses `run_job_context_aware()` | `run_job_context_aware()` | `execute_job` | ✅ **CONSISTENT** |
| `run_local/` | POST | ✅ **GOOD** - Uses `run_job_context_aware(force_local=True)` | `run_job_context_aware()` | `execute_job --force-local` | ✅ **CONSISTENT** |

#### Analysis

These endpoints already use centralized utilities that are shared with management commands. ✅ No refactoring needed!

**`clone/` (Line 564-612)**: ✅ Already uses `clone_job()` utility
**`run/` (Line 620-671)**: ✅ Already uses `run_job_context_aware()` utility
**`run_local/` (Line 679-741)**: ✅ Already uses same utility with `force_local=True`

---

### Category 5: File Operations

| Endpoint | Method | Architecture | Utility Function | Management Command | Status |
|----------|--------|--------------|------------------|-------------------|--------|
| `files/` | GET | ✅ **GOOD** - Uses Django ORM directly | N/A | ❌ No command | ✅ **OK AS IS** |
| `export_job/` | GET | ✅ **GOOD** - Uses `export_project_to_zip()` | `export_project_to_zip()` | `export_job` | ✅ **CONSISTENT** |
| `export_job_file/` | GET | ✅ **GOOD** - Uses `export_job_file()` utility | `export_job_file()` | ❌ No command | ✅ **OK AS IS** |
| `export_job_file_menu/` | GET | ✅ **GOOD** - Uses TaskManager directly | N/A | ❌ No command | ✅ **OK AS IS** |

#### Analysis

These endpoints already use proper utilities or direct ORM access. ✅ No refactoring needed!

---

### Category 6: Advanced Operations

| Endpoint | Method | Architecture | Utility Function | Management Command | Status |
|----------|--------|--------------|------------------|-------------------|--------|
| `object_method/` | POST | ⚠️ **OLD** - Uses `object_method()` | `object_method()` | ❌ No command | **LOW PRIORITY** |
| `container/` | GET | ⚠️ **OLD** - Uses `json_for_job_container()` | `json_for_job_container()` | ❌ No command | **LOW PRIORITY** |
| `what_next/` | GET | ✅ **GOOD** - Uses `get_what_next()` utility | `get_what_next()` | ❌ No command | ✅ **OK AS IS** |
| `dependent_jobs/` | GET | ✅ **GOOD** - Uses `find_dependent_jobs()` utility | `find_dependent_jobs()` | ❌ No command | ✅ **OK AS IS** |

#### Analysis

`object_method/` and `container/` use legacy job_utils functions, but these are specialized endpoints for dynamic object manipulation. Refactoring would provide limited benefit.

**Priority**: LOW

---

### Category 7: Utility/Helper Endpoints

| Endpoint | Method | Architecture | Utility Function | Management Command | Status |
|----------|--------|--------------|------------------|-------------------|--------|
| `digest/` | GET | ⚠️ **OLD** - Uses `digest_param_file()` | `digest_param_file()` | ❌ No command | **LOW PRIORITY** |
| `digest_param_file/` | GET | ⚠️ **OLD** - Uses `digest_param_file()` | `digest_param_file()` | ❌ No command | **LOW PRIORITY** |
| `i2run_command/` | GET | ⚠️ **OLD** - Uses `i2run_for_job()` | `i2run_for_job()` | ❌ No command | **LOW PRIORITY** |
| `preview/` | POST | ⚠️ **OLD** - Uses `preview_job()` | `preview_job()` | ❌ No command | **LOW PRIORITY** |

#### Analysis

These are helper/utility endpoints that don't modify job state. Current implementation is acceptable.

**Priority**: LOW

---

## Refactoring Priority

### 🔴 HIGH PRIORITY (Functional Improvements)

1. **`set_parameter/` endpoint** - Should use new `ccp4x.lib.utils.parameters.set_param.set_parameter()`
   - **Current**: Uses legacy `job_utils.set_parameter` (doesn't use CPluginScript properly)
   - **Should use**: `ccp4x.lib.utils.parameters.set_param.set_parameter()` (CPluginScript + dbHandler)
   - **Benefit**: Proper DB sync, consistent with management command
   - **Impact**: HIGH - This is a core data modification endpoint

2. **`validation/` endpoint** - Should use new `ccp4x.lib.utils.jobs.validate.validate_job()`
   - **Current**: Uses `get_job_container()` + `validate_container()`
   - **Should use**: `ccp4x.lib.utils.jobs.validate.validate_job()` (CPluginScript)
   - **Benefit**: Uses new CErrorReport API, consistent with management command
   - **Impact**: MEDIUM - Affects validation workflow

### 🟡 MEDIUM PRIORITY (Consistency Improvements)

3. **`upload_file_param/` endpoint** - Should have unified utility
   - **Current**: Uses legacy `job_utils.upload_file_param`
   - **Should have**: New utility in `ccp4x.lib.utils.parameters.upload_file.py`
   - **Benefit**: Consistency, proper error handling
   - **Impact**: MEDIUM - File uploads are important but less frequent

4. **`set_context_job/` endpoint** - Should have unified utility
   - **Current**: Uses legacy `job_utils.set_input_by_context_job`
   - **Should have**: New utility in `ccp4x.lib.utils.jobs.context.py`
   - **Benefit**: Consistency with other utilities
   - **Impact**: LOW-MEDIUM - Used for workflow automation

### 🟢 LOW PRIORITY (Optional Consistency)

5. **Report endpoints** (`params_xml/`, `report_xml/`, `diagnostic_xml/`)
   - **Current**: Mix of direct filesystem reads and `make_old_report()`
   - **Could use**: Unified utilities from `ccp4x.lib.utils.jobs.reports`
   - **Benefit**: Consistency, centralized logic
   - **Impact**: LOW - Current approach works well

6. **Advanced endpoints** (`object_method/`, `container/`, etc.)
   - **Current**: Various legacy utilities
   - **Benefit**: Limited - these are specialized operations
   - **Impact**: LOW - Not frequently used

---

## Architecture Mapping

### Current State

```
API Endpoint → Legacy job_utils Function → Direct container manipulation
```

### Target State (CPluginScript Architecture)

```
API Endpoint → Unified Utility Function → get_plugin_with_context() → CPluginScript + dbHandler
                    ↓
           Management Command
```

### Benefits of Target Architecture

1. **Single Source of Truth**: Both API and CLI use same utility functions
2. **Proper DB Sync**: dbHandler ensures database stays synchronized
3. **Type Safety**: Result[T] pattern for consistent error handling
4. **Maintainability**: Change business logic in one place
5. **Testability**: Test utilities once, API and CLI both benefit

---

## Recommended Action Plan

### Phase 1: High-Priority Refactoring (Immediate)

**Step 1**: Refactor `set_parameter/` endpoint
```python
# File: server/ccp4x/api/JobViewSet.py (line 1111)

# OLD:
from ..lib.job_utils.set_parameter import set_parameter

# NEW:
from ..lib.utils.parameters.set_param import set_parameter as set_job_param

@action(detail=True, methods=["post"])
def set_parameter(self, request, pk=None):
    form_data = json.loads(request.body.decode("utf-8"))
    job = models.Job.objects.get(id=pk)
    object_path = form_data["object_path"]
    value = form_data["value"]

    # Use unified utility (CPluginScript architecture)
    result = set_job_param(job, object_path, value)

    if result.success:
        return JsonResponse({"status": "Success", "updated_item": result.data})
    else:
        return JsonResponse({
            "status": "Failed",
            "reason": result.error,
            "details": result.error_details
        }, status=400)
```

**Step 2**: Refactor `validation/` endpoint
```python
# File: server/ccp4x/api/JobViewSet.py (line 1049)

# OLD:
from ..lib.job_utils.get_job_container import get_job_container
from ..lib.job_utils.validate_container import validate_container

# NEW:
from ..lib.utils.jobs.validate import validate_job

@action(detail=True, methods=["get"])
def validation(self, request, pk=None):
    try:
        the_job = models.Job.objects.get(id=pk)

        # Use unified utility (CPluginScript architecture)
        result = validate_job(the_job)

        if result.success:
            error_etree = result.data
            ET.indent(error_etree, " ")
            return Response({"status": "Success", "xml": ET.tostring(error_etree)})
        else:
            return Response({"status": "Failed", "reason": result.error}, status=400)

    except models.Job.DoesNotExist as err:
        logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
        return Response({"status": "Failed", "reason": str(err)}, status=404)
```

### Phase 2: Medium-Priority Consistency (Next)

**Step 3**: Create unified utility for `upload_file_param` (if needed)
**Step 4**: Create unified utility for `set_context_job` (if needed)

### Phase 3: Optional Improvements (Future)

**Step 5**: Refactor report endpoints to use unified utilities
**Step 6**: Consider refactoring advanced endpoints if usage justifies it

---

## Testing Strategy

For each refactored endpoint:

1. **Unit Test**: Test the unified utility function
2. **API Test**: Test the endpoint with real HTTP requests
3. **CLI Test**: Verify management command still works
4. **Integration Test**: Test complete workflow (create → set params → validate → run)

---

## Summary

### Current Status

- ✅ **9 endpoints** already use proper utilities (clone, run, export, etc.)
- ⚠️ **2 HIGH PRIORITY endpoints** need refactoring (set_parameter, validation)
- ⚠️ **2 MEDIUM PRIORITY endpoints** could be improved (upload_file_param, set_context_job)
- 🟢 **11 LOW PRIORITY endpoints** work but could be more consistent

### Immediate Next Steps

1. **Refactor `set_parameter/` endpoint** to use `ccp4x.lib.utils.parameters.set_param`
2. **Refactor `validation/` endpoint** to use `ccp4x.lib.utils.jobs.validate`
3. **Test both endpoints** with real HTTP requests
4. **Document changes** and update API documentation

### Long-term Goal

Achieve **100% consistency** where every endpoint either:
- Uses a unified utility function (shared with management commands), OR
- Has a clear architectural reason for direct implementation (e.g., simple ORM queries)

---

**Analysis Complete** ✅
**Priority Refactorings Identified**: 2 (set_parameter, validation)
**Ready to Implement**: Yes 🚀
