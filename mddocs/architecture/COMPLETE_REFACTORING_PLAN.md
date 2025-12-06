# Complete API Endpoint Refactoring Plan
## Goal: 13/24 → 24/24 Endpoints Using Unified Architecture

**Date**: 2025-10-31
**Current**: 13/24 endpoints use unified architecture
**Target**: 24/24 endpoints (100% consistency)

---

## Quick Status

### ✅ Already Done (13 endpoints)
1. `set_parameter/` - Uses `ccp4x.lib.utils.parameters.set_param` ✅
2. `validation/` - Uses `ccp4x.lib.utils.jobs.validate` ✅
3. `clone/` - Uses `clone_job()` utility ✅
4. `run/` - Uses `run_job_context_aware()` ✅
5. `run_local/` - Uses `run_job_context_aware()` ✅
6. `export_job/` - Uses `export_project_to_zip()` ✅
7. `export_job_file/` - Uses `export_job_file()` ✅
8. `export_job_file_menu/` - Uses TaskManager ✅
9. `files/` - Direct ORM (appropriate) ✅
10. `dependent_jobs/` - Uses `find_dependent_jobs()` ✅
11. `what_next/` - Uses `get_what_next()` ✅
12. `destroy()` (delete) - Uses `delete_job_and_dependents()` ✅
13. Standard REST (list, create, retrieve, update) - Django ORM ✅

### 🔴 To Refactor (11 endpoints)

**Phase 1: Medium Priority - Create Utilities & Refactor (5 endpoints)**
1. `upload_file_param/` - Create unified utility
2. `set_context_job/` - Create unified utility
3. `params_xml/` - Use existing `get_job_params_xml()`
4. `report_xml/` - Use existing `get_job_report_xml()`
5. `diagnostic_xml/` - Use existing `get_job_diagnostic_xml()`

**Phase 2: Low Priority - Refactor with New Utilities (6 endpoints)**
6. `object_method/` - Create utility if beneficial
7. `container/` - Create utility if beneficial
8. `digest/` - Create utility if beneficial
9. `digest_param_file/` - Create utility if beneficial
10. `i2run_command/` - Create utility if beneficial
11. `preview/` - Create utility if beneficial

---

## Refactoring Strategy

### Phase 1: Quick Wins (Report Endpoints)

These already have utilities in `ccp4x.lib.utils.jobs.reports`!

#### 1. `params_xml/` ✅ Utility Exists
**Current**: Direct filesystem read
**Target**: Use `get_job_params_xml(job)` from `ccp4x.lib.utils.jobs.reports`
**Effort**: 5 minutes
**Benefit**: Consistency + better error handling

#### 2. `report_xml/` ✅ Utility Exists
**Current**: Mix of `generate_job_report()` and filesystem
**Target**: Use `get_job_report_xml(job)` from `ccp4x.lib.utils.jobs.reports`
**Effort**: 5 minutes
**Benefit**: Consistency + caching logic centralized

#### 3. `diagnostic_xml/` ✅ Utility Exists
**Current**: Direct filesystem read
**Target**: Use `get_job_diagnostic_xml(job)` from `ccp4x.lib.utils.jobs.reports`
**Effort**: 5 minutes
**Benefit**: Consistency + better error handling

---

### Phase 2: Create New Utilities

#### 4. `upload_file_param/` ⚠️ Needs Utility
**Current**: Uses legacy `job_utils.upload_file_param`
**Target**: Create `ccp4x.lib.utils.parameters.upload_file.py`
**Effort**: 20 minutes (create utility + refactor endpoint)
**Benefit**: Consistency + shares code with potential CLI command

**New Utility Structure**:
```python
# ccp4x/lib/utils/parameters/upload_file.py
def upload_file_param(job: models.Job, request) -> Result[Dict]:
    """Upload file and set as parameter using CPluginScript."""
    # Extract file from request
    # Use get_plugin_with_context()
    # Set file parameter on plugin.container
    # Save and sync to DB
    return Result.ok(data)
```

#### 5. `set_context_job/` ⚠️ Needs Utility
**Current**: Uses legacy `job_utils.set_input_by_context_job`
**Target**: Create `ccp4x.lib.utils.jobs.context.py`
**Effort**: 20 minutes
**Benefit**: Consistency + workflow automation

**New Utility Structure**:
```python
# ccp4x/lib/utils/jobs/context.py
def set_input_by_context(job: models.Job, context_job_uuid: str) -> Result[models.Job]:
    """Set input parameters based on context job outputs."""
    # Get both jobs with context
    # Copy outputs from context job to inputs of target job
    # Use CPluginScript for proper hierarchy
    return Result.ok(job)
```

---

### Phase 3: Helper Endpoints (Optional)

These are less critical but we'll do them for completeness:

#### 6. `object_method/` 🟢 Low Priority
**Current**: Uses legacy `object_method()`
**Target**: Create `ccp4x.lib.utils.objects.method_call.py`
**Effort**: 15 minutes
**Benefit**: Minimal (specialized operation)

#### 7. `container/` 🟢 Low Priority
**Current**: Uses `json_for_job_container()`
**Target**: Create `ccp4x.lib.utils.jobs.container_json.py`
**Effort**: 15 minutes
**Benefit**: Minimal (specialized operation)

#### 8-11. Digest & Helper Endpoints 🟢 Low Priority
**Effort**: 10 minutes each
**Benefit**: Completeness only

---

## Execution Plan

### Step 1: Quick Wins (15 minutes) ✅
- Refactor `params_xml/` to use `get_job_params_xml()`
- Refactor `report_xml/` to use `get_job_report_xml()`
- Refactor `diagnostic_xml/` to use `get_job_diagnostic_xml()`
- **Result**: 16/24 endpoints done

### Step 2: Create Upload File Utility (20 minutes) ✅
- Create `ccp4x.lib.utils.parameters.upload_file.py`
- Refactor `upload_file_param/` endpoint
- **Result**: 17/24 endpoints done

### Step 3: Create Context Job Utility (20 minutes) ✅
- Create `ccp4x.lib.utils.jobs.context.py`
- Refactor `set_context_job/` endpoint
- **Result**: 18/24 endpoints done

### Step 4: Helper Utilities (60 minutes) ✅
- Create utilities for remaining 6 endpoints
- Refactor each endpoint
- **Result**: 24/24 endpoints done! 🎉

---

## Total Estimated Time

- **Phase 1 (Quick Wins)**: 15 minutes
- **Phase 2 (New Utilities)**: 40 minutes
- **Phase 3 (Helpers)**: 60 minutes
- **Testing**: 20 minutes
- **Total**: ~2-2.5 hours

---

## Testing Strategy

After each phase:
1. Run management command tests (if they exist)
2. Quick smoke test of endpoint (curl)
3. Verify no regressions

Final validation:
- Run end-to-end workflow test
- Test a few endpoints via HTTP
- Verify all management commands work

---

## Success Criteria

✅ All 24 endpoints use unified utilities
✅ Consistent Result[T] pattern everywhere
✅ Proper error handling (HTTP status codes)
✅ No regressions (existing tests still pass)
✅ Documentation updated

---

## Let's Go! 🚀

**Starting with**: Phase 1 - Quick Wins (Report Endpoints)
**Target**: 16/24 endpoints in 15 minutes
**Current**: 13/24 endpoints
