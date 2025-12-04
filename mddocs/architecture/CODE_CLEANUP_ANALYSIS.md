# Code Cleanup Analysis: job_utils Refactoring

**Date**: 2025-10-31
**Purpose**: Identify redundant code in `server/ccp4x/lib/job_utils/` that can be safely removed after API endpoint refactoring
**Current Status**: 16/24 API endpoints refactored to use unified architecture

---

## Executive Summary

After refactoring 16 API endpoints to use the unified CPluginScript architecture, we can now identify and remove redundant legacy code. This analysis categorizes all 44 modules in `job_utils/` into:

- **🔴 CRITICAL - Keep**: Core infrastructure (3 modules)
- **🟡 REDIRECT**: Duplicates with replacements in `lib/utils/` (3 modules)
- **🟢 ACTIVE USE**: Still used by API endpoints, no replacement yet (15 modules)
- **🔵 INTERNAL HELPERS**: Can migrate to `lib/utils/` (20 modules)
- **⚫ SAFE TO REMOVE**: Completely unused (3 modules)

**Quick Actions Available**:
- Remove 3 completely unused modules immediately
- Redirect 3 API imports to use `lib/utils/` versions
- Deprecate 4 top-level functions that have been replaced

---

## Category Breakdown

### 🔴 CRITICAL - Keep (Core Infrastructure)

These provide essential infrastructure used across the entire system. **DO NOT REMOVE**.

#### 1. `get_job_plugin.py` ✅ KEEP
**Status**: Core infrastructure - NOT replaced
**Used by**: 22+ files across the codebase

**Key Usage**:
- `lib/utils/plugins/plugin_context.py` - Canonical entry point
- `lib/utils/parameters/set_param.py`
- `lib/utils/files/upload_param.py`
- `async_run_job.py` (2 calls)
- Multiple test files

**Why**: Instantiates plugin with dbHandler attachment and proper context

**Action**: ✅ **KEEP** - Essential for plugin system

---

#### 2. `validate_container.py::getEtree()` ✅ KEEP
**Status**: Essential XML conversion utility
**Used by**: Error handling across API

**Specific API Endpoints**:
- `JobViewSet.set_context_job()` (line 287)
- `JobViewSet.object_method()` (line 354)
- Import at line 47 in JobViewSet

**Why**: Converts `CErrorReport` to XML ElementTree for API responses

**Note**: The `validate_container()` function itself is replaced, but `getEtree()` helper is still needed

**Action**: ✅ **KEEP** - But consider moving to a more appropriate location like `lib/utils/errors/`

---

#### 3. `set_parameter.py::set_parameter_container()` ✅ KEEP
**Status**: Internal helper - still used by new unified utilities
**Used by**:
- `lib/utils/parameters/set_param.py` (wrapper)
- `lib/utils/plugins/plugin_context.py` (line 50)
- `i2run/CCP4i2RunnerBase.py`

**Why**: Low-level container manipulation logic

**Note**: The top-level `set_parameter()` function is replaced, but the internal `set_parameter_container()` helper is still essential

**Action**: ✅ **KEEP** - Internal helper needed by refactored code

---

### 🟡 REDIRECT - Duplicates with Replacements (3 modules)

These have newer versions in `lib/utils/`. Redirect imports and deprecate originals.

#### 1. `upload_file_param.py` 🟡 REDIRECT
**Status**: Duplicate - identical code exists in `lib/utils/files/upload_param.py`
**Current API usage**: `JobViewSet.upload_file_param()` (line 67)

**Replacement**: `ccp4x.lib.utils.files.upload_param.upload_file_param()`

**Action Plan**:
1. Update JobViewSet line 67 to import from `lib.utils.files.upload_param`
2. Add deprecation warning to `job_utils/upload_file_param.py`
3. Mark for removal in future version

**Effort**: 2 minutes

---

#### 2. `set_input_by_context_job.py` 🟡 REDIRECT
**Status**: Replacement exists in `lib/utils/parameters/set_input_by_context.py`
**Current API usage**: `JobViewSet.set_context_job()` (line 48)

**Replacement**: `ccp4x.lib.utils.parameters.set_input_by_context.set_input_by_context_job()`

**Action Plan**:
1. Update JobViewSet line 48 to import from `lib.utils.parameters.set_input_by_context`
2. Add deprecation warning to `job_utils/set_input_by_context_job.py`
3. Mark for removal in future version

**Effort**: 2 minutes

---

#### 3. `get_job_container.py` 🟡 REDIRECT
**Status**: Wrapper exists in `lib/utils/jobs/get_container.py`
**Current usage**:
- `JobViewSet.container()` endpoint (line 68)
- Test files

**Replacement**: `ccp4x.lib.utils.jobs.get_container.get_job_container()`

**Action Plan**:
1. Update JobViewSet line 68 to import from `lib.utils.jobs.get_container`
2. Update test files to use new import
3. Add deprecation warning
4. Mark for removal in future version

**Effort**: 5 minutes (need to update tests too)

---

### 🟢 ACTIVE USE - No Replacement Yet (15 modules)

These are actively used by API endpoints and have no replacement. **Keep for now**.

| Module | Used By | Purpose |
|--------|---------|---------|
| `clone_job.py` | JobViewSet.clone/ | Job cloning |
| `run_job.py` | JobViewSet.run/ | Job execution |
| `context_dependent_run.py` | JobViewSet.run_local/ | Local execution |
| `i2run_for_job.py` | JobViewSet.i2run_command/ | CLI command generation |
| `digest_file.py` | JobViewSet.digest_param_file/ | File digest computation |
| `preview_job.py` | JobViewSet.preview/ | Launch external viewers |
| `preview_file.py` | FileViewSet, ProjectViewSet | File preview |
| `ccp4i2_report.py` | JobViewSet.report_xml/ | Report generation |
| `get_what_next.py` | JobViewSet.what_next/ | Workflow recommendations |
| `find_dependent_jobs.py` | JobViewSet.destroy/ | Dependency tracking |
| `object_method.py` | JobViewSet.object_method/ | Dynamic method calls |
| `json_for_job_container.py` | JobViewSet.container/ | JSON serialization |
| `export_job_file.py` | JobViewSet.export_job_file/ | File export |
| `load_nested_xml.py` | JobViewSet (multiple) | Config loading |
| `list_project.py` | Management command | Project listing |

**Action**: ✅ **KEEP** - All actively used, no replacements available

---

### 🔵 INTERNAL HELPERS - Can Migrate (20 modules)

These are only used internally by other `job_utils` modules. Can be migrated to `lib/utils/` over time.

Many already have equivalents in `lib/utils/`:

| job_utils Module | Equivalent in lib/utils/ | Status |
|------------------|-------------------------|--------|
| `find_objects.py` | `lib/utils/containers/find_objects.py` | ✅ Exists |
| `value_dict_for_object.py` | `lib/utils/parameters/value_dict.py` | ✅ Exists |
| `available_file_name_based_on.py` | `lib/utils/files/available_name.py` | ✅ Exists |
| `detect_file_type.py` | `lib/utils/files/detect_type.py` | ✅ Exists |
| `json_encoder.py` | `lib/utils/containers/json_encoder.py` | ✅ Exists |
| `save_params_for_job.py` | `lib/utils/parameters/save_params.py` | ✅ Moved |
| `gemmi_split_mtz.py` | `lib/utils/formats/gemmi_split_mtz.py` | ✅ Moved |

**Other internal helpers** (no immediate action needed):
- `create_job.py`
- `create_task.py`
- `export_job_mtz_file.py`
- `get_file_by_job_context.py`
- `get_source_reflection_file.py`
- `get_task_tree.py`
- `glean_job_files.py`
- `job_directory.py`
- `mtz_as_dict.py`
- `patch_output_file_paths.py`
- `remove_container_default_values.py`
- `set_output_file_names.py`
- `unset_output_data.py`

**Action**: 🔄 **Gradual migration** - Move to `lib/utils/` as they're encountered during refactoring

---

### ⚫ SAFE TO REMOVE - Completely Unused (3 modules)

These have **ZERO** usages across the entire codebase. Safe to delete immediately.

#### 1. `open_terminal_in_directory.py` ⚫ DELETE
**Status**: 0 usages found
**Verified**: No imports, no references

**Action**: ✅ **SAFE TO REMOVE** immediately

---

#### 2. `parse_cif_ligand_summary.py` ⚫ DELETE
**Status**: 0 usages found
**Verified**: No imports, no references

**Action**: ✅ **SAFE TO REMOVE** immediately

---

#### 3. `import_files.py` ⚫ DELETE
**Status**: 0 usages found
**Verified**: No imports, no references

**Action**: ✅ **SAFE TO REMOVE** immediately

---

## Immediate Action Plan

### Phase 1: Quick Wins (10 minutes) ✅

**Remove completely unused modules**:
```bash
rm server/ccp4x/lib/job_utils/open_terminal_in_directory.py
rm server/ccp4x/lib/job_utils/parse_cif_ligand_summary.py
rm server/ccp4x/lib/job_utils/import_files.py
```

**Impact**: Clean up 3 unused files, reduce confusion

---

### Phase 2: Redirect API Imports (10 minutes) ✅

**Update `JobViewSet.py` imports**:

1. Line 67 - `upload_file_param`:
```python
# OLD
from ..lib.job_utils.upload_file_param import upload_file_param

# NEW
from ..lib.utils.files.upload_param import upload_file_param
```

2. Line 48 - `set_input_by_context_job`:
```python
# OLD
from ..lib.job_utils.set_input_by_context_job import set_input_by_context_job

# NEW
from ..lib.utils.parameters.set_input_by_context import set_input_by_context_job
```

3. Line 68 - `get_job_container`:
```python
# OLD
from ..lib.job_utils.get_job_container import get_job_container

# NEW
from ..lib.utils.jobs.get_container import get_job_container
```

**Impact**: 3 API endpoints now use unified utilities

---

### Phase 3: Add Deprecation Warnings (15 minutes)

Add deprecation warnings to old modules:

```python
# At top of server/ccp4x/lib/job_utils/upload_file_param.py
import warnings
warnings.warn(
    "ccp4x.lib.job_utils.upload_file_param is deprecated. "
    "Use ccp4x.lib.utils.files.upload_param instead.",
    DeprecationWarning,
    stacklevel=2
)
```

Repeat for:
- `set_input_by_context_job.py`
- `get_job_container.py`
- Top-level `set_parameter()` in `set_parameter.py`
- `validate_container()` in `validate_container.py`

**Impact**: Clear migration path for any remaining legacy code

---

### Phase 4: Update Tests (20 minutes)

Update test files to import from `lib/utils/` instead of `job_utils/`:

Files to update:
- `tests/api/test_job_utils.py`
- `tests/lib/test_container_utils.py`
- `tests/lib/test_export_mtz.py`
- `tests/lib/test_report.py`

**Impact**: Tests use modern import paths

---

## Summary Statistics

### Current State
- **Total job_utils modules**: 44
- **Critical (keep)**: 3 (7%)
- **Redirect available**: 3 (7%)
- **Active use (keep)**: 15 (34%)
- **Internal helpers**: 20 (45%)
- **Completely unused**: 3 (7%)

### After Cleanup
- **Removed**: 3 modules (7%)
- **Redirected**: 3 API imports (7%)
- **Deprecated**: 5 functions (11%)
- **Remaining job_utils**: 41 modules (93%)

---

## Dependency Graph

### Top-Level Entry Points (Keep)
```
get_job_plugin.py ← Used by 22+ files
    ↓
validate_container.py::getEtree() ← Error handling
    ↓
set_parameter.py::set_parameter_container() ← Low-level operations
```

### Refactored APIs Now Use lib/utils/
```
OLD: job_utils/set_parameter.py::set_parameter()
NEW: lib/utils/parameters/set_param.py::set_parameter()
    ↳ Still uses: set_parameter_container() [KEEP]

OLD: job_utils/validate_container.py::validate_container()
NEW: lib/utils/jobs/validate.py::validate_job()
    ↳ Still uses: getEtree() [KEEP]

OLD: job_utils/upload_file_param.py
NEW: lib/utils/files/upload_param.py [DUPLICATE - redirect API]

OLD: job_utils/set_input_by_context_job.py
NEW: lib/utils/parameters/set_input_by_context.py [DUPLICATE - redirect API]

OLD: job_utils/get_job_container.py
NEW: lib/utils/jobs/get_container.py [WRAPPER - redirect API]
```

---

## Risk Assessment

### Zero Risk ✅
- Removing 3 completely unused modules
- Adding deprecation warnings

### Very Low Risk ✅
- Redirecting 3 API imports to `lib/utils/` (code is identical)
- Updating test imports

### Low Risk 🟡
- Gradually migrating internal helpers to `lib/utils/`

### Medium Risk 🟠
- Eventually removing deprecated modules after grace period

---

## Recommendations

### Immediate (Today) ✅
1. ✅ Remove 3 completely unused modules
2. ✅ Redirect 3 API imports to `lib/utils/`
3. ✅ Add deprecation warnings

**Total Time**: 35 minutes
**Impact**: Cleaner codebase, clearer migration path

### Short-term (This Week)
4. Update test imports to use `lib/utils/`
5. Document migration guide for any external code

**Total Time**: 20 minutes
**Impact**: Tests use modern patterns

### Long-term (Next Sprint)
6. Create script to detect remaining `job_utils` imports
7. Migrate internal helpers to `lib/utils/` as encountered
8. After 1-2 releases, remove deprecated modules

**Total Time**: 2-3 hours over time
**Impact**: Complete migration to unified architecture

---

## Migration Checklist

### For Each Deprecated Module

- [ ] Identify all current usages (done via this analysis)
- [ ] Confirm replacement exists in `lib/utils/`
- [ ] Add deprecation warning to old module
- [ ] Update API endpoint imports
- [ ] Update test imports
- [ ] Document change in migration guide
- [ ] Wait 1-2 releases
- [ ] Remove deprecated module

---

## Conclusion

**We can safely clean up 7% of job_utils immediately** (3 unused modules) and **redirect another 7%** (3 API imports) with minimal effort.

The remaining modules fall into two categories:
1. **Core infrastructure** (7%) - Keep forever
2. **Active use** (86%) - Migrate gradually as we refactor more endpoints

**Key Insight**: The refactoring work has been successful - we've created a clean unified architecture in `lib/utils/`, and now we can systematically deprecate and remove the legacy code.

**Next Step**: Execute Phase 1-3 of the action plan (35 minutes total) to clean up the codebase immediately.

---

**Analysis Complete** ✅
**Safe Removal Candidates Identified** ✅
**Migration Path Documented** ✅
**Ready for Cleanup** ✅
