# ViewSet Modernization Progress

## Completed Modernizations

### ✅ 1. JobViewSet `container` Endpoint (HIGH PRIORITY)

**Status**: COMPLETE ✅
**Test**: `test_job_container` - **PASSING**

**Changes Made**:
1. Updated imports in `JobViewSet.py`:
   ```python
   from ..lib.utils.plugins.get_plugin import get_job_plugin
   from ..lib.utils.containers.json_encoder import CCP4i2JsonEncoder
   ```

2. Modernized endpoint implementation:
   - Uses `get_job_plugin()` instead of legacy `json_for_job_container()`
   - Direct serialization with `CCP4i2JsonEncoder`
   - Better error handling with proper HTTP status codes
   - Follows pattern from management commands

3. Fixed `get_plugin.py` import path:
   - Changed from `...db.models` to `ccp4i2.db.models`

4. Fixed `get_plugin.py` API compatibility:
   - Added fallback for `loadDataFromXml()` signature variations

5. Fixed `json_encoder.py` imports:
   - Updated from `core.CCP4Data.CData` to `core.base_object.CData`
   - Added proper imports for fundamental types

6. Added `core/base_object/__init__.py`:
   - Exported `CData`, `CContainer`, etc. for easier importing

### ✅ 2. CData CONTENTS Property (INFRASTRUCTURE)

**Status**: COMPLETE ✅

**Added to** `core/base_object/cdata.py`:

```python
@property
def CONTENTS(self):
    """
    Get list of child CData objects or attributes.

    For CContainer: Returns list of child CData objects (children)
    For other CData types: Returns list of attribute names that are CData objects
    """
    from .ccontainer import CContainer

    # For CContainer, return actual children objects
    if isinstance(self, CContainer):
        return list(self.get_children())

    # For other CData types, return list of CData attribute names
    # ... (metadata-aware + fallback implementation)
    return cdata_attributes
```

**Benefits**:
- Unified interface for navigating CData hierarchy
- Compatible with legacy CCP4i2 CONTENTS pattern
- Metadata-aware (uses MetadataRegistry when available)
- Works for both CContainer (children) and other CData (attributes)

### ✅ 3. Dependencies & Configuration

**Added**:
- `scipy>=1.10.0` to requirements.txt
- `numpy>=1.24.0` to requirements.txt
- `from rest_framework import status` to JobViewSet

**Fixed**:
- Test settings now include REST Framework + auth apps
- URL routing configured for tests

## Test Results

### JobViewSet Tests (23 total)
- **15 PASSED** ✅ (65%)
- **5 FAILED** (parameter setting issues - not yet modernized)
- **3 SKIPPED** (intentional - external tools)

**Passing Tests**:
- ✅ test_job_clone
- ✅ test_job_container (NEWLY FIXED!)
- ✅ test_job_delete
- ✅ test_job_dependent_jobs
- ✅ test_job_diagnostic_xml
- ✅ test_job_digest
- ✅ test_job_digest_param_file
- ✅ test_job_files
- ✅ test_job_i2run_command
- ✅ test_job_list
- ✅ test_job_params_xml
- ✅ test_job_report_xml
- ✅ test_job_retrieve
- ✅ test_job_validation
- ✅ test_job_what_next

**Failing Tests** (need modernization):
- ❌ test_job_def_xml
- ❌ test_job_export_job_file_menu
- ❌ test_job_set_parameter_file
- ❌ test_job_set_parameter_simple
- ❌ test_job_upload_file_param

### FileViewSet Tests
- **7 PASSED** ✅
- **1 SKIPPED** (intentional)

### ProjectViewSet Tests
- Not yet fully tested

## Architecture Improvements

### Before (Legacy)
```python
# Old approach
from ..lib.job_utils.json_for_job_container import json_for_job_container

def container(self, request, pk=None):
    job = models.Job.objects.get(id=pk)
    result = json_for_job_container(job)  # Returns dict
    return JsonResponse({"status": "Success", "result": result})
```

**Problems**:
- No Result pattern
- Inconsistent with management commands
- No proper HTTP status codes
- Limited error handling

### After (Modern)
```python
# Modern approach
from ..lib.utils.plugins.get_plugin import get_job_plugin
from ..lib.utils.containers.json_encoder import CCP4i2JsonEncoder

def container(self, request, pk=None):
    job = models.Job.objects.get(id=pk)

    # Use CPluginScript architecture
    plugin = get_job_plugin(job)

    # Serialize with modern encoder
    container_json = json.dumps(plugin.container, cls=CCP4i2JsonEncoder)

    return Response({
        "status": "Success",
        "result": json.loads(container_json)
    })
```

**Benefits**:
- ✅ Uses CPluginScript architecture
- ✅ Consistent with management commands
- ✅ Proper HTTP status codes (404, 500)
- ✅ Better error handling
- ✅ Type-safe with proper imports
- ✅ Database context handled correctly

## Remaining Work

### HIGH PRIORITY
- [ ] `digest` endpoints (FileViewSet & JobViewSet)
- [ ] `create_task` endpoint (ProjectViewSet)
- [ ] `what_next` endpoint (JobViewSet) - **Already works, may just need documentation**

### MEDIUM PRIORITY
- [ ] `upload_file_param` endpoint
- [ ] `object_method` endpoint
- [ ] `set_context_job` endpoint
- [ ] `def_xml` endpoint
- [ ] `directory` endpoint (ProjectViewSet)

### Test Failures to Address
The set_parameter failures are due to tests using non-existent parameters. Need to either:
1. Fix tests to use actual plugin parameters
2. Or modernize the endpoints to handle missing parameters better

## Key Learnings

1. **Import Paths Matter**: The relative imports (`...db`) don't work from all locations - use absolute (`ccp4i2.db`)

2. **API Compatibility**: `loadDataFromXml()` signature varies - need try/except fallback

3. **Module Structure**: `CCP4Data.CData` doesn't exist - import from `core.base_object`

4. **Export Important Classes**: Adding `__init__.py` exports makes importing much easier

5. **CONTENTS Property**: Provides clean, unified interface for CData hierarchy navigation

## Documentation Updates Needed

- [x] Created VIEWSET_MODERNIZATION_MAP.md
- [x] Created MODERNIZATION_EXAMPLE.md
- [x] Created API_ENDPOINT_REFERENCE.md
- [x] Created this MODERNIZATION_PROGRESS.md
- [ ] Update main CLAUDE.md with CONTENTS property info
- [ ] Document CPluginScript patterns

## Next Session

Start with:
1. Modernize `digest` endpoints (files are easier, already have modern utilities)
2. Then tackle `what_next` (already passing, might just need verification)
3. Then `create_task` (important for workflow)

## Commands for Testing

```bash
# Single test
export CCP4I2_ROOT=$CCP4I2_ROOT
pytest server/ccp4i2/tests/api/test_viewsets_comprehensive.py::JobViewSetTests::test_job_container -xvs

# Full JobViewSet suite
pytest server/ccp4i2/tests/api/test_viewsets_comprehensive.py::JobViewSet Tests -v --tb=no

# All API tests
pytest server/ccp4i2/tests/api/ -v
```

## Success Metrics

- **Before**: 14 JobViewSet tests passing
- **After**: 15 JobViewSet tests passing (+1) ✅
- **Target**: 18 passing (all except parameter setting edge cases)

The modernization is working! Each endpoint updated brings us closer to a fully modern, maintainable API built on CPluginScript architecture.
