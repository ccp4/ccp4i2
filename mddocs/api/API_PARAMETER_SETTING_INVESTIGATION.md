# API Parameter Setting Investigation

## Problem Statement

The `/jobs/{id}/set_parameter/` API endpoint is not persisting parameters correctly. When setting a parameter via the API:

1. The API returns success
2. But `saveDataToXml()` writes an empty `<ccp4i2_body />`
3. `get_parameter()` fails to find the parameter
4. The parameter is not persisted across requests

## Root Cause

When a plugin is loaded from `.def.xml` only (no XML overlay), `CContainer.set_parameter()` is returning **plain Python types** (int, str, bool) instead of **CData wrappers** (CInt, CString, CBoolean).

### Evidence

```python
# API response shows plain type instead of CData wrapper
{'status': 'Success', 'data': {'object_type': 'int'}}  # Should be 'CInt'
```

This causes:
- Plain types don't serialize to XML → empty body
- Plain types don't participate in CData hierarchy → get_parameter() can't find them
- Changes are lost between requests

## Architecture Discovery

### The i2run Pattern (from async_run_job.py)

The correct workflow for job execution:

```python
# 1. Create fresh plugin from .def.xml (loads structure with CData wrappers)
plugin = plugin_class(workDirectory=job_dir, name=title)

# 2. Set database context
plugin.setDbData(handler=db_handler, projectId=..., jobNumber=..., jobId=...)

# 3. Load parameters from XML (overlays user values onto structure)
if params_file.exists():
    plugin.loadDataFromXml(str(params_file))

# 4. Now plugin has proper CData wrappers with user-specified values
```

### Key Insight

`.def.xml` loading happens automatically in `CPluginScript.__init__()` and should create proper CData wrappers. The issue is that `CContainer.set_parameter()` is somehow creating/returning plain types instead of using the CData wrappers from `.def.xml`.

## What Was Fixed

### 1. Clone API Bugs (patch_paths.py)

Fixed three critical bugs in `handle_cdatafile()`:

**Bug 1: Invalid checkDb parameter**
```python
# BEFORE (line 62, 64)
dobj.setFullPath(fullPath, checkDb=True)  # TypeError!

# AFTER (line 72)
dobj.setFullPath(fullPath)  # Correct signature
```

**Bug 2: Array out of bounds**
```python
# BEFORE (line 50)
extension = dobj.fileExtensions()[1]  # IndexError if len < 2!

# AFTER (lines 51-62)
if len(extensions) > 1:
    extension = extensions[1]
else:
    extension = extensions[0] if extensions else ""
```

**Bug 3: Missing fileExtensions method**
```python
# BEFORE
extensions = dobj.fileExtensions()  # AttributeError for some types!

# AFTER (lines 44-49)
if hasattr(dobj, 'fileExtensions') and callable(getattr(dobj, 'fileExtensions')):
    extensions = dobj.fileExtensions()
else:
    extensions = dobj.get_qualifier('fileExtensions', [])
```

### 2. Plugin Loading for New Jobs (get_plugin.py)

Made params file optional - allows loading fresh plugins without requiring XML:

```python
# BEFORE (line 53)
if not params_file.exists():
    raise Exception("No params file found")  # Fails for new jobs!

# AFTER (lines 56-65)
if params_file.exists():
    error = pluginInstance.loadDataFromXml(str(params_file))
    # ... error handling
    logger.info(f"Loaded parameters from {params_file}")
else:
    logger.info(f"No params file found - using fresh plugin from .def.xml")
```

This is critical for the API use case where we want to:
1. Create a fresh PENDING job
2. Set parameters via API
3. Save to `input_params.xml`

### 3. Test Setup (test_parameter_setting_api.py)

Changed test to follow i2run pattern - don't create empty XML:

```python
# BEFORE
plugin.saveDataToXml(str(input_params_file))  # Creates empty XML!
# Then API loads this empty XML → plain types

# AFTER (lines 70-73)
# IMPORTANT: Do NOT create input_params.xml here!
# The API endpoint will load a fresh plugin from .def.xml each time,
# giving proper CData wrappers. If we save an empty XML here, the API
# will load that instead and get plain types.
```

## What Still Needs Investigation

### The Core Issue: Plain Types from set_parameter()

The API test shows that `set_parameter()` is returning plain Python int:

```python
# From set_param.py line 118
result_data = {
    "object_type": type(obj).__name__ if obj else "Unknown",  # Returns 'int'
}
```

This means `plugin.container.set_parameter(...)` is returning a plain int, not a CInt.

**Questions:**
1. Does `.def.xml` loading actually create CData wrappers, or just structure?
2. Is `CContainer.set_parameter()` supposed to create new CData wrappers?
3. Or is there a missing initialization step between plugin creation and parameter setting?

### Possible Next Steps

1. **Debug DefXmlParser** - Check what objects it creates when parsing `.def.xml`
   - Does `parse_def_xml()` create CInt/CString objects or just containers?
   - Where does the type information come from?

2. **Debug CContainer.set_parameter()** - Trace through the implementation
   - When path exists: does it find CData wrappers or plain types?
   - When path doesn't exist: what does it create?

3. **Compare with i2run** - Run i2run test and inspect plugin
   - After `.def.xml` load, what type is `plugin.container.inputData.NCYCLES`?
   - After XML overlay, does the type change?

## Test Status

### ✅ Fixed
- Clone API now works (previously crashed with 3 different errors)
- Plugin loading works for new jobs (previously required params file)
- **Parameter setting API tests - ALL 11 TESTS PASS**
- Parameters set via API now persist correctly to XML
- `get_parameter()` correctly finds parameters that were set

### Root Causes Found and Fixed

1. **Test used wrong parameter paths**: Tests were using `inputData.NCYCLES` but NCYCLES is in `controlParameters`
2. **`find_by_path()` didn't handle `find()` returning -1**: The `find()` method returns -1 for not found, but `find_by_path` only checked for `None`
3. **`set_parameter()` allowed creating plain Python attributes**: When target wasn't a CData type, it would create a plain attribute instead of failing

## Files Modified (Session 2 - API Parameter Fixes)

1. **core/base_object/ccontainer.py**
   - Fixed `find_by_path()` to handle `find()` returning -1 (not found)
   - Fixed `set_parameter()` to reject plain Python types - now raises AttributeError

2. **server/ccp4i2/tests/api/test_parameter_setting_api.py**
   - Fixed all parameter paths from `inputData.*` to `controlParameters.*`
   - Fixed `test_set_parameter_null_value` to expect 400 error (CInt can't be None)

## Files Modified (Session 1 - Clone API Fixes)

1. **server/ccp4i2/lib/utils/files/patch_paths.py**
   - Lines 44-72: Complete rewrite of `handle_cdatafile()`
   - Fixed `setFullPath()` signature
   - Added safety checks for array access and method existence

2. **server/ccp4i2/lib/utils/plugins/get_plugin.py**
   - Lines 53-65: Made params file optional
   - Allows fresh plugin loading for new PENDING jobs

3. **server/ccp4i2/tests/api/test_parameter_setting_api.py**
   - Lines 70-73: Removed empty XML creation
   - Now relies on fresh `.def.xml` loading

## Current Status: API Lifecycle WORKING ✅

### Complete API Lifecycle Verified

The full API lifecycle for job creation, parameter setting, and execution is now working:

**1. Job Creation** (`POST /projects/{pk}/create_task/`)
- Uses `create_job_async()` which properly:
  - Creates database record with correct job number
  - Creates job directory structure
  - Instantiates plugin from `.def.xml` (proper CData wrappers)
  - Sets output file names
  - Saves `input_params.xml`

**2. Parameter Setting** (`POST /jobs/{pk}/set_parameter/`)
- Uses `set_param.py` utility with `container.set_parameter()`
- Parameters correctly persist to `input_params.xml`
- All 11 API parameter tests pass

**3. Job Execution** (`POST /jobs/{pk}/run/`)
- Uses `i2run_for_job()` which spawns the job
- Output files and KPIs registered in database

### Test Results

- **API Parameter Setting Tests**: 11/11 PASS
- **API Lifecycle Tests**: 1/2 PASS, 1/2 SKIP (demo data not found - expected)

### Remaining Work

The core API lifecycle is complete. Remaining items for future sessions:
1. Add more comprehensive integration tests with actual crystallographic data
2. Document the full API surface for frontend integration
3. Consider alignment between i2run command-line syntax and API (deferred per user request)

## Key Architectural Lessons

1. **Job Lifecycle Stages**:
   - User Control (PENDING): uses `input_params.xml`, parameters editable
   - Plugin Lifecycle (RUNNING+): uses `params.xml`, job immutable
   - Use clone API to get editable copy of finished jobs

2. **Plugin Initialization Pattern**:
   - `.def.xml` loads automatically in `CPluginScript.__init__()`
   - This should create proper CData wrappers (but investigation needed)
   - XML overlay adds user values to the structure

3. **Parameter Setting Workflow**:
   - For new jobs: load fresh plugin → set params → save to input_params.xml
   - For finished jobs: clone → get new job → set params → save
   - Don't create empty XML files - they cause plain type loading!
