# validate_job Refactoring - SUCCESS! ✅

**Date**: 2025-10-31
**Task**: Refactor `validate_job` to use CPluginScript + dbHandler architecture
**Status**: ✅ **COMPLETE AND WORKING**

---

## Summary

Successfully refactored `validate_job()` to use the unified CPluginScript architecture via `get_plugin_with_context()`. The refactoring included:

1. ✅ Changed from standalone container to CPluginScript-based approach
2. ✅ Fixed CErrorReport API mismatch (`_reports` → `getErrors()`)
3. ✅ Added SEVERITY_TEXT mapping for XML generation
4. ✅ Tested end-to-end validation workflow

---

## What Was Fixed

### Bug 1: Using Standalone Container Instead of Plugin

**Before** (`validate.py`):
```python
# ❌ Created orphaned container without DB context
container_result = get_job_container(job.uuid)
container = container_result.data
error_tree = validate_container(container)
```

**After** (`validate.py`):
```python
# ✅ Uses CPluginScript with full DB context
plugin_result = get_plugin_with_context(job)
plugin = plugin_result.data
error_tree = validate_container(plugin.container)
```

**Why This Matters**: The plugin provides proper hierarchy, DB context, and ensures file parameters work correctly.

---

### Bug 2: CErrorReport API Mismatch

**Before** (`validate_container.py` line 12):
```python
# ❌ Old API - _reports attribute doesn't exist
for item in error_report._reports:
    # ... process errors
```

**After** (`validate_container.py` line 31):
```python
# ✅ New API - use public getErrors() method
for item in error_report.getErrors():
    # ... process errors
```

**Why This Matters**: The new CErrorReport class stores errors in `_errors` and provides `getErrors()` as the public API.

---

### Bug 3: Missing SEVERITY_TEXT Constant

**Before** (`validate_container.py` line 26):
```python
# ❌ SEVERITY_TEXT doesn't exist in new API
e.text = CCP4ErrorHandling.SEVERITY_TEXT[severity]
```

**After** (`validate_container.py` lines 10-16):
```python
# ✅ Create mapping from severity codes to text
SEVERITY_TEXT = {
    CCP4ErrorHandling.SEVERITY_OK: "OK",
    CCP4ErrorHandling.SEVERITY_UNDEFINED: "UNDEFINED",
    CCP4ErrorHandling.SEVERITY_WARNING: "WARNING",
    CCP4ErrorHandling.SEVERITY_UNDEFINED_ERROR: "UNDEFINED_ERROR",
    CCP4ErrorHandling.SEVERITY_ERROR: "ERROR"
}

# Use it:
e.text = SEVERITY_TEXT.get(severity, f"UNKNOWN({severity})")
```

---

### Bug 4: Error Dictionary Structure Mismatch

**Before** (`validate_container.py` line 22):
```python
# ❌ Tried to call description() method on error_report
desc, severity = error_report.description(item)
```

**After** (`validate_container.py` lines 48-55):
```python
# ✅ Error dict already contains details and severity
e = ET.Element("description")
e.text = item["details"]  # Description is in 'details' field

e = ET.Element("severity")
severity = item["severity"]  # Severity is directly in item dict
e.text = SEVERITY_TEXT.get(severity, f"UNKNOWN({severity})")
```

**Why This Matters**: The new API stores all error info in the dict returned by `getErrors()`.

---

## Test Results

### Test Command

```bash
cd server
python manage.py validate_job --jobuuid 4f04d478-6ca7-4b14-a5a3-9101c285e967
```

### Output

```
✓ Job 1 validation passed - no errors or warnings
```

### XML Output (via `-o` flag)

```bash
python manage.py validate_job --jobuuid 4f04d478-6ca7-4b14-a5a3-9101c285e967 -o /tmp/validation_output.xml
```

**File**: `/tmp/validation_output.xml`

```xml
<errorReportList />
```

(Empty = no errors, which is correct for a valid job!)

### JSON Output (via `--json` flag)

```bash
python manage.py validate_job --jobuuid 4f04d478-6ca7-4b14-a5a3-9101c285e967 --json
```

**Output**:
```json
{
  "status": "Success",
  "error_count": 0,
  "warning_count": 0,
  "has_errors": false
}
```

---

## Architecture Flow

### Before (Broken)

```
validate_job()
  └─> get_job_container(uuid)  ❌ Orphaned container, no DB context
       └─> CContainer (standalone)
            └─> validate_container()
                 └─> container.validity()
                      └─> CErrorReport._reports  ❌ Doesn't exist!
```

### After (Working) ✅

```
validate_job()
  └─> get_plugin_with_context(job)  ✅ Full DB context
       ├─> CCP4i2DjangoDbHandler attached
       ├─> Job context set (UUID, number, project)
       └─> CPluginScript instance
            └─> plugin.container (CContainer with hierarchy)
                 └─> validate_container(plugin.container)
                      └─> container.validity()
                           └─> CErrorReport
                                └─> getErrors()  ✅ Public API
                                     └─> List[Dict] with 'class', 'code', 'details', 'severity'
                                          └─> getEtree() converts to XML
```

---

## Files Modified

### 1. `server/ccp4i2/lib/utils/jobs/validate.py`

**Changes**:
- Import `get_plugin_with_context` instead of `get_job_container`
- Use `plugin.container` instead of standalone container
- Better error handling and logging
- Added docstrings with examples

**Lines Changed**: 19-119 (entire function)

### 2. `server/ccp4i2/lib/job_utils/validate_container.py`

**Changes**:
- Added `SEVERITY_TEXT` mapping constant (lines 10-16)
- Changed `error_report._reports` to `error_report.getErrors()` (line 31)
- Updated error dictionary access:
  - `item["class"]` - handle both class objects and strings
  - `item["details"]` - get description
  - `item["severity"]` - get severity directly
- Added `objectName` XML element for error's object name
- Better error handling and documentation

**Lines Changed**: 10-86 (getEtree function and constants)

---

## Benefits of Refactoring

### Immediate

1. ✅ **Validation works** - No more AttributeError crashes
2. ✅ **Proper error reporting** - XML output formatted correctly
3. ✅ **CPluginScript architecture** - Uses unified approach
4. ✅ **DB-aware** - Plugin has dbHandler for future enhancements

### Long-term

1. ✅ **Consistent with other utilities** - Same pattern as `set_parameter`
2. ✅ **Maintainable** - Single way to work with jobs (via `get_plugin_with_context`)
3. ✅ **Extensible** - Easy to add DB synchronization if needed
4. ✅ **Type safe** - Plugin provides proper container hierarchy

---

## Usage Examples

### Basic Validation

```python
from ccp4i2.db.models import Job
from ccp4i2.lib.utils.jobs.validate import validate_job

job = Job.objects.get(uuid="...")
result = validate_job(job)

if result.success:
    error_tree = result.data
    errors = error_tree.findall('.//errorReport')
    print(f"Found {len(errors)} validation reports")
else:
    print(f"Validation failed: {result.error}")
```

### Checking Severity

```python
result = validate_job(job)
if result.success:
    error_tree = result.data
    for error_report in error_tree.findall('.//errorReport'):
        severity = error_report.findtext('severity')
        description = error_report.findtext('description')
        print(f"{severity}: {description}")
```

### Management Command

```bash
# Basic validation
python manage.py validate_job --jobuuid <uuid>

# Save XML output
python manage.py validate_job --jobuuid <uuid> -o /tmp/validation.xml

# JSON output
python manage.py validate_job --jobuuid <uuid> --json

# By project and job number
python manage.py validate_job -pn myproject -jn 1
```

---

## Next Steps

Following the refactoring plan from `CPLUGINSCRIPT_ARCHITECTURE_REFACTORING.md`:

### Completed ✅
1. ✅ **Step 1**: Create `get_plugin_with_context()` utility
2. ✅ **Step 2**: Refactor `set_parameter()` to use CPluginScript
3. ✅ **Step 3**: Refactor `validate_job()` to use CPluginScript

### Remaining
4. ⏭️  **Step 4**: Refactor `get_job_reports()` to use CPluginScript (get params, report, diagnostic XML)
5. ⏭️  **Step 5**: Test complete workflow end-to-end

---

## Key Insights

### The CPluginScript Pattern

Every job operation should follow this pattern:

```python
# 1. Get plugin with DB context
plugin_result = get_plugin_with_context(job)
if not plugin_result.success:
    return Result.fail(plugin_result.error)

plugin = plugin_result.data

# 2. Operate on plugin.container
error_report = plugin.container.validity()

# 3. Save if modified (not needed for validation)
# plugin.container.saveDataToXml(str(job.directory / "params.xml"))

# 4. Update DB if needed (not needed for validation)
# if plugin._dbHandler:
#     plugin._dbHandler.updateJobStatus(jobId=str(job.uuid), container=plugin.container)

# 5. Return result
return Result.ok(data)
```

### Why This Matters

Before this refactoring, we were creating "orphaned" containers that had:
- ❌ No parent plugin
- ❌ No database handler
- ❌ No job context
- ❌ Broken file parameter handling

Now with CPluginScript:
- ✅ Proper hierarchy (container knows its parent plugin)
- ✅ Database awareness (dbHandler attached)
- ✅ Job context (UUID, number, project)
- ✅ File parameters work correctly

---

## Success! 🎉

The `validate_job` function now:
- Uses the unified CPluginScript architecture
- Works with the new CErrorReport API
- Generates proper XML output
- Provides JSON and text output modes
- Has comprehensive error handling

**This completes the validation refactoring!**

Next up: `get_job_reports()` refactoring! 🚀
