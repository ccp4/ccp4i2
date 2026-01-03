# makeHklin Implementation - Complete

**Date**: 2025-10-27
**Status**: ✅ COMPLETED

## Summary

Successfully implemented complete MTZ file merging functionality with both modern Pythonic API (`makeHklinGemmi`) and backward-compatible legacy API (`makeHklin`) for CCP4 crystallographic workflows.

## Architecture Overview

```
┌──────────────────────────────────────────────────────────┐
│ Legacy Scripts (CCP4i2)                                  │
│ • Use makeHklin(['HKLIN1', ['HKLIN2', flag]])          │
│ • Backward-compatible with old API                      │
└────────────────┬─────────────────────────────────────────┘
                 │
                 ▼
┌──────────────────────────────────────────────────────────┐
│ CPluginScript.makeHklin (Adapter)                       │
│ • Translates legacy [name, flag] syntax                 │
│ • Temporarily overrides contentFlag attributes          │
│ • Returns CErrorReport for compatibility                │
└────────────────┬─────────────────────────────────────────┘
                 │
                 ▼
┌──────────────────────────────────────────────────────────┐
│ CPluginScript.makeHklinGemmi (High-Level API)           │
│ • Resolves container names to CMiniMtzDataFile objects  │
│ • Reads contentFlag attribute from file objects         │
│ • Uses CONTENT_SIGNATURE_LIST[contentFlag-1] for cols   │
│ • Builds column_mapping dictionaries                    │
│ • Returns Path to output file                           │
└────────────────┬─────────────────────────────────────────┘
                 │
                 ▼
┌──────────────────────────────────────────────────────────┐
│ CCP4Utils.merge_mtz_files (Pure Gemmi Utility)          │
│ • NO CData dependencies                                  │
│ • Pure gemmi-based MTZ manipulation                     │
│ • Takes column_mapping: {input_label -> output_label}   │
│ • Handles conflicts with merge strategies               │
└──────────────────────────────────────────────────────────┘
```

## Key Design Principles

### 1. Separation of Concerns
- **Low-level utility** (`merge_mtz_files`): Pure gemmi operations, completely CData-agnostic
- **High-level integration** (`makeHklinGemmi`): CData-aware, handles contentFlag resolution
- **Backward compatibility** (`makeHklin`): Legacy API wrapper for existing code

### 2. Clean API Design
- Modern API uses explicit `column_mapping` dictionaries
- No implicit column resolution in low-level utilities
- Column logic lives where CData knowledge exists

### 3. Type Safety
- Uses Path objects for file paths
- Strong typing in function signatures
- Clear error reporting with CErrorReport

## Implementation Details

### 1. merge_mtz_files (core/CCP4Utils.py:19-222)

**Pure gemmi utility - NO CData dependencies**

```python
def merge_mtz_files(
    input_specs: List[dict],
    output_path: Union[str, Path],
    merge_strategy: str = 'first'
) -> Path:
    """
    Merge multiple MTZ files using gemmi.

    Args:
        input_specs: [
            {
                'path': str/Path,
                'column_mapping': {input_label: output_label}
            }
        ]
        output_path: Where to write merged file
        merge_strategy: 'first' | 'last' | 'error' | 'rename'

    Returns:
        Path to output file
    """
```

**Features**:
- H,K,L matching between files
- Space group and unit cell validation
- Conflict resolution strategies
- Comprehensive error messages
- MTZ history tracking

### 2. makeHklinGemmi (core/CCP4PluginScript.py:610-734)

**Modern Pythonic API with CData integration**

```python
def makeHklinGemmi(
    self,
    file_objects: list,
    output_name: str = 'hklin',
    merge_strategy: str = 'first'
) -> Path:
    """
    Merge normalized mini-MTZ files (new Pythonic API).

    Args:
        file_objects: List of specs:
            - str: 'HKLIN1' (uses object's contentFlag)
            - dict: {'name': 'HKLIN1', 'rename': {'F': 'F_NAT'}}
        output_name: Base name for output (default: 'hklin')
        merge_strategy: Conflict resolution strategy

    Returns:
        Path to output file
    """
```

**Example Usage**:
```python
# Simple merge
hklin = self.makeHklinGemmi(['HKLIN1', 'FREERFLAG'])

# With column renaming
hklin = self.makeHklinGemmi([
    {
        'name': 'HKLIN1',
        'rename': {'F': 'F_native', 'SIGF': 'SIGF_native'}
    },
    'HKLIN2'
])
```

**Column Resolution Logic**:
```python
# For each file object:
file_obj = getattr(self.inputData, 'HKLIN1')
content_flag = int(file_obj.contentFlag)  # e.g., 4 (FMEAN)
columns = file_obj.CONTENT_SIGNATURE_LIST[content_flag - 1]  # ['F', 'SIGF']

# Build mapping
column_mapping = {
    'F': 'F_native',      # With rename
    'SIGF': 'SIGF_native'
}
```

### 3. makeHklin (core/CCP4PluginScript.py:736-844)

**Backward-compatible legacy API wrapper**

```python
def makeHklin(
    self,
    miniMtzsIn: list,
    hklin: str = 'hklin'
) -> CErrorReport:
    """
    Merge mini-MTZ files (backward-compatible legacy API).

    Args:
        miniMtzsIn: List of:
            - str: 'HKLIN1' (uses object's contentFlag)
            - [str, int]: ['HKLIN2', explicit_contentFlag]
        hklin: Output base name

    Returns:
        CErrorReport (empty if successful)
    """
```

**Example Usage** (old CCP4i2 code):
```python
# Simple merge
error = self.makeHklin(['HKLIN1', 'FREERFLAG'])

# Override contentFlag
error = self.makeHklin([
    'HKLIN1',
    ['HKLIN2', CObsDataFile.CONTENT_FLAG_FPAIR]  # Override to FPAIR
])
```

**ContentFlag Override Implementation**:
```python
# Build list of overrides
overrides = []  # (file_obj, original_flag)

for item in miniMtzsIn:
    if isinstance(item, [str, int]):
        name, explicit_flag = item
        file_obj = getattr(self.inputData, name)

        # Save and override
        overrides.append((file_obj, file_obj.contentFlag))
        file_obj.contentFlag = explicit_flag

try:
    # Call new API with overrides active
    self.makeHklinGemmi(file_objects, hklin, 'first')
finally:
    # Restore all original contentFlags
    for file_obj, original_flag in overrides:
        file_obj.contentFlag = original_flag
```

## Testing

### Test Suite: tests/test_cpluginscript_makehklin.py

**19 comprehensive integration tests covering**:

#### makeHklinGemmi Tests (9 tests)
1. ✅ `test_simple_merge_with_string_names` - Basic merge with contentFlag resolution
2. ✅ `test_merge_with_rename` - Column renaming functionality
3. ✅ `test_custom_output_name` - Custom output file naming
4. ✅ `test_merge_strategy_error` - Error on column conflicts
5. ✅ `test_merge_strategy_rename` - Auto-renaming conflicts (F, F_1, F_2)
6. ✅ `test_file_not_found_error` - Error when container name not found
7. ✅ `test_no_content_signature_list_error` - Error for non-MTZ objects
8. ✅ `test_invalid_content_flag_error` - Invalid contentFlag handling
9. ✅ `test_no_path_error` - Error when file has no path set

#### makeHklin Tests (7 tests)
10. ✅ `test_simple_merge_with_string_names` - Legacy API basic merge
11. ✅ `test_custom_output_name` - Legacy API custom naming
12. ✅ `test_override_content_flag` - ContentFlag override with [name, flag]
13. ✅ `test_file_not_found_error` - Error handling in legacy API
14. ✅ `test_invalid_content_flag_error` - Invalid flag in override
15. ✅ `test_invalid_item_format_error` - Invalid item format detection
16. ✅ `test_contentflag_restoration` - Verify flags restored after override

#### Edge Case Tests (3 tests)
17. ✅ `test_empty_file_list` - Error on empty file list
18. ✅ `test_single_file_merge` - Single file "merge"
19. ✅ `test_partial_rename` - Rename only some columns

**Test Results**: All 19 tests passing
**Full Suite**: 151 tests passing, 26 skipped

### Mock Objects for Testing

```python
class MockCObsDataFile(MockMiniMtzDataFile):
    """Mock with real CONTENT_SIGNATURE_LIST metadata."""

    CONTENT_FLAG_IPAIR = 1  # Anomalous Is
    CONTENT_FLAG_FPAIR = 2  # Anomalous SFs
    CONTENT_FLAG_IMEAN = 3  # Mean Is
    CONTENT_FLAG_FMEAN = 4  # Mean SFs

    CONTENT_SIGNATURE_LIST = [
        ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'],
        ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus'],
        ['I', 'SIGI'],
        ['F', 'SIGF']
    ]
```

## Files Modified/Created

### Created
1. ✅ `core/CCP4Utils.py` (223 lines) - Pure gemmi utility
2. ✅ `tests/test_merge_mtz_files.py` (270 lines) - 10 low-level tests
3. ✅ `tests/test_cpluginscript_makehklin.py` (441 lines) - 19 integration tests
4. ✅ `MERGE_MTZ_REFACTORING_SUMMARY.md` - API refactoring documentation
5. ✅ `MAKEHKLIN_IMPLEMENTATION_COMPLETE.md` - This document

### Modified
1. ✅ `core/CCP4PluginScript.py` - Added makeHklinGemmi (lines 610-734) and makeHklin (lines 736-844)
2. ✅ `migration/CData/cdata.json` - Added MTZ metadata (CONTENT_FLAGS, SUBTYPES)
3. ✅ `migration/CData/stub_generator.py` - Enhanced with MTZ metadata rendering
4. ✅ `core/cdata_stubs/CCP4XtalData.py` - Regenerated with class constants

## Usage Examples

### Example 1: Simple Merge (Modern API)
```python
class MyRefineTask(CPluginScript):
    def process(self):
        # Merge observed data with FreeR flags
        hklin_path = self.makeHklinGemmi([
            'HKLIN1',      # CObsDataFile with F, SIGF
            'FREERFLAG'    # CFreeRDataFile with FREER
        ])

        # hklin_path = /workdir/hklin.mtz
        # Contains: H, K, L, F, SIGF, FREER
```

### Example 2: Column Renaming (Modern API)
```python
class MyTask(CPluginScript):
    def process(self):
        # Merge native and derivative with renamed columns
        hklin_path = self.makeHklinGemmi([
            {
                'name': 'HKLIN_NATIVE',
                'rename': {'F': 'F_NAT', 'SIGF': 'SIGF_NAT'}
            },
            {
                'name': 'HKLIN_DERIV',
                'rename': {'F': 'F_DER', 'SIGF': 'SIGF_DER'}
            },
            'FREERFLAG'
        ])

        # Output contains:
        # H, K, L, F_NAT, SIGF_NAT, F_DER, SIGF_DER, FREER
```

### Example 3: ContentFlag Override (Legacy API)
```python
class MyTask(CPluginScript):
    def process(self):
        # Legacy code that overrides contentFlag
        error = self.makeHklin([
            'HKLIN1',  # Uses HKLIN1.contentFlag (e.g., 4=FMEAN)
            ['HKLIN2', CObsDataFile.CONTENT_FLAG_FPAIR]  # Force FPAIR format
        ])

        if error.maxSeverity() > SEVERITY_OK:
            print(f"Error: {error.report()}")
```

### Example 4: Conflict Resolution
```python
class MyTask(CPluginScript):
    def process(self):
        # Handle column conflicts with auto-renaming
        hklin_path = self.makeHklinGemmi(
            ['HKLIN1', 'HKLIN2', 'HKLIN3'],  # All have F, SIGF
            merge_strategy='rename'  # Creates F, F_1, F_2
        )
```

## MTZ Metadata Integration

### Class Constants in Stubs

All CMiniMtzDataFile subclasses now have class-level constants:

```python
# In core/cdata_stubs/CCP4XtalData.py

class CObsDataFileStub(CMiniMtzDataFileStub):
    """An MTZ experimental data file"""

    # Subtype constants
    SUBTYPE_OBSERVED = 1  # observed data
    SUBTYPE_DERIVED = 2   # derived data
    SUBTYPE_REFERENCE = 3 # reference data

    # Content flag constants
    CONTENT_FLAG_IPAIR = 1  # Anomalous Is
    CONTENT_FLAG_FPAIR = 2  # Anomalous SFs
    CONTENT_FLAG_IMEAN = 3  # Mean Is
    CONTENT_FLAG_FMEAN = 4  # Mean SFs

    # Metadata arrays (0-indexed)
    CONTENT_ANNOTATION = ['Anomalous Is', 'Anomalous SFs', 'Mean Is', 'Mean SFs']
    CONTENT_SIGNATURE_LIST = [
        ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'],
        ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus'],
        ['I', 'SIGI'],
        ['F', 'SIGF']
    ]
```

### Usage in Code
```python
# Get columns for a file object
file_obj = self.inputData.HKLIN1  # CObsDataFile instance
content_flag = int(file_obj.contentFlag)  # 4 (FMEAN)
columns = file_obj.CONTENT_SIGNATURE_LIST[content_flag - 1]  # ['F', 'SIGF']
annotation = file_obj.CONTENT_ANNOTATION[content_flag - 1]  # 'Mean SFs'
```

## Error Handling

### Error Reporting with CErrorReport
```python
error = self.makeHklin(['HKLIN1', 'NONEXISTENT'])

if error.maxSeverity() > SEVERITY_OK:
    # Print formatted error report
    print(error.report())
    # ERROR in CPluginScript 'hklin': File object 'NONEXISTENT' not found (code 205)

    # Get error count
    print(f"Found {error.count()} errors")

    # Get raw error list
    for err in error.getErrors():
        print(f"Class: {err['class']}, Code: {err['code']}")
```

### Common Error Codes

| Code | Meaning |
|------|---------|
| 200 | General MTZ merging error |
| 205 | File object not found in containers |
| 206 | Invalid contentFlag value |
| 207 | Invalid miniMtzsIn item format |

### Exceptions from merge_mtz_files
```python
try:
    result = merge_mtz_files(input_specs, output_path)
except FileNotFoundError as e:
    print(f"MTZ file not found: {e}")
except ValueError as e:
    print(f"Column or validation error: {e}")
except MtzMergeError as e:
    print(f"Gemmi operation failed: {e}")
```

## Performance Characteristics

- **File I/O**: One read per input file, one write for output
- **Memory**: Loads full reflection arrays into memory (gemmi)
- **Column Matching**: Gemmi handles H,K,L matching efficiently
- **Typical Use**: 2-3 input files with 10k-100k reflections

## Future Enhancements

### Potential Improvements
1. **Streaming**: For very large MTZ files (>1M reflections)
2. **Parallel Processing**: Multi-threaded column copying
3. **Advanced Merging**: Wavelength-aware merging for MAD data
4. **Validation**: Check column types match expected types
5. **History Preservation**: Merge history from input files

### Extension Points
- Custom merge strategies via callback functions
- Column transformation functions (scaling, etc.)
- Flexible output formats (beyond MTZ)

## Lessons Learned

### Key Design Decisions
1. **CData Agnosticism**: Keeping low-level utilities pure makes them reusable and testable
2. **ContentFlag Override**: Deferred restoration until after API call prevents premature cleanup
3. **Mock Testing**: Creating mock objects with real metadata structure enables comprehensive testing
4. **Error API**: Using CErrorReport for legacy compatibility while providing modern exceptions

### Challenges Overcome
1. **gemmi API**: Required numpy arrays, not lists
2. **ContentFlag Timing**: Initially restored too early in try/finally blocks
3. **CErrorReport API**: No `isNoError()` method, used `maxSeverity()` instead
4. **Import Paths**: CInt in fundamental_types, not base_classes

## Related Documentation

- `MERGE_MTZ_REFACTORING_SUMMARY.md` - API refactoring details
- `MAKEHKLIN_REFACTOR_PLAN.md` - Original implementation plan
- `MTZ_METADATA_EXTRACTION_SUMMARY.md` - Metadata extraction process
- `API_HARMONIZATION_PLAN.md` - Overall API design philosophy

## Verification Checklist

- ✅ All 19 integration tests passing
- ✅ Full test suite passing (151 tests)
- ✅ Both modern and legacy APIs functional
- ✅ Column resolution uses CONTENT_SIGNATURE_LIST correctly
- ✅ ContentFlag override works and restores properly
- ✅ Error handling comprehensive
- ✅ Documentation complete with examples
- ✅ No CData dependencies in low-level utility
- ✅ Backward compatibility maintained
- ✅ Code follows project conventions

## Conclusion

The makeHklin implementation is **complete and production-ready**. It provides:

1. **Clean Architecture**: Proper separation between low-level utilities and high-level integration
2. **Dual APIs**: Modern Pythonic API + backward-compatible legacy wrapper
3. **Comprehensive Testing**: 19 integration tests + 10 low-level utility tests
4. **Type Safety**: Strong typing, clear error messages
5. **Maintainability**: Well-documented, follows project patterns

The implementation successfully bridges old CCP4i2 code with modern Python practices while maintaining full backward compatibility.
