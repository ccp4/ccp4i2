# merge_mtz_files API Refactoring Summary

**Date**: 2025-10-27
**Status**: ✅ COMPLETED

## Changes Made

### API Redesign - CData Agnosticism

Refactored `merge_mtz_files()` in `core/CCP4Utils.py` to be **completely agnostic** to CMiniMtzDataFile, CONTENT_SIGNATURE_LIST, and all CData conventions.

#### Before (CData-aware):
```python
merge_mtz_files(
    input_specs=[
        {
            'path': '/data/native.mtz',
            'columns': ['F', 'SIGF'],              # Just column names
            'rename': {'F': 'F_NAT', 'SIGF': 'SIGF_NAT'}  # Separate rename dict
        }
    ],
    output_path='/data/merged.mtz'
)
```

#### After (CData-agnostic):
```python
merge_mtz_files(
    input_specs=[
        {
            'path': '/data/native.mtz',
            'column_mapping': {                    # Unified mapping
                'F': 'F_NAT',                      # input_label -> output_label
                'SIGF': 'SIGF_NAT'
            }
        }
    ],
    output_path='/data/merged.mtz'
)
```

### Key Benefits

1. **Pure Utility** - No knowledge of CMiniMtzDataFile, CONTENT_FLAGS, or CONTENT_SIGNATURE_LIST
2. **Cleaner API** - Single `column_mapping` dict instead of separate `columns` and `rename`
3. **Explicit Mapping** - Each input column explicitly mapped to output name
4. **Separation of Concerns** - Column resolution logic moves to `makeHklinGemmi` where it belongs

### Where Column Resolution Happens

**Now (Correct):**
```python
# In makeHklinGemmi (CPluginScript method):
file_obj = self.inputData.HKLIN1  # CObsDataFile instance
content_flag = int(file_obj.contentFlag)
columns = file_obj.CONTENT_SIGNATURE_LIST[content_flag - 1]
# columns = ['F', 'SIGF']  # From class metadata

# Build column_mapping for merge_mtz_files
column_mapping = {col: col for col in columns}  # Or with renaming
```

**Not in merge_mtz_files** - It just receives the final `input_label -> output_label` mapping.

## Files Modified

### 1. `core/CCP4Utils.py`
- Changed `input_specs` format from `{path, columns, rename}` to `{path, column_mapping}`
- Removed `CONTENT_FLAG_COLUMNS_BY_CLASS` constant (moved logic to makeHklinGemmi)
- Updated history generation to reflect new API
- Updated docstrings with new examples

### 2. `tests/test_merge_mtz_files.py`
- Updated all 10 tests to use `column_mapping` instead of `columns` + `rename`
- All tests passing (10/10)

## Test Results

```
============================= test session starts ==============================
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_merge_two_files_simple PASSED
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_merge_with_rename PASSED
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_merge_strategy_error_on_conflict PASSED
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_merge_strategy_first_keeps_first PASSED
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_merge_strategy_rename_auto_renames PASSED
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_file_not_found PASSED
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_column_not_found PASSED
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_empty_input_specs PASSED
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_missing_required_keys PASSED
tests/test_merge_mtz_files.py::TestMergeMtzFiles::test_preserves_spacegroup_and_cell PASSED

============================== 10 passed in 0.05s ===============================

All project tests: 132 passed, 26 skipped, 2 warnings
```

## Architecture Diagram

```
┌─────────────────────────────────────────────────┐
│ CPluginScript.makeHklinGemmi                    │
│                                                 │
│ • Resolves container.HKLIN1 -> CObsDataFile    │
│ • Reads contentFlag attribute                   │
│ • Looks up CONTENT_SIGNATURE_LIST[flag-1]       │
│ • Builds column_mapping dict                    │
│ • Calls merge_mtz_files with paths + mapping   │
└─────────────────────┬───────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────┐
│ core/CCP4Utils.merge_mtz_files                  │
│                                                 │
│ • Pure gemmi utility (NO CData knowledge)       │
│ • Takes: path + column_mapping dict             │
│ • Uses gemmi to merge MTZ files                 │
│ • Returns: Path to output file                  │
└─────────────────────────────────────────────────┘
```

## Next Steps

1. ✅ `merge_mtz_files` refactored (DONE)
2. ⏭️ Implement `makeHklinGemmi` in CPluginScript
   - Resolve container names to CMiniMtzDataFile instances
   - Extract `contentFlag` attribute
   - Use `instance.CONTENT_SIGNATURE_LIST[int(contentFlag) - 1]` to get columns
   - Build `column_mapping` dictionaries
   - Call `merge_mtz_files`
3. ⏭️ Implement `makeHklin` backward-compatible wrapper
4. ⏭️ Integration tests

## Example Usage (Planned)

```python
# In makeHklinGemmi:
def makeHklinGemmi(self, file_objects, output_name='hklin'):
    input_specs = []

    for file_spec in file_objects:
        # Resolve container name to CMiniMtzDataFile instance
        if isinstance(file_spec, str):
            file_obj = getattr(self.inputData, file_spec)

            # Get columns from contentFlag
            content_flag = int(file_obj.contentFlag)
            columns = file_obj.CONTENT_SIGNATURE_LIST[content_flag - 1]

            # Build column_mapping (identity mapping by default)
            column_mapping = {col: col for col in columns}
        else:
            # Dict with explicit columns/renaming
            file_obj = getattr(self.inputData, file_spec['name'])
            column_mapping = file_spec.get('column_mapping', ...)

        input_specs.append({
            'path': file_obj.getFullPath(),
            'column_mapping': column_mapping
        })

    # Call pure utility
    output_path = merge_mtz_files(
        input_specs=input_specs,
        output_path=self.workDirectory / f"{output_name}.mtz"
    )
    return output_path
```

## Rationale

This refactoring enforces **clean separation of concerns**:

- **Low-level (merge_mtz_files)**: Pure MTZ manipulation with gemmi
- **High-level (makeHklinGemmi)**: CData integration, contentFlag resolution, CONTENT_SIGNATURE_LIST lookup

This makes the code:
- **More testable** - Can test merge_mtz_files without any CData infrastructure
- **More reusable** - merge_mtz_files can be used by non-CData code
- **More maintainable** - Clear boundary between layers
- **Easier to understand** - Each function has one clear responsibility
