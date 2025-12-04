# Refactoring: Move _get_conversion_output_path to CDataFile Base Class

**Date**: 2025-10-27
**Status**: ✅ COMPLETE

## Overview

Moved the `_get_conversion_output_path()` helper method from **CMiniMtzDataFile** (MTZ-specific) to **CDataFile** (base class) to enable file format conversions for all file types, not just MTZ files.

## Motivation

The original implementation in `CMiniMtzDataFile` was MTZ-specific but the logic was actually generic:
- Calculate output filename: `{input_stem}_as_{CONTENT_TYPE}{extension}`
- Try writing to input directory (if writable), else use work_directory
- Handle permissions gracefully

This same pattern can be used for:
- **CPdbDataFile**: PDB ↔ mmCIF conversion
- **Future file types**: Any format conversion needs

## Changes Made

### 1. Added Method to CDataFile Base Class

**File**: `core/base_object/base_classes.py`

Added generic `_get_conversion_output_path()` with enhanced signature:

```python
def _get_conversion_output_path(
    self,
    target_content_type: str,
    target_extension: Optional[str] = None,  # NEW: Optional extension override
    work_directory: Optional[Any] = None
) -> str:
    """Calculate output path for converted file (generic for all file types).

    Args:
        target_content_type: Name of target format (e.g., 'FMEAN', 'MMCIF')
        target_extension: Optional extension override (e.g., '.cif').
                        If None, uses input file's extension.
        work_directory: Fallback directory if input dir not writable

    Examples:
        # MTZ: Preserves .mtz extension
        obs_file._get_conversion_output_path('FMEAN')
        # → '/data/input_as_FMEAN.mtz'

        # PDB→mmCIF: Changes to .cif extension
        pdb_file._get_conversion_output_path('MMCIF', target_extension='.cif')
        # → '/data/model_as_MMCIF.cif'
    """
```

**Key Enhancement**: Added optional `target_extension` parameter to support conversions that change file extensions (PDB→mmCIF).

### 2. Removed Duplicate from CMiniMtzDataFile

**File**: `core/CCP4XtalData.py` (line 532)

Removed the MTZ-specific version and added comment:
```python
# Note: _get_conversion_output_path() is now in CDataFile base class
```

### 3. Updated All MTZ Conversion Method Calls

Fixed all calls to use keyword argument for `work_directory`:

**Before** (would assign work_directory to target_extension by mistake):
```python
output_path = self._get_conversion_output_path('FMEAN', work_directory)
```

**After**:
```python
output_path = self._get_conversion_output_path('FMEAN', work_directory=work_directory)
```

**Files changed**:
- `CObsDataFile.as_IPAIR()` (line 668)
- `CObsDataFile.as_FPAIR()` (line 696)
- `CObsDataFile.as_IMEAN()` (line 719)
- `CObsDataFile.as_FMEAN()` (line 742)
- `CPhsDataFile.as_HL()` (line 810)
- `CPhsDataFile.as_PHIFOM()` (line 833)

### 4. Added Missing Import

**File**: `core/base_object/base_classes.py` (lines 8, 22)

Added `Optional` to imports:
```python
from typing import Any, Dict, Optional  # Added Optional
```

## Testing

All tests pass:
```bash
$ python -m pytest tests/test_cpluginscript_makehklin.py -v
==================== 19 passed in 0.13s ====================

$ python -m pytest tests/ -v
==================== 151 passed, 26 skipped, 2 warnings in 1.46s ====================
```

No regressions introduced.

## Benefits

### 1. Reusability
The method is now available to **all 44 CDataFile descendants**, including:
- MTZ files (CObsDataFile, CPhsDataFile, CMapCoeffsDataFile, CFreeRDataFile)
- Coordinate files (CPdbDataFile)
- Any future file types needing conversion

### 2. Flexibility
The new `target_extension` parameter supports conversions that change file extensions:

```python
# MTZ: Same extension
obs.as_FMEAN()  # input.mtz → input_as_FMEAN.mtz

# PDB→mmCIF: Different extension
pdb.as_MMCIF()  # model.pdb → model_as_MMCIF.cif (with target_extension='.cif')
```

### 3. Consistency
All file conversions now use the same naming pattern and directory selection logic.

## Usage Examples

### MTZ Conversion (No Extension Change)

```python
class CObsDataFile(CObsDataFileStub):
    def as_FMEAN(self, work_directory: Optional[Any] = None) -> str:
        """Convert to Mean Structure Factors format."""
        output_path = self._get_conversion_output_path('FMEAN', work_directory=work_directory)

        # Conversion logic here...
        # Input:  /data/obs.mtz (IPAIR format)
        # Output: /data/obs_as_FMEAN.mtz (FMEAN format)

        return output_path
```

### PDB→mmCIF Conversion (Extension Change)

```python
class CPdbDataFile(CPdbDataFileStub):
    def as_MMCIF(self, work_directory: Optional[Any] = None) -> str:
        """Convert PDB to mmCIF format."""
        output_path = self._get_conversion_output_path(
            'MMCIF',
            target_extension='.cif',  # Change extension to .cif
            work_directory=work_directory
        )

        import gemmi
        structure = gemmi.read_structure(str(self.getFullPath()))
        gemmi.write_minimal_cif(structure, str(output_path))

        # Input:  /data/model.pdb
        # Output: /data/model_as_MMCIF.cif

        return output_path
```

## Implementation Details

### Directory Selection Logic

The method tries these locations in order:

1. **Input file's directory** (preferred)
   - Tests writability by creating temporary `.write_test_{id}` file
   - Deletes test file immediately
   - Returns this path if successful

2. **Plugin's work_directory** (fallback)
   - Used if input directory not writable
   - Passed as parameter to conversion method

3. **Input directory** (last resort)
   - May fail at actual write time if truly not writable
   - Better to fail with clear error than guess wrong location

### Extension Handling

```python
# Use target_extension if provided, otherwise preserve input extension
extension = target_extension if target_extension else input_path.suffix

# Examples:
# MTZ: target_extension=None → uses '.mtz' from input
# PDB: target_extension='.cif' → uses '.cif' override
```

## Related Documentation

- **MTZ_CONVERSION_SYSTEM.md** - Overall MTZ conversion architecture
- **CONTENT_FLAGS_SUBTYPES_COMPLETE.md** - Complete metadata for all file types
- **MAKEHKLIN_IMPLEMENTATION_COMPLETE.md** - makeHklin system

## Future Work

### CPdbDataFile Conversions

Now that the base helper is available, implement PDB↔mmCIF conversion:

```python
class CPdbDataFile(CPdbDataFileStub):
    def as_PDB(self, work_directory: Optional[Any] = None) -> str:
        """Convert to PDB format."""
        output_path = self._get_conversion_output_path(
            'PDB',
            target_extension='.pdb',
            work_directory=work_directory
        )

        import gemmi
        structure = gemmi.read_structure(str(self.getFullPath()))
        gemmi.write_pdb(structure, str(output_path))
        return output_path

    def as_MMCIF(self, work_directory: Optional[Any] = None) -> str:
        """Convert to mmCIF format."""
        output_path = self._get_conversion_output_path(
            'MMCIF',
            target_extension='.cif',
            work_directory=work_directory
        )

        import gemmi
        structure = gemmi.read_structure(str(self.getFullPath()))
        gemmi.write_minimal_cif(structure, str(output_path))
        return output_path
```

## Summary

✅ **Method moved to CDataFile base class**
✅ **Enhanced with target_extension parameter**
✅ **All MTZ conversion calls updated**
✅ **All tests passing (151 passed)**
✅ **No regressions**
✅ **Ready for CPdbDataFile and future file type conversions**

The refactoring makes the conversion system truly generic and extensible for all file types in the codebase.
