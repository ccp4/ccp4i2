# MTZ File Content Type Conversion System

**Date**: 2025-10-27
**Status**: ✅ STUB METHODS IMPLEMENTED (Conversion logic pending)

## Overview

The makeHklin system now supports **automatic content type conversion** for MTZ files. When a file's current `contentFlag` doesn't match the requested format, the system automatically converts the file to the target format before merging.

## Key Concept: Content Flags as Conversion Requests

In the legacy CCP4i2 API, the syntax `['HKLIN2', target_flag]` has two meanings:

1. **If `file.contentFlag == target_flag`**: Use file as-is (no conversion)
2. **If `file.contentFlag != target_flag`**: Convert file to target format first

This is fundamentally different from just "reading different columns" - it's actual **data transformation**.

## Base Class Helper Method

All file conversion methods use a shared helper in the **CDataFile base class**:

```python
class CDataFile(CData):
    def _get_conversion_output_path(
        self,
        target_content_type: str,
        target_extension: Optional[str] = None,
        work_directory: Optional[Any] = None
    ) -> str:
        """Calculate output path for converted file (generic for all file types).

        Args:
            target_content_type: Name of target format (e.g., 'FMEAN', 'MMCIF')
            target_extension: Optional extension override (e.g., '.cif').
                            If None, uses input file's extension.
            work_directory: Fallback directory if input dir not writable

        Examples:
            # MTZ conversion (preserves .mtz extension)
            obs_file._get_conversion_output_path('FMEAN')
            # → '/data/input_as_FMEAN.mtz'

            # PDB to mmCIF (changes extension)
            pdb_file._get_conversion_output_path('MMCIF', target_extension='.cif')
            # → '/data/model_as_MMCIF.cif'
        """
```

This generic helper is available to **all CDataFile descendants**, including:
- MTZ files (CObsDataFile, CPhsDataFile, etc.)
- Coordinate files (CPdbDataFile for PDB↔mmCIF)
- Any future file conversion needs

## Conversion Methods

Each `CMiniMtzDataFile` subclass that supports multiple content types now has `as_CONTENTTYPE()` methods:

### CObsDataFile (4 conversion methods)

```python
class CObsDataFile(CObsDataFileStub):
    """Observed data file with 4 possible formats."""

    def as_IPAIR(self, work_directory: Optional[Any] = None) -> str:
        """Convert to Anomalous Intensities format.

        Columns: Iplus, SIGIplus, Iminus, SIGIminus
        """

    def as_FPAIR(self, work_directory: Optional[Any] = None) -> str:
        """Convert to Anomalous Structure Factors format.

        Columns: Fplus, SIGFplus, Fminus, SIGFminus
        """

    def as_IMEAN(self, work_directory: Optional[Any] = None) -> str:
        """Convert to Mean Intensities format.

        Columns: I, SIGI
        """

    def as_FMEAN(self, work_directory: Optional[Any] = None) -> str:
        """Convert to Mean Structure Factors format.

        Columns: F, SIGF
        """
```

### CPhsDataFile (2 conversion methods)

```python
class CPhsDataFile(CPhsDataFileStub):
    """Phase data file with 2 possible formats."""

    def as_HL(self, work_directory: Optional[Any] = None) -> str:
        """Convert to Hendrickson-Lattman coefficients.

        Columns: HLA, HLB, HLC, HLD
        """

    def as_PHIFOM(self, work_directory: Optional[Any] = None) -> str:
        """Convert to Phase + Figure of Merit format.

        Columns: PHI, FOM
        """
```

### CMapCoeffsDataFile (1 method - no conversion needed)

```python
class CMapCoeffsDataFile(CMapCoeffsDataFileStub):
    """Map coefficients file - only one format."""

    def as_FPHI(self, work_directory: Optional[Any] = None) -> str:
        """Return path to file (no conversion needed).

        Columns: F, PHI
        """
```

### CFreeRDataFile

No conversion methods needed - only has one content type (`FREER`).

## Output File Naming Convention

When conversion is needed, the output file follows this pattern:

```
{inputroot}_as_{CONTENT_TYPE}.mtz
```

**Examples**:
- Input: `/data/obs.mtz`, Target: `IPAIR` → `/data/obs_as_IPAIR.mtz`
- Input: `/data/native.mtz`, Target: `FMEAN` → `/data/native_as_FMEAN.mtz`

## Output Directory Selection

The conversion methods use this priority:

1. **Input file's directory** (if writable)
   - Test by creating temporary `.write_test_{id}` file
   - Delete test file immediately

2. **Plugin's workDirectory** (fallback)
   - Used if input directory not writable
   - Passed as parameter to conversion method

3. **Input directory** (last resort)
   - May fail at write time if truly not writable

## makeHklin Integration

The `makeHklin` method automatically handles conversions:

```python
def makeHklin(self, miniMtzsIn: list, hklin: str = 'hklin') -> CErrorReport:
    """
    For each [name, target_flag] in miniMtzsIn:

    1. Lookup file object
    2. Check if int(file_obj.contentFlag) == target_flag
    3. If NOT equal:
       a. Map target_flag to method name (e.g., 2 -> 'as_FPAIR')
       b. Call file_obj.as_FPAIR(self.workDirectory)
       c. Get converted file path
       d. Create temporary file object pointing to converted file
       e. Set temp_file_obj.contentFlag = target_flag
       f. Add to inputData with temp name
       g. Use temp name for merge
    4. If equal:
       - Use original file (no conversion)
    """
```

## Conversion Flow Example

### Scenario: Convert FMEAN to FPAIR

```python
# File has contentFlag=4 (FMEAN: F, SIGF)
# But we request contentFlag=2 (FPAIR: Fplus, SIGFplus, Fminus, SIGFminus)

error = self.makeHklin([
    ['HKLIN1', CObsDataFile.CONTENT_FLAG_FPAIR]
])
```

**What happens**:

1. **Detection**: `int(HKLIN1.contentFlag) = 4, target_flag = 2, 4 != 2` → Need conversion

2. **Method Resolution**:
   ```python
   # Map flag 2 to content type name
   target_name = _get_content_flag_name(HKLIN1, 2)  # Returns 'FPAIR'

   # Get method
   method = getattr(HKLIN1, 'as_FPAIR')  # HKLIN1.as_FPAIR()
   ```

3. **Conversion Call**:
   ```python
   converted_path = HKLIN1.as_FPAIR(self.workDirectory)
   # Returns: /workdir/input_as_FPAIR.mtz (or NotImplementedError for now)
   ```

4. **Temporary File Object**:
   ```python
   temp_file = CObsDataFile(parent=self.inputData, name='_converted_HKLIN1_0')
   temp_file.baseName = 'input_as_FPAIR.mtz'
   temp_file.relPath = str(workdir)
   temp_file.contentFlag = CInt(2)  # FPAIR

   self.inputData._converted_HKLIN1_0 = temp_file
   ```

5. **Merge**:
   ```python
   self.makeHklinGemmi(['_converted_HKLIN1_0'], ...)
   ```

6. **Cleanup**:
   ```python
   delattr(self.inputData, '_converted_HKLIN1_0')
   ```

## Content Flag Name Mapping

The `_get_content_flag_name()` helper dynamically maps flag values to names:

```python
def _get_content_flag_name(self, file_obj, content_flag: int) -> str:
    """
    Search file_obj's class for CONTENT_FLAG_* constants.

    For CObsDataFile:
        1 -> 'IPAIR'
        2 -> 'FPAIR'
        3 -> 'IMEAN'
        4 -> 'FMEAN'
    """
    for attr_name in dir(file_obj.__class__):
        if attr_name.startswith('CONTENT_FLAG_'):
            if getattr(file_obj.__class__, attr_name) == content_flag:
                return attr_name.replace('CONTENT_FLAG_', '')

    raise ValueError(f"No content flag name found for value {content_flag}")
```

## Error Codes

New error codes for conversion system:

| Code | Description |
|------|-------------|
| 208 | Conversion method not found on file object |
| 209 | Conversion not yet implemented (NotImplementedError) |
| 210 | Generic conversion error |

## Current Status: Stub Methods

**All conversion methods are currently STUBS** that raise `NotImplementedError`:

```python
def as_IPAIR(self, work_directory: Optional[Any] = None) -> str:
    output_path = self._get_conversion_output_path('IPAIR', work_directory)

    # TODO: Implement actual conversion logic using gemmi
    raise NotImplementedError(
        f"Conversion to IPAIR format not yet implemented. "
        f"Would output to: {output_path}"
    )
```

## Future Implementation: Conversion Logic

When implementing conversion methods, they should:

### 1. Read Input File
```python
import gemmi

input_path = self.getFullPath()
mtz = gemmi.read_mtz_file(str(input_path))
```

### 2. Perform Data Transformation

**Example: FMEAN → FPAIR** (mean SFs to anomalous SFs)
```python
# Get F, SIGF columns
F = mtz.column_with_label('F')
SIGF = mtz.column_with_label('SIGF')

# Transform to anomalous pairs
# This is a simplification - real conversion is more complex
F_data = F.array
SIGF_data = SIGF.array

Fplus = F_data  # Simplified
SIGFplus = SIGF_data
Fminus = F_data
SIGFminus = SIGF_data
```

**Example: FPAIR → FMEAN** (anomalous SFs to mean SFs)
```python
# Get anomalous pair columns
Fplus = mtz.column_with_label('Fplus').array
Fminus = mtz.column_with_label('Fminus').array

# Calculate mean
F = (Fplus + Fminus) / 2.0
SIGF = ...  # Proper error propagation
```

### 3. Create Output MTZ
```python
out_mtz = gemmi.Mtz(with_base=True)
out_mtz.spacegroup = mtz.spacegroup
out_mtz.set_cell_for_all(mtz.cell)

# Copy H, K, L
hkl_data = mtz.array[:, :3]
out_mtz.set_data(hkl_data)

# Add converted columns
out_mtz.add_dataset('converted')
out_mtz.add_column('Fplus', 'G')
out_mtz.add_column('SIGFplus', 'L')
# ... add all required columns

out_mtz.update_reso()
out_mtz.write_to_file(str(output_path))
```

### 4. Return Path
```python
return str(output_path)
```

## Testing

### Test: Conversion Triggered
```python
def test_override_content_flag_triggers_conversion(plugin_script):
    """Test that contentFlag mismatch triggers conversion."""
    # File has contentFlag=4 (FMEAN)
    error = plugin_script.makeHklin([
        ['HKLIN1', CObsDataFile.CONTENT_FLAG_FPAIR]  # Request flag=2 (FPAIR)
    ])

    # Should get NotImplementedError for now
    assert error.maxSeverity() > SEVERITY_OK
    assert 'not yet implemented' in error.report().lower()
```

### Test: No Conversion When Flags Match
```python
def test_no_conversion_when_flags_match(plugin_script):
    """Test that matching flags skip conversion."""
    # File has contentFlag=4 (FMEAN)
    error = plugin_script.makeHklin([
        ['HKLIN1', CObsDataFile.CONTENT_FLAG_FMEAN]  # Same flag - no conversion
    ])

    assert error.maxSeverity() == SEVERITY_OK  # Success
```

## Files Modified

### Created/Modified:
1. ✅ `core/base_object/base_classes.py`
   - Added `_get_conversion_output_path()` helper method to **CDataFile base class**
   - Generic implementation supports all file types with optional `target_extension` parameter
   - Enables future conversions for CPdbDataFile (PDB↔mmCIF) and other file types

2. ✅ `core/CCP4XtalData.py`
   - Added `as_IPAIR()`, `as_FPAIR()`, `as_IMEAN()`, `as_FMEAN()` to `CObsDataFile`
   - Added `as_HL()`, `as_PHIFOM()` to `CPhsDataFile`
   - Added `as_FPHI()` to `CMapCoeffsDataFile`
   - All methods use `_get_conversion_output_path(work_directory=work_directory)` from base class

3. ✅ `core/CCP4PluginScript.py`
   - Updated `makeHklin()` to detect and handle conversions
   - Added `_get_content_flag_name()` helper method
   - Added temporary file object creation for converted files
   - Added cleanup in finally block

4. ✅ `tests/test_cpluginscript_makehklin.py`
   - Updated `test_override_content_flag` → `test_override_content_flag_triggers_conversion`
   - Updated `test_contentflag_restoration` → `test_no_conversion_when_flags_match`
   - Both tests now reflect conversion behavior

## Common Conversions

### Intensity ↔ Structure Factor

**Intensities (I, SIGI) → Structure Factors (F, SIGF)**:
```
F = sqrt(I)
SIGF = SIGI / (2 * sqrt(I))
```

**Structure Factors (F, SIGF) → Intensities (I, SIGI)**:
```
I = F²
SIGI = 2 * F * SIGF
```

### Anomalous Pairs ↔ Mean Values

**Anomalous (Fplus, Fminus) → Mean (F)**:
```
F = (Fplus + Fminus) / 2
SIGF = sqrt((SIGFplus² + SIGFminus²) / 4)
```

**Mean (F) → Anomalous (Fplus, Fminus)**:
- Cannot be done uniquely without additional data
- Requires anomalous signal information
- May create symmetric pairs as approximation

### Hendrickson-Lattman ↔ Phase/FOM

**HL (HLA, HLB, HLC, HLD) → PHIFOM (PHI, FOM)**:
```
# Complex calculation involving HL coefficient algebra
# See Hendrickson-Lattman 1970 paper
```

**PHIFOM (PHI, FOM) → HL (HLA, HLB, HLC, HLD)**:
- Can approximate but loses information
- HLA = FOM * cos(PHI), rest zeros (simplified)

## Architecture Integration

```
Legacy API Call
    ↓
makeHklin
    ↓
[name, target_flag] detected
    ↓
Check: current_flag != target_flag?
    ↓ YES
Conversion Needed
    ↓
1. Map target_flag → method name (as_CONTENTTYPE)
2. Call file_obj.as_CONTENTTYPE(workDirectory)
3. Get converted_path
4. Create temp file object
5. Set contentFlag = target_flag
6. Add to inputData
    ↓
makeHklinGemmi (uses converted file)
    ↓
merge_mtz_files (pure gemmi)
```

## Usage Examples

### Example 1: Convert FMEAN to FPAIR
```python
class MyTask(CPluginScript):
    def process(self):
        # Input file has F, SIGF (contentFlag=4, FMEAN)
        # Request Fplus, SIGFplus, Fminus, SIGFminus (contentFlag=2, FPAIR)

        error = self.makeHklin([
            ['HKLIN1', CObsDataFile.CONTENT_FLAG_FPAIR]
        ])

        # System automatically:
        # 1. Detects contentFlag mismatch (4 != 2)
        # 2. Calls HKLIN1.as_FPAIR(self.workDirectory)
        # 3. Creates /workdir/HKLIN1_as_FPAIR.mtz
        # 4. Uses converted file for merging
```

### Example 2: No Conversion Needed
```python
class MyTask(CPluginScript):
    def process(self):
        # File already has correct format
        error = self.makeHklin([
            ['HKLIN1', CObsDataFile.CONTENT_FLAG_FMEAN]  # Matches file's flag
        ])

        # No conversion - uses original file directly
```

### Example 3: Mixed Conversion
```python
class MyTask(CPluginScript):
    def process(self):
        error = self.makeHklin([
            'HKLIN1',  # Use as-is (contentFlag not specified)
            ['HKLIN2', CObsDataFile.CONTENT_FLAG_FPAIR],  # May convert
            'FREERFLAG'  # Use as-is
        ])
```

## Next Steps

1. **Implement Conversion Logic**: Replace NotImplementedError with actual gemmi-based conversions
2. **Add Unit Tests**: Test each conversion method independently
3. **Crystallographic Validation**: Ensure converted data is scientifically correct
4. **Performance Optimization**: Cache converted files to avoid re-conversion
5. **Documentation**: Add examples to user documentation

## Related Documentation

- `MAKEHKLIN_IMPLEMENTATION_COMPLETE.md` - Overall makeHklin system
- `MERGE_MTZ_REFACTORING_SUMMARY.md` - API design rationale
- `MTZ_METADATA_EXTRACTION_SUMMARY.md` - CONTENT_FLAG system

## Conclusion

The conversion system is **architecturally complete** with stub methods in place. The framework handles:

- ✅ Automatic detection of conversion needs
- ✅ Dynamic method resolution
- ✅ Temporary file object creation
- ✅ Output path calculation
- ✅ Error handling
- ✅ Cleanup of temporary objects

The only remaining work is implementing the actual crystallographic transformation logic in each `as_CONTENTTYPE()` method.
