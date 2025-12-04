# makeHklin Refactoring Plan

## Overview

Refactor the CPluginScript `makeHklin` method to use the gemmi crystallographic toolkit, with two implementations:
1. **makeHklinGemmi** - New, Pythonic API (in CPluginScript)
2. **makeHklin** - Backward-compatible with old syntax (in CPluginScript)
3. **merge_mtz_files** - Low-level gemmi utility (in core/CCP4Utils.py)

## Architecture Understanding (REVISED)

### Data Flow

```
.def.xml defines containers:
    <file name="HKLIN1" class="CObsDataFile" contentFlag="1" .../>
    <file name="FREERFLAG" class="CFreeRDataFile" contentFlag="2" .../>
           ↓
Plugin loads into container:
    self.inputData.HKLIN1 (CObsDataFile instance)
    self.inputData.FREERFLAG (CFreeRDataFile instance)
           ↓
User/code sets filesystem paths:
    self.inputData.HKLIN1.setFullPath("/data/native.mtz")
    self.inputData.FREERFLAG.setFullPath("/data/free.mtz")
           ↓
makeHklin called with container names:
    self.makeHklin(['HKLIN1', 'FREERFLAG'])
           ↓
Method resolves:
    1. Lookup: self.inputData.HKLIN1 → CObsDataFile object
    2. contentFlag: HKLIN1.contentFlag → CONTENT_FLAG_FCOL (value 1)
    3. Columns: CONTENT_FLAG_FCOL → ['F', 'SIGF'] (normalized/standard)
    4. Path: HKLIN1.getFullPath() → "/data/native.mtz"
           ↓
Calls low-level utility:
    merge_mtz_files([
        {'path': '/data/native.mtz', 'columns': ['F', 'SIGF']},
        {'path': '/data/free.mtz', 'columns': ['FreeR_flag']}
    ])
           ↓
Uses gemmi to merge MTZ files
```

### Key Classes

- **CMiniMtzDataFile** - Base class for mini-MTZ files
  - Subclasses: `CObsDataFile`, `CFreeRDataFile`, `CAnomalousDataFile`, etc.
  - Inherits from: `CDataFile` (which provides `getFullPath()`, `setFullPath()`)
  - Has attribute: `contentFlag` (integer indicating data type)

- **contentFlag** - Metadata indicating standard column layout
  - Each flag value implies specific normalized column names
  - Files are (or will be) normalized to have these exact columns
  - Example: `CONTENT_FLAG_FCOL = 1` → file MUST have `['F', 'SIGF']` columns

## Current (Old) API Analysis

From `rstdocs/source/developers/pipelines.rst` (line 806+):

```python
makeHklin(miniMtzsIn, hklin='hklin')
```

**Parameters:**
- `miniMtzsIn`: A list of either:
  - String: name of attribute in `container.inputData` (e.g., 'HKLIN1')
    - Uses the object's own `contentFlag` to determine columns
  - Sublist: `[name, CONTENT_FLAG_xxxx]` - name + explicit contentFlag
    - Overrides the object's `contentFlag` for column selection
- `hklin`: Output MTZ basename (default: 'hklin')

**Issues with old API:**
- ❌ Tortuous mixing of strings and sublists in same argument
- ❌ Content flags are magic constants (CONTENT_FLAG_xxxx)
- ❌ No explicit control over column renaming
- ❌ Hard to understand without extensive documentation

## New Architecture: Three Layers

### Layer 1: Low-Level Utility (core/CCP4Utils.py)

```python
def merge_mtz_files(
    input_specs: List[dict],
    output_path: Union[str, Path],
    merge_strategy: str = 'first'
) -> Path:
    """
    Low-level MTZ merging using gemmi (no CData dependencies).

    Pure utility function that knows nothing about CPluginScript,
    containers, or contentFlags. Just merges MTZ files.

    Args:
        input_specs: List of dicts with:
            {
                'path': str/Path,              # filesystem path to MTZ
                'columns': List[str],          # which columns to copy
                'rename': Dict[str, str]       # optional column renaming
            }
        output_path: Where to write merged MTZ
        merge_strategy: How to handle column conflicts
            - 'first': Keep column from first file (default)
            - 'last': Keep column from last file
            - 'error': Raise error on conflicts
            - 'rename': Auto-rename conflicts (F, F_1, F_2, ...)

    Returns:
        Path: Full path to created file

    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If column conflict and strategy='error'
        GemmiError: If gemmi operations fail

    Example:
        >>> merge_mtz_files(
        ...     input_specs=[
        ...         {
        ...             'path': '/data/native.mtz',
        ...             'columns': ['F', 'SIGF'],
        ...             'rename': {'F': 'F_NAT', 'SIGF': 'SIGF_NAT'}
        ...         },
        ...         {
        ...             'path': '/data/free.mtz',
        ...             'columns': ['FreeR_flag']
        ...         }
        ...     ],
        ...     output_path='/data/merged.mtz',
        ...     merge_strategy='first'
        ... )
        Path('/data/merged.mtz')
    """
```

### Layer 2: New Pythonic API (CPluginScript.makeHklinGemmi)

```python
def makeHklinGemmi(
    self,
    file_objects: List[Union[str, dict]],
    output_name: str = 'hklin',
    merge_strategy: str = 'first'
) -> Path:
    """
    Merge normalized mini-MTZ files from container (new Pythonic API).

    Works with container attribute names, resolves to CDataFile objects,
    uses contentFlag to determine which columns to copy.

    Args:
        file_objects: List of specifications, each is either:

            - str: Attribute name in self.inputData or self.outputData
                   Example: 'HKLIN1'
                   → Looks up self.inputData.HKLIN1
                   → Uses its contentFlag to determine columns automatically
                   → Extracts filesystem path via getFullPath()

            - dict: Explicit specification
                   {
                       'name': str,              # attribute name (required)
                       'columns': List[str],     # optional - override contentFlag
                       'rename': Dict[str, str]  # optional - rename in output
                   }

        output_name: Base name for output (no extension, no path)
                     Writes to: self.workDirectory / f"{output_name}.mtz"

        merge_strategy: How to handle column conflicts

    Returns:
        Path: Full path to created HKLIN file

    Raises:
        AttributeError: If container attribute not found
        ValueError: If contentFlag unknown or column conflict
        FileNotFoundError: If MTZ file doesn't exist

    Example 1 - Simple merge using contentFlags:
        >>> # HKLIN1 is CObsDataFile with contentFlag=CONTENT_FLAG_FCOL
        >>> # FREERFLAG is CFreeRDataFile with contentFlag=CONTENT_FLAG_FREER
        >>> hklin = self.makeHklinGemmi(['HKLIN1', 'FREERFLAG'])
        >>> # Result: merged.mtz with columns ['H','K','L','F','SIGF','FreeR_flag']

    Example 2 - Override columns:
        >>> hklin = self.makeHklinGemmi([
        ...     {'name': 'HKLIN1'},  # uses contentFlag automatically
        ...     {
        ...         'name': 'HKLIN2',
        ...         'columns': ['F', 'SIGF'],
        ...         'rename': {'F': 'F_derivative', 'SIGF': 'SIGF_derivative'}
        ...     }
        ... ])

    Example 3 - Complex with anomalous data:
        >>> hklin = self.makeHklinGemmi(
        ...     file_objects=[
        ...         'HKLIN_NATIVE',      # F, SIGF
        ...         'HKLIN_DERIVATIVE',  # F(+), SIGF(+), F(-), SIGF(-)
        ...         'FREER_FLAG'         # FreeR_flag
        ...     ],
        ...     output_name='phaser_input',
        ...     merge_strategy='error'  # fail if column conflicts
        ... )
    """
```

### Layer 3: Backward-Compatible API (CPluginScript.makeHklin)

```python
def makeHklin(
    self,
    miniMtzsIn: List[Union[str, List]],
    hklin: str = 'hklin'
) -> CErrorReport:
    """
    Merge mini-MTZ files (backward-compatible legacy API).

    Wrapper around makeHklinGemmi that provides backward compatibility
    with the old CCP4i2 API. Converts old syntax to new syntax internally.

    Args:
        miniMtzsIn: List of either:
            - str: Attribute name in self.inputData
                   → Uses object's own contentFlag

            - [str, int]: [attribute_name, explicit_contentFlag]
                   → Overrides object's contentFlag for column selection

        hklin: Base name for output file (no extension)

    Returns:
        CErrorReport: Error report (empty if successful)

    Example 1 - Simple (uses objects' contentFlags):
        >>> error = self.makeHklin(['HKLIN1', 'FREERFLAG'])

    Example 2 - Override HKLIN2's contentFlag:
        >>> error = self.makeHklin([
        ...     'HKLIN1',                           # uses HKLIN1.contentFlag
        ...     ['HKLIN2', CONTENT_FLAG_IOBS]       # treat as intensities
        ... ])
    """
```

## contentFlag System ✅ EXTRACTED

### Standard Column Mappings

**Extracted from `/Users/nmemn/Developer/ccp4i2/core/CCP4XtalData.py`**

The system defines four CMiniMtzDataFile subclasses, each with their own content flags:

#### CObsDataFile - Observed Crystallographic Data

```python
# Content flags for observed data
CONTENT_FLAG_IPAIR = 1  # Anomalous Is
CONTENT_FLAG_FPAIR = 2  # Anomalous SFs
CONTENT_FLAG_IMEAN = 3  # Mean Is
CONTENT_FLAG_FMEAN = 4  # Mean SFs

# Standard column names for each content flag
CONTENT_SIGNATURE_LIST = [
    ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'],  # IPAIR
    ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus'],  # FPAIR
    ['I', 'SIGI'],                                  # IMEAN
    ['F', 'SIGF']                                   # FMEAN
]

# Subtypes
SUBTYPE_OBSERVED = 1   # observed data
SUBTYPE_DERIVED = 2    # derived data
SUBTYPE_REFERENCE = 3  # reference data
```

#### CPhsDataFile - Phase Data

```python
# Content flags for phase data
CONTENT_FLAG_HL = 1      # Hendrickson-Lattmann coeffs
CONTENT_FLAG_PHIFOM = 2  # Phi,FOM

# Standard column names
CONTENT_SIGNATURE_LIST = [
    ['HLA', 'HLB', 'HLC', 'HLD'],  # HL
    ['PHI', 'FOM']                  # PHIFOM
]

# Subtypes
SUBTYPE_UNBIASED = 1  # unbiased data
SUBTYPE_BIASED = 2    # biased data
```

#### CMapCoeffsDataFile - Map Coefficients

```python
# Content flags for map coefficients
CONTENT_FLAG_FPHI = 1  # FPhi

# Standard column names
CONTENT_SIGNATURE_LIST = [
    ['F', 'PHI']  # FPHI
]

# Subtypes
SUBTYPE_NORMAL = 1           # normal map
SUBTYPE_DIFFERENCE = 2       # difference map
SUBTYPE_ANOM_DIFFERENCE = 3  # anomalous difference map
```

#### CFreeRDataFile - Free R Flags

```python
# Content flags (implicit, only one type)
CONTENT_FLAG_FREER = 1  # FreeR

# Standard column names
CONTENT_SIGNATURE_LIST = [
    ['FREER']  # FREER
]

# No subtypes defined
```

### Complete Mapping Table

This metadata has been added to `cdata.json` and is available for code generation:

```python
# For makeHklin implementation, create a lookup by class
CONTENT_FLAG_COLUMNS_BY_CLASS = {
    'CObsDataFile': {
        1: ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'],
        2: ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus'],
        3: ['I', 'SIGI'],
        4: ['F', 'SIGF']
    },
    'CPhsDataFile': {
        1: ['HLA', 'HLB', 'HLC', 'HLD'],
        2: ['PHI', 'FOM']
    },
    'CMapCoeffsDataFile': {
        1: ['F', 'PHI']
    },
    'CFreeRDataFile': {
        1: ['FREER']
    }
}
```

## Implementation Strategy

### Phase 0: Research contentFlags ✅ COMPLETED

**Status: COMPLETED**

1. ✅ **Extracted from old ccp4i2**:
   - All CONTENT_FLAG constant values for 4 CMiniMtzDataFile subclasses
   - Complete column name mappings (CONTENT_SIGNATURE_LIST)
   - CMiniMtzDataFile subclass list: CObsDataFile, CPhsDataFile, CMapCoeffsDataFile, CFreeRDataFile
   - SUBTYPE constants and descriptions

2. ✅ **Added to cdata.json**:
   - All classes now have `CONTENT_FLAGS` metadata with values, annotations, and column lists
   - All classes now have `SUBTYPES` metadata (where applicable)
   - See script: `migration/CData/add_mtz_metadata.py`

3. ⚠️ **Normalization** (STILL NEEDS CLARIFICATION):
   - Are files already normalized when passed to makeHklin?
   - If not, should normalizer be implemented before makeHklin?
   - User mentioned: "These files will have been normalised (maybe something we need to do soon)"

### Phase 1: Core Utility (merge_mtz_files) ✅ COMPLETED

Location: `core/CCP4Utils.py`

**Status**: ✅ Implemented and tested (10/10 tests passing)

**API Design**:
```python
def merge_mtz_files(
    input_specs: List[dict],  # [{'path': str, 'column_mapping': {in: out}}]
    output_path: Union[str, Path],
    merge_strategy: str = 'first'
) -> Path
```

**Key Features**:
1. ✅ **Completely CData-agnostic** - No knowledge of CMiniMtzDataFile
2. ✅ **Column mapping** - Explicit input_label -> output_label mapping
3. ✅ **Conflict resolution** - Strategies: first, last, error, rename
4. ✅ **Validation** - Space group and cell compatibility checks
5. ✅ **H,K,L matching** - Gemmi handles reflection matching across files
6. ✅ **Comprehensive error handling** - FileNotFoundError, ValueError, MtzMergeError

**Implementation Details**:
1. ✅ **Input Processing** - Validates required keys, checks file existence
2. ✅ **MTZ Reading** - Uses gemmi.read_mtz_file(), extracts metadata from first file
3. ✅ **Column Merging** - Creates output with H,K,L, copies columns with gemmi.copy_column()
4. ✅ **Output Generation** - Updates resolution, adds history, writes to file

**Tests**: `tests/test_merge_mtz_files.py` (10 tests, all passing)

### Phase 2: CPluginScript Integration (makeHklinGemmi)

Location: `core/CCP4PluginScript.py`

**Status**: ⏭️ NEXT TO IMPLEMENT

**API Design**:
```python
def makeHklinGemmi(
    self,
    file_objects: List[Union[str, dict]],
    output_name: str = 'hklin',
    merge_strategy: str = 'first'
) -> Path:
    """Merge normalized mini-MTZ files (new Pythonic API)."""
```

**Implementation Steps**:

1. **Container Lookup Helper**
   ```python
   def _lookup_file_object(self, name: str) -> CDataFile:
       """Look up file object in inputData or outputData."""
       # Use existing container.child_name attribute access
       if hasattr(self.inputData, name):
           return getattr(self.inputData, name)
       elif hasattr(self.outputData, name):
           return getattr(self.outputData, name)
       else:
           raise AttributeError(f"No file '{name}' in inputData/outputData")
   ```

2. **Column Resolution from contentFlag** ⭐ KEY CHANGE
   ```python
   for spec in file_objects:
       # Lookup container object
       file_obj = self._lookup_file_object(name)

       # Get filesystem path
       path = file_obj.getFullPath()
       if not path:
           raise ValueError(f"{name} has no path set")

       # Determine columns from CONTENT_SIGNATURE_LIST ⭐
       content_flag = int(file_obj.contentFlag)
       columns = file_obj.CONTENT_SIGNATURE_LIST[content_flag - 1]

       # Build column_mapping (identity by default, or with renaming)
       if 'rename' in spec:
           column_mapping = {col: spec['rename'].get(col, col) for col in columns}
       else:
           column_mapping = {col: col for col in columns}

       # Build spec for merge_mtz_files
       input_specs.append({
           'path': path,
           'column_mapping': column_mapping
       })
   ```

3. **Call Utility**
   ```python
   output_path = self.workDirectory / f"{output_name}.mtz"
   return merge_mtz_files(input_specs, output_path, merge_strategy)
   ```

### Phase 3: Backward Compatibility (makeHklin)

Location: `core/CCP4PluginScript.py`

1. **Parse Old Syntax**
   ```python
   file_objects = []
   for item in miniMtzsIn:
       if isinstance(item, str):
           # Simple name - use object's contentFlag
           file_objects.append(item)
       else:
           # [name, explicit_flag] - override contentFlag
           name, explicit_flag = item
           columns = CONTENT_FLAG_COLUMNS[explicit_flag]
           file_objects.append({
               'name': name,
               'columns': columns
           })
   ```

2. **Call New API**
   ```python
   try:
       path = self.makeHklinGemmi(file_objects, output_name=hklin)
       # Success - return empty error report
       return CErrorReport()
   except Exception as e:
       error = CErrorReport()
       error.append(
           klass=self.__class__.__name__,
           code=200,
           details=f"Error merging MTZ files: {e}",
           name=hklin
       )
       return error
   ```

## Key gemmi Operations

```python
# Reading
mtz = gemmi.read_mtz_file(path)

# Creating output
out_mtz = gemmi.Mtz(with_base=True)  # Creates with H,K,L columns
out_mtz.spacegroup = mtz.spacegroup
out_mtz.set_cell_for_all(mtz.cell)

# Adding dataset (if needed)
out_mtz.add_dataset('merged_data')

# Copying columns (with their sigma columns)
src_col = mtz.column_with_label('F')
out_mtz.copy_column(-1, src_col, trailing_cols=['SIGF'])

# Renaming column
out_mtz.columns[-2].label = 'F_derivative'
out_mtz.columns[-1].label = 'SIGF_derivative'

# Update metadata and write
out_mtz.update_reso()
out_mtz.history = ['Merged by makeHklin']
out_mtz.write_to_file(str(output_path))
```

## Error Handling

```python
# Custom exception for MTZ operations
class MtzMergeError(Exception):
    """Errors during MTZ merging"""
    pass

# Error codes for CErrorReport
ERROR_MTZ_FILE_NOT_FOUND = 200
ERROR_MTZ_MERGE_FAILED = 201
ERROR_MTZ_INCOMPATIBLE_CELLS = 202
ERROR_MTZ_COLUMN_NOT_FOUND = 203
ERROR_MTZ_COLUMN_CONFLICT = 204
ERROR_CONTAINER_OBJECT_NOT_FOUND = 205
ERROR_UNKNOWN_CONTENT_FLAG = 206

# Usage in makeHklin wrapper
try:
    self.makeHklinGemmi(...)
except FileNotFoundError as e:
    error.append(
        klass=self.__class__.__name__,
        code=ERROR_MTZ_FILE_NOT_FOUND,
        details=str(e),
        name=hklin
    )
except MtzMergeError as e:
    error.append(
        klass=self.__class__.__name__,
        code=ERROR_MTZ_MERGE_FAILED,
        details=str(e),
        name=hklin
    )
```

## Testing Strategy

### Unit Tests (test_merge_mtz_files.py)

Test the low-level utility without CData dependencies:

```python
def test_merge_two_files_simple():
    """Test basic merge of two MTZ files."""
    result = merge_mtz_files(
        input_specs=[
            {'path': 'test1.mtz', 'columns': ['F', 'SIGF']},
            {'path': 'test2.mtz', 'columns': ['FreeR_flag']}
        ],
        output_path='merged.mtz'
    )
    assert result.exists()
    # Verify output has all columns

def test_merge_with_rename():
    """Test column renaming."""
    result = merge_mtz_files(
        input_specs=[
            {
                'path': 'test1.mtz',
                'columns': ['F', 'SIGF'],
                'rename': {'F': 'F_nat'}
            }
        ],
        output_path='out.mtz'
    )
    # Verify renamed column

def test_merge_strategy_error():
    """Test that conflict detection works."""
    with pytest.raises(ValueError, match="Column conflict"):
        merge_mtz_files(
            input_specs=[
                {'path': 'test1.mtz', 'columns': ['F']},
                {'path': 'test2.mtz', 'columns': ['F']}
            ],
            output_path='out.mtz',
            merge_strategy='error'
        )
```

### Integration Tests (test_makehklin.py)

Test CPluginScript methods with mock containers:

```python
def test_makehklin_gemmi_simple():
    """Test makeHklinGemmi with simple string list."""
    plugin = create_test_plugin_with_files()
    hklin = plugin.makeHklinGemmi(['HKLIN1', 'FREERFLAG'])
    assert hklin.exists()

def test_makehklin_backward_compat():
    """Test old makeHklin API."""
    plugin = create_test_plugin_with_files()
    error = plugin.makeHklin(['HKLIN1', 'FREERFLAG'])
    assert not error  # empty error report

def test_makehklin_with_content_flag_override():
    """Test [name, flag] syntax."""
    plugin = create_test_plugin_with_files()
    error = plugin.makeHklin([
        'HKLIN1',
        ['HKLIN2', CONTENT_FLAG_IOBS]
    ])
    assert not error
```

## Timeline Estimate

- **Phase 0** (Research contentFlags): ~~2-4 hours~~ ✅ **COMPLETED**
- Phase 1 (merge_mtz_files): 2-3 hours
- Phase 2 (makeHklinGemmi): 1-2 hours
- Phase 3 (makeHklin wrapper): 1 hour
- Testing: 2 hours
- Documentation: 30 min
- **Total**: **6-8 hours** (research completed)

## Critical Blockers ⚠️

**Status: MOSTLY RESOLVED**

1. ✅ Architecture clarified (DONE)
2. ✅ **CONTENT_FLAG constants extracted** from old ccp4i2 (DONE)
3. ✅ **Column name mappings documented** for each flag (DONE)
4. ✅ **CMiniMtzDataFile subclasses listed** and understood (DONE)
5. ⚠️ **Normalization status clarified** - do we need to implement it first? (OPTIONAL - can assume normalized)

## Questions Resolved ✅

1. ✅ **Container-level API** - Methods work with attribute names (e.g., 'HKLIN1')
2. ✅ **Path resolution** - Via CDataFile.getFullPath()
3. ✅ **contentFlag usage** - Determines which columns to copy
4. ✅ **Separation of concerns** - Low-level gemmi in CCP4Utils.py, integration in CPluginScript
5. ✅ **CONTENT_FLAG constants** - All extracted and documented (see above)
6. ✅ **Column mappings** - Complete mappings in CONTENT_SIGNATURE_LIST
7. ✅ **CMiniMtzDataFile subclasses** - 4 classes: CObsDataFile, CPhsDataFile, CMapCoeffsDataFile, CFreeRDataFile

## Questions Still Open ⚠️

1. ⚠️ **contentFlag access** - Is it `file_obj.contentFlag` or `file_obj.contentFlag.value` (CInt attribute)?
2. ⚠️ **Normalization** - Can we assume files are already normalized? Or implement check/normalizer?
3. ⚠️ **Dataset handling** - Single dataset vs. multiple datasets in output MTZ?

## Next Steps

1. ✅ ~~RESEARCH~~ - Extract CONTENT_FLAG info from old ccp4i2 (DONE)
2. ✅ ~~DOCUMENT~~ - Create complete reference table (DONE)
3. ✅ ~~Implement Phase 1~~ - merge_mtz_files utility (DONE - 10/10 tests passing)
4. ✅ ~~Refactor API~~ - Change to column_mapping for CData-agnosticism (DONE)
5. 🎯 **NEXT: Implement Phase 2** (makeHklinGemmi in CPluginScript)
6. ⏳ Implement Phase 3 (makeHklin wrapper in CPluginScript)
7. ⏳ Write integration tests
8. ⏳ Document and create examples

---

**Author**: Claude Code
**Date**: 2025-10-27
**Status**: Phase 1 Complete - Ready for Phase 2 (makeHklinGemmi)
**Last Updated**: 2025-10-27 (merge_mtz_files implemented, tested, and refactored for CData-agnosticism)
