# Complete CONTENT_FLAGS and SUBTYPES Metadata for CDataFile Descendants

**Date**: 2025-10-27
**Status**: ✅ COMPLETE - All metadata extracted and added to cdata.json

## Executive Summary

This document provides a comprehensive listing of all CDataFile descendant classes with CONTENT_FLAGS and/or SUBTYPES metadata. A systematic scan of the CCP4i2 codebase identified **6 classes across 3 source files** with these metadata types.

All metadata has been:
- ✅ Extracted from CCP4i2 source files
- ✅ Added to `migration/CData/cdata.json`
- ✅ Propagated to generated stub files in `core/cdata_stubs/`
- ✅ Verified with full test suite (151 tests passed)

## Complete Class Inventory

### Overview

| Class | Source File | Content Flags | Subtypes | Purpose |
|-------|-------------|---------------|----------|---------|
| CObsDataFile | CCP4XtalData.py | 4 | 3 | Observed diffraction data (MTZ) |
| CPhsDataFile | CCP4XtalData.py | 2 | 2 | Phase probability data (MTZ) |
| CMapCoeffsDataFile | CCP4XtalData.py | 1 | 3 | Map coefficients (MTZ) |
| CFreeRDataFile | CCP4XtalData.py | 1 | 0 | Free R flags (MTZ) |
| CPdbDataFile | CCP4ModelData.py | 2 | 5 | Protein coordinate files (PDB/mmCIF) |
| CCootHistoryDataFile | CCP4CootData.py | 0 | 2 | Coot modeling history files |

### Total CDataFile Descendants

The cdata.json contains **44 total CDataFile descendant classes**:

```
CAbInitioSeqsDataFile, CAminoAcidsDataFile, CCalcEDDataFile, CCifDataFile,
CComplementaryDataFile, CCootHistoryDataFile, CDataFileContent,
CDataSequenceDataFile, CDicDataFile, CDocDataFile, CFreeRDataFile,
CGFacDataFile, CLogDataFile, CMapCoeffsDataFile, CMapDataFile,
CMinimMtzDataFile, CMiniMtzDataFile, CObsDataFile, COptionsDataFile,
CPdbDataFile, CPhsDataFile, CPirDataFile, CPngDataFile, CProjectDataFile,
CReflectionDataFile, CSequenceChainsDataFile, CStructureDescriptionDataFile,
CTextDataFile, CXMLDataFile, CDiffractionProjectDataFile, CXrayDataFile,
CPictureDataFile, CGnuplotDataFile, CAmberParmsDataFile, CAmberPrepDataFile,
CFastaDataFile, CFileOrPathDataFile, CPictureFile, CPlotDataFile,
CSmallMolDescriptionDataFile, CBigStructureDataFile,
CBigStructurePictureDataFile, CExternalCommandDataFile, CAmberRstDataFile
```

Of these 44, only **6 have CONTENT_FLAGS or SUBTYPES metadata**.

## Detailed Class Metadata

### 1. CObsDataFile (Observed Diffraction Data)

**Source**: `/Users/nmemn/Developer/ccp4i2/core/CCP4XtalData.py`
**Purpose**: MTZ files containing experimental diffraction data

#### Subtypes

| Constant | Value | Description |
|----------|-------|-------------|
| SUBTYPE_OBSERVED | 1 | observed data |
| SUBTYPE_DERIVED | 2 | derived data |
| SUBTYPE_REFERENCE | 3 | reference data |

#### Content Flags

| Constant | Value | Annotation | Columns |
|----------|-------|------------|---------|
| CONTENT_FLAG_IPAIR | 1 | Anomalous Is | Iplus, SIGIplus, Iminus, SIGIminus |
| CONTENT_FLAG_FPAIR | 2 | Anomalous SFs | Fplus, SIGFplus, Fminus, SIGFminus |
| CONTENT_FLAG_IMEAN | 3 | Mean Is | I, SIGI |
| CONTENT_FLAG_FMEAN | 4 | Mean SFs | F, SIGF |

#### Generated Stub Constants

```python
class CObsDataFileStub(CMiniMtzDataFile):
    # Subtype constants
    SUBTYPE_OBSERVED = 1  # observed data
    SUBTYPE_DERIVED = 2  # derived data
    SUBTYPE_REFERENCE = 3  # reference data

    # Content flag constants
    CONTENT_FLAG_IPAIR = 1  # Anomalous Is
    CONTENT_FLAG_FPAIR = 2  # Anomalous SFs
    CONTENT_FLAG_IMEAN = 3  # Mean Is
    CONTENT_FLAG_FMEAN = 4  # Mean SFs

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = [
        'Anomalous Is',
        'Anomalous SFs',
        'Mean Is',
        'Mean SFs'
    ]

    # Content signatures - column names for each format
    CONTENT_SIGNATURE_LIST = [
        ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'],
        ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus'],
        ['I', 'SIGI'],
        ['F', 'SIGF']
    ]
```

#### Conversion Methods

```python
def as_IPAIR(self, work_directory: Optional[Any] = None) -> str
def as_FPAIR(self, work_directory: Optional[Any] = None) -> str
def as_IMEAN(self, work_directory: Optional[Any] = None) -> str
def as_FMEAN(self, work_directory: Optional[Any] = None) -> str
```

**Status**: Stub methods implemented, raise NotImplementedError

---

### 2. CPhsDataFile (Phase Probability Data)

**Source**: `/Users/nmemn/Developer/ccp4i2/core/CCP4XtalData.py`
**Purpose**: MTZ files containing phase probability information

#### Subtypes

| Constant | Value | Description |
|----------|-------|-------------|
| SUBTYPE_UNBIASED | 1 | unbiased data |
| SUBTYPE_BIASED | 2 | biased data |

#### Content Flags

| Constant | Value | Annotation | Columns |
|----------|-------|------------|---------|
| CONTENT_FLAG_HL | 1 | Hendrickson-Lattmann coeffs | HLA, HLB, HLC, HLD |
| CONTENT_FLAG_PHIFOM | 2 | Phi,FOM | PHI, FOM |

#### Generated Stub Constants

```python
class CPhsDataFileStub(CMiniMtzDataFile):
    # Subtype constants
    SUBTYPE_UNBIASED = 1  # unbiased data
    SUBTYPE_BIASED = 2  # biased data

    # Content flag constants
    CONTENT_FLAG_HL = 1  # Hendrickson-Lattmann coeffs
    CONTENT_FLAG_PHIFOM = 2  # Phi,FOM

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = [
        'Hendrickson-Lattmann coeffs',
        'Phi,FOM'
    ]

    # Content signatures
    CONTENT_SIGNATURE_LIST = [
        ['HLA', 'HLB', 'HLC', 'HLD'],
        ['PHI', 'FOM']
    ]
```

#### Conversion Methods

```python
def as_HL(self, work_directory: Optional[Any] = None) -> str
def as_PHIFOM(self, work_directory: Optional[Any] = None) -> str
```

**Status**: Stub methods implemented, raise NotImplementedError

---

### 3. CMapCoeffsDataFile (Map Coefficients)

**Source**: `/Users/nmemn/Developer/ccp4i2/core/CCP4XtalData.py`
**Purpose**: MTZ files containing map coefficients for electron density calculation

#### Subtypes

| Constant | Value | Description |
|----------|-------|-------------|
| SUBTYPE_NORMAL | 1 | normal map |
| SUBTYPE_DIFFERENCE | 2 | difference map |
| SUBTYPE_ANOM_DIFFERENCE | 3 | anomalous difference map |

#### Content Flags

| Constant | Value | Annotation | Columns |
|----------|-------|------------|---------|
| CONTENT_FLAG_FPHI | 1 | FPhi | F, PHI |

#### Generated Stub Constants

```python
class CMapCoeffsDataFileStub(CMiniMtzDataFile):
    # Subtype constants
    SUBTYPE_NORMAL = 1  # normal map
    SUBTYPE_DIFFERENCE = 2  # difference map
    SUBTYPE_ANOM_DIFFERENCE = 3  # anomalous difference map

    # Content flag constants
    CONTENT_FLAG_FPHI = 1  # FPhi

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = ['FPhi']

    # Content signatures
    CONTENT_SIGNATURE_LIST = [['F', 'PHI']]
```

#### Conversion Methods

```python
def as_FPHI(self, work_directory: Optional[Any] = None) -> str
```

**Note**: Only one content type, so this method just returns the file path (no conversion needed)

---

### 4. CFreeRDataFile (Free R Flags)

**Source**: `/Users/nmemn/Developer/ccp4i2/core/CCP4XtalData.py`
**Purpose**: MTZ files containing Free R flags for cross-validation

#### Subtypes

None defined.

#### Content Flags

| Constant | Value | Annotation | Columns |
|----------|-------|------------|---------|
| CONTENT_FLAG_FREER | 1 | FreeR | FREER |

#### Generated Stub Constants

```python
class CFreeRDataFileStub(CMiniMtzDataFile):
    # Content flag constants
    CONTENT_FLAG_FREER = 1  # FreeR

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = ['FreeR']

    # Content signatures
    CONTENT_SIGNATURE_LIST = [['FREER']]
```

#### Conversion Methods

None needed (only one content type).

---

### 5. CPdbDataFile (Protein Coordinate Files)

**Source**: `/Users/nmemn/Developer/ccp4i2/core/CCP4ModelData.py`
**Purpose**: Protein/nucleic acid coordinate files in PDB or mmCIF format

#### Subtypes

| Constant | Value | Description |
|----------|-------|-------------|
| SUBTYPE_UNKNOWN | 0 | unknown |
| SUBTYPE_MODEL | 1 | model |
| SUBTYPE_HOMOLOG | 2 | homolog |
| SUBTYPE_FRAGMENT | 3 | fragment |
| SUBTYPE_HEAVY_ATOMS | 4 | heavy atoms |

#### Content Flags

| Constant | Value | Annotation | Columns |
|----------|-------|------------|---------|
| CONTENT_FLAG_PDB | 1 | PDB format | N/A (not MTZ) |
| CONTENT_FLAG_MMCIF | 2 | mmCIF format | N/A (not MTZ) |

#### Generated Stub Constants

```python
class CPdbDataFileStub(CDataFile):
    # Subtype constants
    SUBTYPE_UNKNOWN = 0  # unknown
    SUBTYPE_MODEL = 1  # model
    SUBTYPE_HOMOLOG = 2  # homolog
    SUBTYPE_FRAGMENT = 3  # fragment
    SUBTYPE_HEAVY_ATOMS = 4  # heavy atoms

    # Content flag constants
    CONTENT_FLAG_PDB = 1  # PDB format
    CONTENT_FLAG_MMCIF = 2  # mmCIF format

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = ['PDB format', 'mmCIF format']
```

#### Conversion Methods

Not implemented. PDB↔mmCIF conversion would use a different pattern (likely gemmi-based utilities rather than the MTZ conversion system).

---

### 6. CCootHistoryDataFile (Coot Modeling History)

**Source**: `/Users/nmemn/Developer/ccp4i2/core/CCP4CootData.py`
**Purpose**: Coot modeling session history files (.scm Scheme scripts)

#### Subtypes

| Constant | Value | Description |
|----------|-------|-------------|
| SUBTYPE_INITIAL | 1 | Coot 0-state.scm |
| SUBTYPE_HISTORY | 2 | Coot history.scm |

#### Content Flags

None defined.

#### Generated Stub Constants

```python
class CCootHistoryDataFileStub(CDataFile):
    # Subtype constants
    SUBTYPE_INITIAL = 1  # Coot 0-state.scm
    SUBTYPE_HISTORY = 2  # Coot history.scm
```

#### Conversion Methods

Not applicable (Scheme script files, not binary data).

---

## Extraction Methodology

### Discovery Process

1. **Comprehensive Scan**: Used grep to search all Python files in `/Users/nmemn/Developer/ccp4i2/core/` for:
   - Classes with `CONTENT_FLAG_*` constants
   - Classes with `SUBTYPE_*` constants

2. **Source Files Examined**:
   ```
   /Users/nmemn/Developer/ccp4i2/core/CCP4XtalData.py    (MTZ data classes)
   /Users/nmemn/Developer/ccp4i2/core/CCP4ModelData.py   (Model/coordinate classes)
   /Users/nmemn/Developer/ccp4i2/core/CCP4CootData.py    (Coot integration classes)
   ```

3. **Validation**: Cross-referenced findings against cdata.json MRO to confirm all classes are CDataFile descendants

### Metadata Extraction Script

Created `$CCP4I2_ROOT/migration/CData/add_all_content_metadata.py` to:
- Define complete metadata for all 6 classes
- Update cdata.json with CONTENT_FLAGS and SUBTYPES
- Preserve existing metadata (CONTENTS, QUALIFIERS, etc.)
- Generate statistics and verification report

Script execution:
```bash
cd $CCP4I2_ROOT/migration/CData
python add_all_content_metadata.py
```

Output:
```
======================================================================
Adding CONTENT_FLAGS and SUBTYPES metadata to cdata.json
======================================================================

✓ CObsDataFile found in cdata.json
  ✓ Added 3 subtypes
  ✓ Added 4 content flags

✓ CPhsDataFile found in cdata.json
  ✓ Added 2 subtypes
  ✓ Added 2 content flags

✓ CMapCoeffsDataFile found in cdata.json
  ✓ Added 3 subtypes
  ✓ Added 1 content flags

✓ CFreeRDataFile found in cdata.json
  ✓ Added 1 content flags

✓ CPdbDataFile found in cdata.json
  ✓ Added 5 subtypes
  ✓ Added 2 content flags

✓ CCootHistoryDataFile found in cdata.json
  ✓ Added 2 subtypes

Classes processed: 6/6
Classes with SUBTYPES added: 5
Classes with CONTENT_FLAGS added: 5
```

---

## Files Modified

### 1. Metadata Source

**`$CCP4I2_ROOT/migration/CData/cdata.json`**

Added complete CONTENT_FLAGS and SUBTYPES metadata to 6 classes:

```json
{
  "classes": {
    "CObsDataFile": {
      "SUBTYPES": {
        "SUBTYPE_OBSERVED": {"value": 1, "description": "observed data"},
        "SUBTYPE_DERIVED": {"value": 2, "description": "derived data"},
        "SUBTYPE_REFERENCE": {"value": 3, "description": "reference data"}
      },
      "CONTENT_FLAGS": {
        "CONTENT_FLAG_IPAIR": {
          "value": 1,
          "annotation": "Anomalous Is",
          "columns": ["Iplus", "SIGIplus", "Iminus", "SIGIminus"]
        },
        ...
      }
    },
    ...
  }
}
```

### 2. Generated Stub Files

**`$CCP4I2_ROOT/core/cdata_stubs/CCP4XtalData.py`**
Updated with class constants for:
- CObsDataFileStub
- CPhsDataFileStub
- CMapCoeffsDataFileStub
- CFreeRDataFileStub

**`$CCP4I2_ROOT/core/cdata_stubs/CCP4ModelData.py`**
Updated with class constants for:
- CPdbDataFileStub

**`$CCP4I2_ROOT/core/cdata_stubs/CCP4CootData.py`**
Updated with class constants for:
- CCootHistoryDataFileStub

### 3. Implementation Files

**`$CCP4I2_ROOT/core/CCP4XtalData.py`**

Added conversion methods to CMiniMtzDataFile and subclasses:
- `CMiniMtzDataFile._get_conversion_output_path()` (base helper method)
- `CObsDataFile.as_IPAIR()`, `as_FPAIR()`, `as_IMEAN()`, `as_FMEAN()`
- `CPhsDataFile.as_HL()`, `as_PHIFOM()`
- `CMapCoeffsDataFile.as_FPHI()`

### 4. Plugin System

**`$CCP4I2_ROOT/core/CCP4PluginScript.py`**

Updated makeHklin to handle content type conversions:
- `_get_content_flag_name()` helper method
- Conversion detection logic
- Temporary file object creation
- Cleanup in finally block

### 5. Tests

**`$CCP4I2_ROOT/tests/test_cpluginscript_makehklin.py`**

Updated tests to verify conversion system:
- `test_override_content_flag_triggers_conversion()` - Expects NotImplementedError
- `test_no_conversion_when_flags_match()` - No conversion when flags match

---

## Usage Patterns

### Using Content Flags in Stubs

```python
from core.cdata_stubs.CCP4XtalData import CObsDataFileStub

# Access constants
flag = CObsDataFileStub.CONTENT_FLAG_FMEAN  # 4
subtype = CObsDataFileStub.SUBTYPE_OBSERVED  # 1
annotation = CObsDataFileStub.CONTENT_ANNOTATION[3]  # "Mean SFs"
columns = CObsDataFileStub.CONTENT_SIGNATURE_LIST[3]  # ['F', 'SIGF']
```

### Using in CPluginScript

```python
class MyTask(CPluginScript):
    def process(self):
        # File has F, SIGF (FMEAN format)
        # Request conversion to FPAIR (anomalous)
        error = self.makeHklin([
            ['HKLIN1', CObsDataFile.CONTENT_FLAG_FPAIR]
        ])

        # System automatically detects contentFlag mismatch and calls:
        # HKLIN1.as_FPAIR(self.workDirectory)
```

### Checking Content Type at Runtime

```python
obs_file = self.inputData.HKLIN1

# Get current content flag
current_flag = int(obs_file.contentFlag)  # e.g., 4

# Get human-readable annotation
if 1 <= current_flag <= len(obs_file.CONTENT_ANNOTATION):
    annotation = obs_file.CONTENT_ANNOTATION[current_flag - 1]
    print(f"File contains: {annotation}")  # "Mean SFs"

# Get expected columns
columns = obs_file.CONTENT_SIGNATURE_LIST[current_flag - 1]
print(f"Columns: {columns}")  # ['F', 'SIGF']
```

### Checking Subtypes

```python
pdb_file = self.inputData.XYZIN

# Check subtype
if hasattr(pdb_file, 'subtype'):
    if int(pdb_file.subtype) == CPdbDataFile.SUBTYPE_MODEL:
        print("This is a refined model")
    elif int(pdb_file.subtype) == CPdbDataFile.SUBTYPE_HOMOLOG:
        print("This is a homology model")
```

---

## Verification

### Test Results

Full test suite run: **151 tests passed, 26 skipped, 2 warnings**

```bash
$ python -m pytest $CCP4I2_ROOT/tests/ -v --tb=short
...
================= 151 passed, 26 skipped, 2 warnings in 1.56s ==================
```

Key tests passing:
- ✅ All 19 makeHklin integration tests
- ✅ Stub instantiation tests
- ✅ CPdbDataFile import and usage tests
- ✅ Full def.xml workflow tests
- ✅ XML serialization tests

### Stub Verification

Verified that generated stubs contain all expected constants:

```bash
$ grep -A 20 "class CObsDataFileStub" core/cdata_stubs/CCP4XtalData.py
class CObsDataFileStub(CMiniMtzDataFile):
    # Subtype constants
    SUBTYPE_OBSERVED = 1  # observed data
    SUBTYPE_DERIVED = 2  # derived data
    SUBTYPE_REFERENCE = 3  # reference data

    # Content flag constants
    CONTENT_FLAG_IPAIR = 1  # Anomalous Is
    CONTENT_FLAG_FPAIR = 2  # Anomalous SFs
    CONTENT_FLAG_IMEAN = 3  # Mean Is
    CONTENT_FLAG_FMEAN = 4  # Mean SFs

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = ['Anomalous Is', 'Anomalous SFs', 'Mean Is', 'Mean SFs']

    # Content signatures
    CONTENT_SIGNATURE_LIST = [['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'], ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus'], ['I', 'SIGI'], ['F', 'SIGF']]
```

---

## Next Steps

### 1. Implement MTZ Conversion Logic

Replace NotImplementedError stubs with actual gemmi-based conversions:

**Priority Conversions**:
- `CObsDataFile.as_FMEAN()` - Most commonly needed (I→F, IPAIR→FMEAN)
- `CObsDataFile.as_FPAIR()` - For anomalous analysis
- `CPhsDataFile.as_PHIFOM()` - HL→Phase/FOM conversion

**Implementation Pattern**:
```python
def as_FMEAN(self, work_directory: Optional[Any] = None) -> str:
    """Convert to Mean Structure Factors format."""
    output_path = self._get_conversion_output_path('FMEAN', work_directory)

    import gemmi
    mtz = gemmi.read_mtz_file(str(self.getFullPath()))

    # Conversion logic based on current contentFlag
    current_flag = int(self.contentFlag)

    if current_flag == self.CONTENT_FLAG_IMEAN:
        # I, SIGI → F, SIGF
        # F = sqrt(I)
        # SIGF = SIGI / (2 * sqrt(I))
        ...
    elif current_flag == self.CONTENT_FLAG_IPAIR:
        # Iplus, Iminus → F, SIGF
        # Calculate mean intensities first, then convert
        ...

    # Write output MTZ
    out_mtz.write_to_file(str(output_path))
    return str(output_path)
```

### 2. Add Unit Tests for Conversions

Create test suite for each conversion method:

```python
def test_fmean_to_fpair_conversion():
    """Test conversion from mean to anomalous structure factors."""
    obs = CObsDataFile()
    obs.setFullPath("/data/mean.mtz")
    obs.contentFlag = CInt(CObsDataFile.CONTENT_FLAG_FMEAN)

    converted = obs.as_FPAIR()

    # Verify output file exists
    assert Path(converted).exists()

    # Verify columns
    mtz = gemmi.read_mtz_file(converted)
    assert 'Fplus' in [col.label for col in mtz.columns]
    assert 'Fminus' in [col.label for col in mtz.columns]
```

### 3. Consider Non-MTZ Conversions

**CPdbDataFile**: PDB↔mmCIF conversion
- Could use gemmi's built-in conversion: `gemmi.read_structure()` / `gemmi.write_pdb()` / `gemmi.write_minimal_cif()`
- May want separate method pattern (not tied to MTZ conversion system)

**Example**:
```python
def as_MMCIF(self, work_directory: Optional[Any] = None) -> str:
    """Convert PDB to mmCIF format."""
    import gemmi

    structure = gemmi.read_structure(str(self.getFullPath()))
    output_path = self._calculate_output_path('mmcif', work_directory)

    gemmi.write_minimal_cif(structure, str(output_path))
    return str(output_path)
```

### 4. Performance Optimization

**Conversion Caching**: Avoid re-converting the same file multiple times
```python
class CMiniMtzDataFile:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._conversion_cache = {}

    def as_FMEAN(self, work_directory=None):
        cache_key = ('FMEAN', work_directory)
        if cache_key in self._conversion_cache:
            cached_path = self._conversion_cache[cache_key]
            if Path(cached_path).exists():
                return cached_path

        # Perform conversion
        output_path = self._convert_to_fmean(work_directory)
        self._conversion_cache[cache_key] = output_path
        return output_path
```

### 5. Documentation Enhancement

- Add crystallographic theory behind each conversion
- Document precision/accuracy implications
- Provide examples from real workflows
- Add troubleshooting guide for common conversion errors

---

## Related Documentation

- **MTZ_CONVERSION_SYSTEM.md** - Detailed architecture of the MTZ conversion system
- **MAKEHKLIN_IMPLEMENTATION_COMPLETE.md** - Overall makeHklin system
- **MERGE_MTZ_REFACTORING_SUMMARY.md** - API design rationale
- **MTZ_METADATA_EXTRACTION_SUMMARY.md** - Original CONTENT_FLAG extraction

---

## Appendix: CDataFile Descendants Without Metadata

The following 38 CDataFile descendants do **not** have CONTENT_FLAGS or SUBTYPES:

```
CAbInitioSeqsDataFile, CAminoAcidsDataFile, CCalcEDDataFile, CCifDataFile,
CComplementaryDataFile, CDataFileContent, CDataSequenceDataFile, CDicDataFile,
CDocDataFile, CGFacDataFile, CLogDataFile, CMapDataFile, CMinimMtzDataFile,
COptionsDataFile, CPirDataFile, CPngDataFile, CProjectDataFile,
CReflectionDataFile, CSequenceChainsDataFile, CStructureDescriptionDataFile,
CTextDataFile, CXMLDataFile, CDiffractionProjectDataFile, CXrayDataFile,
CPictureDataFile, CGnuplotDataFile, CAmberParmsDataFile, CAmberPrepDataFile,
CFastaDataFile, CFileOrPathDataFile, CPictureFile, CPlotDataFile,
CSmallMolDescriptionDataFile, CBigStructureDataFile,
CBigStructurePictureDataFile, CExternalCommandDataFile, CAmberRstDataFile,
CDataFile (base class)
```

These classes either:
- Have only one implicit format (no need for content flags)
- Don't support format variations (e.g., text files, logs)
- Use file extensions to distinguish types instead of internal flags

---

## Summary Statistics

- **Total CDataFile descendants**: 44 classes
- **Classes with metadata**: 6 classes (13.6%)
- **Total SUBTYPES defined**: 15 across 5 classes
- **Total CONTENT_FLAGS defined**: 10 across 5 classes
- **MTZ-related classes**: 4 of 6 (66.7%)
- **Non-MTZ classes**: 2 of 6 (33.3%)

**Completion Status**: ✅ 100%
- All metadata extracted
- All stubs regenerated
- All tests passing
- Documentation complete
