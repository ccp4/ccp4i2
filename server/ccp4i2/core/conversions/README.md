# Crystallographic Data Format Converters

This directory contains specialized converter modules for transforming crystallographic data between different formats. The converter architecture follows a consistent pattern that separates conversion logic from data model classes.

## Architecture Pattern

### Thin Wrapper Design

Each converter module uses the **thin wrapper pattern**:

1. **Data Model Classes** (e.g., `CObsDataFile`, `CMapCoeffsDataFile`)
   - Contain minimal wrapper methods like `as_FMEAN()`, `as_FPHI()`
   - Delegate to converter modules for implementation
   - Maintain intuitive API: `obs_file.as_FMEAN()`

2. **Converter Modules** (e.g., `ObsDataConverter`, `PhaseDataConverter`)
   - Contain all conversion logic
   - Static methods take file instance as first parameter
   - Can be used independently: `ObsDataConverter.to_fmean(obs_file)`

### Benefits

✅ **Separation of Concerns**: Data structure vs conversion algorithms
✅ **Maintainability**: All conversion logic in one place per data type
✅ **Testability**: Can test converters with mock file objects
✅ **API Stability**: User-facing methods remain unchanged
✅ **Documentation**: Comprehensive docs in converter modules
✅ **Extensibility**: Clear pattern for adding new converters

## Converter Modules

### 1. ObsDataConverter (`obs_data_converter.py`)

**Status**: ✅ Fully Implemented

**Purpose**: MTZ observation data format conversions

**Formats**:
- **IPAIR (1)**: Anomalous intensities (I+, SIGI+, I-, SIGI-)
- **FPAIR (2)**: Anomalous structure factors (F+, SIGF+, F-, SIGF-)
- **IMEAN (3)**: Mean intensities (I, SIGI)
- **FMEAN (4)**: Mean structure factors (F, SIGF)

**Conversion Matrix**:
```
                TO
      IPAIR FPAIR IMEAN FMEAN
FROM
IPAIR   ✓     ✓     ✓     ✓
FPAIR   ✗     ✓     ✗     ✓
IMEAN   ✗     ✗     ✓     ✓
FMEAN   ✗     ✗     ✗     ✓
```

**Implementation Methods**:
- `IPAIR → FPAIR/IMEAN/FMEAN`: Uses ctruncate (French-Wilson conversion)
- `IMEAN → FMEAN`: Uses ctruncate (French-Wilson conversion)
- `FPAIR → FMEAN`: Uses gemmi + numpy (inverse-variance weighted mean)

**Usage**:
```python
from ccp4i2.core.CCP4XtalData import CObsDataFile

# Via thin wrapper (recommended)
obs_file = CObsDataFile("data.mtz")
fmean_path = obs_file.as_FMEAN(work_directory="./output")

# Or directly
from ccp4i2.core.conversions import ObsDataConverter
fmean_path = ObsDataConverter.to_fmean(obs_file, work_directory="./output")
```

**Tests**: `tests/test_obs_conversions.py` - All passing (4/4)

---

### 2. PhaseDataConverter (`phase_data_converter.py`)

**Status**: ✅ Fully Implemented

**Purpose**: Phase data format conversions using gemmi

**Formats**:
- **HL (1)**: Hendrickson-Lattman coefficients (HLA, HLB, HLC, HLD)
- **PHIFOM (2)**: Phase + Figure of Merit (PHI, FOM)
- **FPHI (1)**: Structure factors + Phase (F, PHI)

**Conversion Matrix** (CPhsDataFile):
```
           TO
       HL  PHIFOM
FROM
HL      ✓     ✓
PHIFOM  ✓     ✓
```

**Implementation Methods**:
- `HL → PHIFOM`: Uses CCP4 chltofom plugin for accurate conversion
  - Validated reference implementation from CCP4 suite
  - Calculates best phase and FOM from HL coefficient probability distribution
  - Post-processes with gemmi to extract only converted columns (ensures correct contentFlag)
- `PHIFOM → HL`: Uses CCP4 chltofom plugin for accurate conversion
  - Converts PHI/FOM to HL coefficient representation
  - Post-processes with gemmi to extract only converted columns (ensures correct contentFlag)
- `FPHI` conversions: Pass-through (CMapCoeffsDataFile has only one format)

**Note**: The chltofom plugin outputs both original and converted columns by default.
We post-process with gemmi to extract only the desired columns so that `setContentFlag()`
correctly identifies the output file's format. The module also includes experimental
gemmi-based implementations using numerical methods, but these are not used in
production due to FOM calculation differences (correlation 0.08 vs chltofom).
The experimental code is preserved for future investigation and refinement.

**Usage**:
```python
from ccp4i2.core.CCP4XtalData import CPhsDataFile

# HL → PHIFOM conversion
hl_file = CPhsDataFile("phases.mtz")
phifom_path = hl_file.as_PHIFOM(work_directory="./output")

# PHIFOM → HL conversion
phifom_file = CPhsDataFile("phases.mtz")
hl_path = phifom_file.as_HL(work_directory="./output")
```

**Tests**: `tests/test_phase_conversions.py` - All passing (4/4)
- HL → PHIFOM conversion
- PHIFOM → HL conversion
- Round-trip conversion (HL → PHIFOM → HL)
- Numerical accuracy validation

**References**:
- Read, R.J. (1986). *Acta Cryst.* **A42**, 140-149.
- Hendrickson & Lattman (1970). *Acta Cryst.* **B26**, 136-143.
- gemmi documentation: `docs/hkl.rst` (HL coefficient handling)

---

### 3. ModelConverter (`model_converter.py`)

**Status**: 🚧 Stub Implementation (TODO)

**Purpose**: Macromolecular model format conversions

**Formats**:
- **PDB**: Protein Data Bank format (legacy column-based)
- **mmCIF**: Macromolecular Crystallographic Information File (modern structured)

**Conversion Matrix**:
```
           TO
       PDB  mmCIF
FROM
PDB     ✓     ✓
mmCIF   ✓     ✓
```

**Implementation Approach**:
- Use **gemmi** library for robust format conversion
- `gemmi.read_structure()` auto-detects format
- `structure.write_pdb()` for PDB output
- `structure.make_mmcif_document().write_file()` for mmCIF output

**Future Implementation**:
```python
import gemmi

def to_pdb(model_file, work_directory=None):
    input_path = model_file.getFullPath()
    output_path = model_file._get_conversion_output_path('pdb',
                   target_extension='.pdb', work_directory=work_directory)

    structure = gemmi.read_structure(input_path)
    structure.write_pdb(output_path)

    return output_path
```

**Conversion Considerations**:
- **PDB limitations**: 9999 atom limit, single-character chain IDs
- **mmCIF → PDB**: May require truncating long names, chain ID mapping
- **PDB → mmCIF**: Straightforward, all information preserved
- **Metadata**: Preserve crystal symmetry, B-factors, occupancies
- **Validation**: Consider adding gemmi-based validation methods

**Usage** (when implemented):
```python
from ccp4i2.core.CCP4ModelData import CPdbDataFile

pdb_file = CPdbDataFile("model.pdb")
cif_path = pdb_file.as_mmcif(work_directory="./output")
```

---

## Adding New Converters

To add a new converter module:

### 1. Create Converter Module

Create `core/conversions/new_converter.py`:

```python
"""
Converter for [data type] format transformations.

[Description of formats and conversions]
"""

from typing import Optional, Any

class NewConverter:
    """Static converter class for [data type] transformations."""

    @staticmethod
    def to_format_a(data_file, work_directory: Optional[Any] = None) -> str:
        """Convert to format A."""
        # Implementation here
        pass

    @staticmethod
    def to_format_b(data_file, work_directory: Optional[Any] = None) -> str:
        """Convert to format B."""
        # Implementation here
        pass
```

### 2. Update `__init__.py`

Add to `core/conversions/__init__.py`:

```python
from .new_converter import NewConverter

__all__ = [
    'ObsDataConverter',
    'PhaseDataConverter',
    'ModelConverter',
    'NewConverter',  # Add here
]
```

### 3. Add Thin Wrappers

Add wrapper methods to data file class:

```python
class CNewDataFile(CNewDataFileStub):

    def as_format_a(self, work_directory=None):
        """Convert to format A."""
        from ccp4i2.core.conversions import NewConverter
        return NewConverter.to_format_a(self, work_directory=work_directory)
```

### 4. Write Tests

Create `tests/test_new_conversions.py`:

```python
def test_format_a_to_format_b(tmp_path):
    file = CNewDataFile("input.ext")
    output = file.as_format_b(work_directory=str(tmp_path))
    assert Path(output).exists()
```

## Testing

Run converter tests:

```bash
# All conversion tests
pytest tests/test_obs_conversions.py -v

# Specific converter
pytest tests/test_obs_conversions.py::test_ipair_to_fmean -v

# With output
pytest tests/test_obs_conversions.py -v -s
```

## Dependencies

- **ctruncate**: CCP4 program for French-Wilson conversion (IPAIR/IMEAN → FMEAN/FPAIR)
- **chltofom**: CCP4 program for phase data conversions (HL ↔ PHIFOM)
- **gemmi**: Python library for crystallographic file I/O
  - MTZ file handling
  - PDB/mmCIF conversions
  - Crystallographic calculations
- **numpy**: Numerical computations (weighted averaging, experimental phase calculations)

## Future Work

### High Priority
1. **Implement PDB ↔ mmCIF conversions** (ModelConverter)
   - Use gemmi library (already in dependencies)
   - Add validation methods
   - Handle edge cases (long names, chain IDs)

### Medium Priority
2. **Add MAP file conversions**
   - Create MapConverter module
   - CCP4 MAP ↔ other map formats
   - May use gemmi or CCP4 programs

3. **Expand test coverage for all converters**
   - Add edge case tests (missing data, unusual formats)
   - Test error handling
   - Performance benchmarks for large files

4. **Add conversion validation**
   - Verify output correctness
   - Compare statistics before/after conversion
   - Warn about information loss

### Low Priority
5. **Optimization**
   - Cache converted files
   - Parallel conversions
   - Progress reporting for large files

## References

- **CCP4**: https://www.ccp4.ac.uk/
- **gemmi**: https://gemmi.readthedocs.io/
- **French & Wilson** (1978). *Acta Cryst.* **A34**, 517-525. (French-Wilson method)
- **Read** (1986). *Acta Cryst.* **A42**, 140-149. (Phase probability distributions)

---

*Last updated: 2025-01-XX*
*Maintainer: CCP4i2 project*
