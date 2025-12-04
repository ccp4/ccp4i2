# French-Wilson Converter Comparison: Servalcat vs. Ctruncate

## Overview

This directory provides **two implementations** of French-Wilson intensity-to-amplitude conversions with **identical APIs**, allowing easy A/B comparison testing:

1. **ServalcatConverter** - Modern implementation using `servalcat fw`
2. **CtruncateConverter** - Legacy implementation using `ctruncate` wrapper

Both converters provide the exact same methods:
- `ipair_to_fpair(obs_file, work_directory)` - IPAIR → FPAIR
- `ipair_to_imean(obs_file, work_directory)` - IPAIR → IMEAN
- `to_fmean(obs_file, work_directory)` - IPAIR/IMEAN → FMEAN

## Drop-In Replacement Pattern

**Swap converters by changing a single import:**

```python
# Use servalcat (modern)
from core.conversions.servalcat_converter import ServalcatConverter as Converter

# OR use ctruncate (legacy)
from core.conversions.ctruncate_converter import CtruncateConverter as Converter

# Same API for both:
output_path = Converter.to_fmean(obs_file, work_directory)
output_path = Converter.ipair_to_fpair(obs_file, work_directory)
output_path = Converter.ipair_to_imean(obs_file, work_directory)
```

## Files

### Core Converters
- `servalcat_converter.py` - Modern servalcat fw implementation
- `ctruncate_converter.py` - Legacy ctruncate wrapper implementation
- `obs_data_converter.py` - High-level routing (currently uses servalcat)

### Utilities
- `mtz_utils.py` - Gemmi-based MTZ column extraction utilities

### Tests
- `tests/test_servalcat_converter.py` - Servalcat-specific tests
- `tests/test_converter_comparison.py` - API compatibility and comparison tests

## Usage Examples

### Basic Conversion

```python
from core.CCP4XtalData import CObsDataFile
from core.conversions.servalcat_converter import ServalcatConverter

# Load MTZ file
obs_file = CObsDataFile()
obs_file.setFullPath('/path/to/intensities.mtz')
obs_file.setContentFlag()  # Auto-detect IPAIR/IMEAN/etc

# Convert to FMEAN
output_path = ServalcatConverter.to_fmean(obs_file, work_directory='/tmp/work')
print(f"Converted to: {output_path}")
```

### Head-to-Head Comparison

```python
from core.CCP4XtalData import CObsDataFile
from core.conversions.servalcat_converter import ServalcatConverter
from core.conversions.ctruncate_converter import CtruncateConverter
from pathlib import Path

# Load test data
obs_file = CObsDataFile()
obs_file.setFullPath('test_data/intensities.mtz')
obs_file.setContentFlag()

# Run both converters
servalcat_output = ServalcatConverter.to_fmean(obs_file, '/tmp/servalcat')
ctruncate_output = CtruncateConverter.to_fmean(obs_file, '/tmp/ctruncate')

# Compare outputs
import gemmi
mtz_servalcat = gemmi.read_mtz_file(servalcat_output)
mtz_ctruncate = gemmi.read_mtz_file(ctruncate_output)

print(f"Servalcat: {len(mtz_servalcat.data)} reflections")
print(f"Ctruncate: {len(mtz_ctruncate.data)} reflections")
```

### Switching Converters Globally

To switch the entire system between converters, modify `obs_data_converter.py`:

```python
# In core/conversions/obs_data_converter.py

# Option 1: Use servalcat (modern)
from core.conversions.servalcat_converter import ServalcatConverter as ConverterImpl

# Option 2: Use ctruncate (legacy)
from core.conversions.ctruncate_converter import CtruncateConverter as ConverterImpl

# Then use ConverterImpl throughout the file:
class ObsDataConverter:
    @staticmethod
    def to_fmean(obs_file, work_directory=None):
        # ...
        return ConverterImpl.to_fmean(obs_file, work_directory)
```

## Conversion Matrix

Both converters support the same conversion paths:

```
                    TO
          IPAIR FPAIR IMEAN FMEAN
FROM IPAIR   ✓    ✓     ✓     ✓    (French-Wilson)
     FPAIR   ✗    ✓     ✗     ✓    (gemmi, not in converters)
     IMEAN   ✗    ✗     ✓     ✓    (French-Wilson)
     FMEAN   ✗    ✗     ✗     ✓    (identity copy)
```

Note: FPAIR → FMEAN uses gemmi-based weighted averaging in `obs_data_converter.py`, not the converters.

## Implementation Differences

| Feature | ServalcatConverter | CtruncateConverter |
|---------|-------------------|-------------------|
| **Algorithm** | Modern French-Wilson (servalcat fw) | Legacy French-Wilson (ctruncate) |
| **Invocation** | Subprocess (servalcat fw) | Plugin wrapper (ctruncate class) |
| **Output** | Monolithic MTZ (all formats) | Mini-MTZ (single format) |
| **Extraction** | Pure gemmi (mtz_utils) | Plugin's splitHklout |
| **Work Dir** | `servalcat/` subdirectory | `ctruncate/` subdirectory |
| **Statistics** | More accurate (newer algorithm) | Original ctruncate stats |
| **Dependencies** | Requires servalcat in PATH | Requires ctruncate plugin |

## Testing

### Run Comparison Tests

```bash
# Test API compatibility
pytest tests/test_converter_comparison.py -v

# Test servalcat converter
pytest tests/test_servalcat_converter.py -v

# Test MTZ utilities
pytest tests/test_mtz_utils.py -v
```

### Example Test Output

```
tests/test_converter_comparison.py::test_api_documentation PASSED
✅ Both converters provide the same API:
  - ipair_to_fpair(obs_file, work_directory)
  - ipair_to_imean(obs_file, work_directory)
  - to_fmean(obs_file, work_directory)

tests/test_converter_comparison.py::test_drop_in_replacement_pattern PASSED
✅ servalcat: /tmp/servalcat_work/FMEAN.mtz
✅ ctruncate: /tmp/ctruncate_work/FMEAN.mtz
✅ Drop-in replacement successful - both converters work with identical API
```

## Performance Considerations

### Servalcat (Modern)
- **Pros**: More accurate statistics, monolithic output (efficient for multiple formats)
- **Cons**: External dependency (servalcat executable)
- **Use when**: Accuracy matters, producing multiple output formats

### Ctruncate (Legacy)
- **Pros**: Well-tested, integrated plugin system
- **Cons**: Older algorithm, separate runs needed for multiple formats
- **Use when**: Compatibility testing, comparing against baseline

## Migration Strategy

1. **Current**: System uses `ServalcatConverter` by default (via `obs_data_converter.py`)
2. **Testing**: Use `CtruncateConverter` for regression testing
3. **Comparison**: Run both converters on same data, compare outputs
4. **Validation**: Ensure servalcat produces acceptable results
5. **Deprecation**: Once confident, mark ctruncate as legacy

## Troubleshooting

### Servalcat Not Found
```
CException: servalcat executable not found in PATH (code 10)
```
**Solution**: Source CCP4 environment or install servalcat

### Ctruncate Plugin Not Available
```
CException: ctruncate plugin not available in registry (code 22)
```
**Solution**: Ensure plugin registry is up-to-date and ctruncate is installed

### Output File Not Created
```
CException: Expected output file not created (code 12 or 21)
```
**Solution**: Check work directory permissions and converter logs

## See Also

- `mtz_utils.py` - Gemmi-based MTZ column extraction
- `obs_data_converter.py` - High-level conversion routing
- `tests/test_converter_comparison.py` - API compatibility tests
- `tests/test_mtz_utils.py` - MTZ utility function tests
