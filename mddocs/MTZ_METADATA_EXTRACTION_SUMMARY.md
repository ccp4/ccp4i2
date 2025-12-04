# MTZ Metadata Extraction Summary

**Date**: 2025-10-27
**Task**: Extract MTZ-specific metadata from old CCP4i2 and add to cdata.json

## What Was Extracted

From `/Users/nmemn/Developer/ccp4i2/core/CCP4XtalData.py`, extracted class-level metadata for all CMiniMtzDataFile subclasses:

### 1. CObsDataFile (Observed Crystallographic Data)

**Content Flags:**
- `CONTENT_FLAG_IPAIR = 1` - Anomalous Is → Columns: `[Iplus, SIGIplus, Iminus, SIGIminus]`
- `CONTENT_FLAG_FPAIR = 2` - Anomalous SFs → Columns: `[Fplus, SIGFplus, Fminus, SIGFminus]`
- `CONTENT_FLAG_IMEAN = 3` - Mean Is → Columns: `[I, SIGI]`
- `CONTENT_FLAG_FMEAN = 4` - Mean SFs → Columns: `[F, SIGF]`

**Subtypes:**
- `SUBTYPE_OBSERVED = 1` - observed data
- `SUBTYPE_DERIVED = 2` - derived data
- `SUBTYPE_REFERENCE = 3` - reference data

### 2. CPhsDataFile (Phase Data)

**Content Flags:**
- `CONTENT_FLAG_HL = 1` - Hendrickson-Lattmann coeffs → Columns: `[HLA, HLB, HLC, HLD]`
- `CONTENT_FLAG_PHIFOM = 2` - Phi,FOM → Columns: `[PHI, FOM]`

**Subtypes:**
- `SUBTYPE_UNBIASED = 1` - unbiased data
- `SUBTYPE_BIASED = 2` - biased data

### 3. CMapCoeffsDataFile (Map Coefficients)

**Content Flags:**
- `CONTENT_FLAG_FPHI = 1` - FPhi → Columns: `[F, PHI]`

**Subtypes:**
- `SUBTYPE_NORMAL = 1` - normal map
- `SUBTYPE_DIFFERENCE = 2` - difference map
- `SUBTYPE_ANOM_DIFFERENCE = 3` - anomalous difference map

### 4. CFreeRDataFile (Free R Flags)

**Content Flags:**
- `CONTENT_FLAG_FREER = 1` - FreeR → Columns: `[FREER]`

**Subtypes:** (None)

## Files Modified

### 1. `/Users/nmemn/Developer/cdata-codegen/migration/CData/cdata.json`

Added two new top-level keys to each CMiniMtzDataFile subclass:

```json
{
  "CObsDataFile": {
    "CONTENT_FLAGS": {
      "CONTENT_FLAG_IPAIR": {
        "value": 1,
        "annotation": "Anomalous Is",
        "columns": ["Iplus", "SIGIplus", "Iminus", "SIGIminus"]
      },
      ...
    },
    "SUBTYPES": {
      "SUBTYPE_OBSERVED": {
        "value": 1,
        "description": "observed data"
      },
      ...
    }
  }
}
```

### 2. `/Users/nmemn/Developer/cdata-codegen/migration/CData/add_mtz_metadata.py`

Created script to automate metadata addition. Can be re-run if needed.

**Usage:**
```bash
python3 migration/CData/add_mtz_metadata.py
```

### 3. `/Users/nmemn/Developer/cdata-codegen/MAKEHKLIN_REFACTOR_PLAN.md`

Updated refactoring plan with:
- Complete content flag definitions
- Column name mappings
- Marked Phase 0 (research) as completed
- Updated blockers and questions sections
- Status changed to "Ready for Implementation"

## How This Metadata Will Be Used

### In Code Generation

The stub generator can now create class constants:

```python
@cdata_class(...)
class CObsDataFile(CMiniMtzDataFile):
    # Content flags
    CONTENT_FLAG_IPAIR = 1
    CONTENT_FLAG_FPAIR = 2
    CONTENT_FLAG_IMEAN = 3
    CONTENT_FLAG_FMEAN = 4

    # Subtypes
    SUBTYPE_OBSERVED = 1
    SUBTYPE_DERIVED = 2
    SUBTYPE_REFERENCE = 3
```

### In makeHklin Implementation

The metadata enables automatic column lookup:

```python
def makeHklinGemmi(self, file_objects, ...):
    for spec in file_objects:
        # Lookup file object
        file_obj = self.inputData.HKLIN1  # CObsDataFile instance

        # Get its content flag
        content_flag = file_obj.contentFlag  # e.g., 4 (FMEAN)

        # Lookup columns from metadata
        class_name = file_obj.__class__.__name__  # "CObsDataFile"
        columns = CONTENT_FLAG_COLUMNS[class_name][content_flag]
        # columns = ['F', 'SIGF']

        # Extract from MTZ and merge
        ...
```

## Key Insights

1. **Four distinct CMiniMtzDataFile subclasses** - Each handles different types of crystallographic data
2. **contentFlag is integer** - Stored as CInt attribute, values start at 1
3. **Column names are standardized** - Each contentFlag implies specific, normalized column names
4. **Subtypes add semantic meaning** - Distinguish observed vs. derived data, biased vs. unbiased phases, etc.
5. **Already in cdata.json** - Classes were already defined, just missing this metadata

## What's Next

With metadata extraction complete, the makeHklin implementation can proceed:

1. **Phase 1**: Implement `merge_mtz_files()` utility in `core/CCP4Utils.py`
   - Pure gemmi operations, no CData dependencies
   - Takes file paths and column names

2. **Phase 2**: Implement `makeHklinGemmi()` in `CPluginScript`
   - Container-level API
   - Uses metadata to map contentFlag → columns

3. **Phase 3**: Implement `makeHklin()` backward-compatible wrapper
   - Translates old syntax to new API

## Verification

```bash
# Verify metadata was added
python3 << 'EOF'
import json
with open('migration/CData/cdata.json', 'r') as f:
    data = json.load(f)
    obs = data['classes']['CObsDataFile']
    print(f"CONTENT_FLAGS: {len(obs['CONTENT_FLAGS'])} defined")
    print(f"SUBTYPES: {len(obs['SUBTYPES'])} defined")
    print(f"FMEAN columns: {obs['CONTENT_FLAGS']['CONTENT_FLAG_FMEAN']['columns']}")
EOF
```

**Output:**
```
CONTENT_FLAGS: 4 defined
SUBTYPES: 3 defined
FMEAN columns: ['F', 'SIGF']
```

## References

- **Source**: `/Users/nmemn/Developer/ccp4i2/core/CCP4XtalData.py` (lines 3081-3472)
- **Updated**: `/Users/nmemn/Developer/cdata-codegen/migration/CData/cdata.json`
- **Plan**: `/Users/nmemn/Developer/cdata-codegen/MAKEHKLIN_REFACTOR_PLAN.md`
- **Script**: `/Users/nmemn/Developer/cdata-codegen/migration/CData/add_mtz_metadata.py`
