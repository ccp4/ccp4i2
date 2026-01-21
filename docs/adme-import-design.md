# ADME Data Import System - Design Document

## Overview

This document describes a system for importing ADME (Absorption, Distribution, Metabolism, Excretion) assay data from external vendors into the CCP4i2 Assays data model. The initial implementation targets data from NCU (outsourced CRO), but the architecture is designed to be extensible to other vendors.

## Problem Statement

Outsourced ADME studies return results in Excel files with vendor-specific formats. Each file contains:
- Multiple sheets (Signature, Summary, Materials, Study Design, Bioanalytical Method, raw data)
- Pre-calculated results on a "Summary" sheet
- Multiple tables per sheet with different structures
- Data for multiple compounds, potentially with multiple species or assay conditions per compound

We need to:
1. Parse these disparate formats reliably
2. Map compound identifiers to existing Compound records in the registry
3. Store results in the `AnalysisResult.results` JSON field
4. Maintain traceability to source files

## Vendor: NCU

### File Naming Convention
```
ADME-NCU-{AssayType}-{YYYYMMDD}.xlsx
```

### Supported Assay Types

| Code | Full Name | Description |
|------|-----------|-------------|
| LM | Liver Microsome Stability | Metabolic stability in liver microsomes (human/mouse) |
| BS | Blood/Serum Stability | Stability in blood or serum |
| GSH | GSH Stability | Glutathione stability (reactive metabolite screening) |
| Caco-2 | Caco-2 Permeability | Intestinal permeability using Caco-2 cell monolayers |

### Common Sheet Structure

All NCU files share this sheet structure:
- `Signature` - Sign-off page
- `Summary` - **Key data extraction target**
- `Materials` - Reagents and materials used
- `Study Design` - Protocol parameters
- `Bioanalytical Method` - LC-MS/MS method details
- Additional sheets with raw data (assay-specific)

### Summary Sheet Structure

The Summary sheet follows a consistent pattern:
- Row 0: "Data Summary"
- Row 1: "Table 1. {description}"
- Row 2: Column headers
- Row 3+: Data rows
- Tables separated by blank rows or "Note:" rows

---

## Assay-Specific Schemas

### 1. Liver Microsome Stability (LM)

**Source file example**: `ADME-NCU-LM-20231218.xlsx`

**Additional sheets**: `HLM` (Human), `MLM` (Mouse) with raw analytical data

#### Table 1: Summary Results
| Column | Field Name | Type | Unit |
|--------|------------|------|------|
| Compound ID | `compound_id` | string | - |
| Species | `species` | enum | Human/Mouse |
| in vitro t1/2 (min) | `t1_2_min` | float | min |
| in vitro CLint (μL/min/mg) | `clint_ul_min_mg` | float | μL/min/mg protein |
| Scale-up CLint (mL/min/kg) | `clint_scaled_ml_min_kg` | float | mL/min/kg |
| Predicted Hepatic CLH (mL/min/kg) | `clh_predicted_ml_min_kg` | float | mL/min/kg |
| Hepatic Extraction Ratio (ER) | `extraction_ratio` | float | 0-1 |

**Note**: Compound ID spans multiple rows when multiple species tested. Use forward-fill.

#### Table 2: Time Course Data
| Column | Field Name | Type |
|--------|------------|------|
| Compound ID | `compound_id` | string |
| Species | `species` | enum |
| Assay Format | `assay_format` | enum (+Cofactors/-Cofactors) |
| 0.5 min | `remaining_pct[0]` | float |
| 5 min | `remaining_pct[1]` | float |
| 15 min | `remaining_pct[2]` | float |
| 30 min | `remaining_pct[3]` | float |
| 60 min | `remaining_pct[4]` | float |

**Time points**: [0.5, 5, 15, 30, 60] minutes

**Special values**:
- `BLOD` = Below Limit of Detection -> store as `null` with flag
- `-` = Not measured -> store as `null`

#### JSON Output Schema
```json
{
    "assay_type": "liver_microsome_stability",
    "KPI": "t1_2_min",
    "species": "Human",
    "t1_2_min": 53.46,
    "clint_ul_min_mg": 25.95,
    "clint_scaled_ml_min_kg": 26.67,
    "clh_predicted_ml_min_kg": 11.75,
    "extraction_ratio": 0.56,
    "time_course": {
        "time_points_min": [0.5, 5, 15, 30, 60],
        "with_cofactors": {
            "remaining_pct": [100, 89.48, 79.02, 68.12, 44.76]
        },
        "without_cofactors": {
            "remaining_pct": [100, null, null, null, 104.85]
        }
    },
    "flags": [],
    "source": {
        "vendor": "NCU",
        "file": "ADME-NCU-LM-20231218.xlsx",
        "date": "2023-12-18"
    }
}
```

---

### 2. Blood/Serum Stability (BS)

**Source file example**: `ADME-NCU-BS-20240906.xlsx`

#### Table 1: Stability Results
| Column | Field Name | Type | Unit |
|--------|------------|------|------|
| Compound ID | `compound_id` | string | - |
| Species | `species` | enum | Mouse/Human/etc |
| 0 min | `remaining_pct[0]` | float | % |
| 15 min | `remaining_pct[1]` | float | % |
| 30 min | `remaining_pct[2]` | float | % |
| 60 min | `remaining_pct[3]` | float | % |
| 120 min | `remaining_pct[4]` | float | % |
| t1/2 (min) | `t1_2_min` | float/string | min |

**Time points**: [0, 15, 30, 60, 120] minutes

**Special values**:
- `∞` = Infinite half-life (stable) -> store as `null` with flag `"stable"`

#### JSON Output Schema
```json
{
    "assay_type": "blood_serum_stability",
    "KPI": "t1_2_min",
    "species": "Mouse",
    "t1_2_min": 49.93,
    "time_course": {
        "time_points_min": [0, 15, 30, 60, 120],
        "remaining_pct": [100, 82.53, 72.29, 47.74, 19.02]
    },
    "flags": [],
    "source": {
        "vendor": "NCU",
        "file": "ADME-NCU-BS-20240906.xlsx",
        "date": "2024-09-06"
    }
}
```

---

### 3. GSH Stability (Glutathione)

**Source file example**: `ADME-NCU-GSH Stability-20241106.xlsx`

#### Table 1: GSH Stability Results
| Column | Field Name | Type | Unit |
|--------|------------|------|------|
| Compound ID | `compound_id` | string | - |
| Assay Format | `assay_format` | enum | +GSH/-GSH |
| 0 min | `remaining_pct[0]` | float | % |
| 15 min | `remaining_pct[1]` | float | % |
| 30 min | `remaining_pct[2]` | float | % |
| 60 min | `remaining_pct[3]` | float | % |
| 120 min | `remaining_pct[4]` | float | % |
| t1/2 (min) | `t1_2_min` | float/string | min |

**Time points**: [0, 15, 30, 60, 120] minutes

**Special values**:
- `> 686.01` = Half-life exceeds measurement window -> store numeric part with flag
- `-` = Not measured (for -GSH time points) -> store as `null`

#### JSON Output Schema
```json
{
    "assay_type": "gsh_stability",
    "KPI": "t1_2_min",
    "t1_2_min": 55.76,
    "t1_2_min_no_gsh": null,
    "time_course": {
        "time_points_min": [0, 15, 30, 60, 120],
        "with_gsh": {
            "remaining_pct": [100, 80.25, 63.50, 41.89, 22.40]
        },
        "without_gsh": {
            "remaining_pct": [100, null, null, null, 99.45]
        }
    },
    "gsh_reactive": true,
    "flags": [],
    "source": {
        "vendor": "NCU",
        "file": "ADME-NCU-GSH Stability-20241106.xlsx",
        "date": "2024-11-06"
    }
}
```

**Derived field**: `gsh_reactive` = true if t1/2 with GSH is significantly lower than without GSH

---

### 4. Caco-2 Permeability

**Source file example**: `ADME-NCU-Caco-2 Permeability-20231219.xlsx`

#### Table 1: Permeability Results
| Column | Field Name | Type | Unit |
|--------|------------|------|------|
| Compound ID | `compound_id` | string | - |
| Papp (A-B) | `papp_ab` | float | 10^-6 cm/s |
| Papp (B-A) | `papp_ba` | float | 10^-6 cm/s |
| Efflux Ratio | `efflux_ratio` | float | - |
| Recovery (%) AP-BL | `recovery_ab_pct` | float | % |
| Recovery (%) BL-AP | `recovery_ba_pct` | float | % |

#### Table 2: Monolayer Integrity
| Column | Field Name | Type | Unit |
|--------|------------|------|------|
| Compound ID | `compound_id` | string | - |
| TEER A-B | `teer_ab` | float | Ohm x cm^2 |
| TEER B-A | `teer_ba` | float | Ohm x cm^2 |
| LY Leakage A-B (%) | `ly_leakage_ab_pct` | float | % |
| LY Leakage B-A (%) | `ly_leakage_ba_pct` | float | % |

#### JSON Output Schema
```json
{
    "assay_type": "caco2_permeability",
    "KPI": "efflux_ratio",
    "papp_ab": 0.56,
    "papp_ba": 4.36,
    "papp_unit": "1e-6 cm/s",
    "efflux_ratio": 7.79,
    "recovery_ab_pct": 36.81,
    "recovery_ba_pct": 82.01,
    "monolayer_integrity": {
        "teer_ab": 701.06,
        "teer_ba": 757.54,
        "teer_unit": "ohm_cm2",
        "ly_leakage_ab_pct": 0.15,
        "ly_leakage_ba_pct": 0.13
    },
    "permeability_class": "low",
    "efflux_substrate": true,
    "flags": [],
    "source": {
        "vendor": "NCU",
        "file": "ADME-NCU-Caco-2 Permeability-20231219.xlsx",
        "date": "2023-12-19"
    }
}
```

**Derived fields**:
- `permeability_class`: "high" if Papp A-B > 10, "medium" if 1-10, "low" if < 1
- `efflux_substrate`: true if efflux_ratio > 2

---

## Data Model Integration

### Mapping to Existing Models

The ADME import will use existing Assays models:

```
Assay (one per import file)
  |-- protocol -> Protocol (ADME-specific, e.g., "NCU Liver Microsome Stability")
  |-- data_file -> Original Excel file
  +-- data_series -> [DataSeries, ...]
                      |-- compound -> Compound (matched by ID)
                      |-- compound_name -> Original ID from file
                      +-- analysis -> AnalysisResult
                                      +-- results -> JSON (schema above)
```

### Protocol Records

Create Protocol records for each assay type:

| Protocol Name | Slug | fitting_method |
|---------------|------|----------------|
| NCU Liver Microsome Stability | `ncu-lm` | `None` |
| NCU Blood/Serum Stability | `ncu-bs` | `None` |
| NCU GSH Stability | `ncu-gsh` | `None` |
| NCU Caco-2 Permeability | `ncu-caco2` | `None` |

`fitting_method = None` indicates pre-calculated results (no curve fitting needed).

### DataSeries Records

For assays with multiple conditions (species, +/- cofactors), create **separate DataSeries** for each:

**Example for LM assay with Human and Mouse**:
- DataSeries 1: Compound X, Human, with cofactors
- DataSeries 2: Compound X, Human, without cofactors
- DataSeries 3: Compound X, Mouse, with cofactors
- DataSeries 4: Compound X, Mouse, without cofactors

**Alternative**: Single DataSeries per compound with all conditions nested in JSON.

**Recommendation**: Use **one DataSeries per compound per species**, with +/- conditions nested. This balances granularity with query simplicity.

---

## Implementation Plan

### Phase 1: Core Parser Infrastructure

1. **Create parser module**: `apps/compounds/assays/importers/`
   ```
   importers/
   |-- __init__.py
   |-- base.py          # Base parser class
   |-- ncu/
   |   |-- __init__.py
   |   |-- base.py      # NCU common utilities
   |   |-- lm.py        # Liver Microsome parser
   |   |-- bs.py        # Blood/Serum Stability parser
   |   |-- gsh.py       # GSH Stability parser
   |   +-- caco2.py     # Caco-2 Permeability parser
   +-- registry.py      # Parser registry for auto-detection
   ```

2. **Base parser class**:
   ```python
   class ADMEParser:
       vendor: str
       assay_type: str

       def detect(self, filepath: Path) -> bool:
           """Return True if this parser can handle the file."""

       def parse(self, filepath: Path) -> list[ParsedResult]:
           """Parse file and return structured results."""

       def validate(self, results: list[ParsedResult]) -> list[ValidationError]:
           """Validate parsed results."""
   ```

3. **ParsedResult dataclass**:
   ```python
   @dataclass
   class ParsedResult:
       compound_id: str
       species: str | None
       assay_format: str | None
       results: dict  # JSON-serializable
       source_row: int
   ```

### Phase 2: Django Management Command

Create `apps/compounds/assays/management/commands/import_adme.py`:

```bash
# Preview mode (dry run)
python manage.py import_adme /path/to/ADME-NCU-LM-20231218.xlsx --preview

# Import with compound matching
python manage.py import_adme /path/to/ADME-NCU-LM-20231218.xlsx

# Import with auto-create for unmatched compounds
python manage.py import_adme /path/to/ADME-NCU-LM-20231218.xlsx --create-compounds

# Specify protocol explicitly
python manage.py import_adme /path/to/file.xlsx --protocol ncu-lm
```

**Features**:
- Auto-detect assay type from filename
- Preview mode shows what would be imported
- Compound matching by ID (exact match, then fuzzy)
- Option to create placeholder compounds for unmatched IDs
- Transaction-wrapped import (all or nothing)
- Detailed logging of import results

### Phase 3: Compound Matching

1. **Exact match**: Look up `Compound` by `formatted_id` or `name`
2. **Normalized match**: Strip prefixes, normalize case
3. **Fuzzy match**: Levenshtein distance for typos (with confirmation)
4. **Control compounds**: Skip known controls (Verapamil, Metoprolol, Propantheline, etc.)

```python
CONTROL_COMPOUNDS = {
    'verapamil', 'metoprolol', 'digoxin', 'propantheline',
    'afatinib', 'warfarin', 'testosterone'
}
```

### Phase 4: Protocol Setup

Data migration to create Protocol records:

```python
def create_adme_protocols(apps, schema_editor):
    Protocol = apps.get_model('assays', 'Protocol')

    protocols = [
        {
            'name': 'NCU Liver Microsome Stability',
            'slug': 'ncu-lm',
            'description': 'Metabolic stability in human/mouse liver microsomes',
            'fitting_method': None,
            'fitting_parameters': {
                'vendor': 'NCU',
                'assay_type': 'liver_microsome_stability',
                'kpi': 't1_2_min'
            }
        },
        # ... other protocols
    ]
```

### Phase 5: Frontend Display (Optional)

Extend the Assays frontend to display ADME results appropriately:
- Recognize `assay_type` in results JSON
- Render appropriate visualization (table vs. time-course chart)
- Show derived classifications (permeability class, efflux substrate, etc.)

---

## Special Value Handling

| Value | Meaning | Storage | Flag |
|-------|---------|---------|------|
| `BLOD` | Below Limit of Detection | `null` | `"below_lod"` |
| `-` | Not measured | `null` | - |
| `inf` | Infinite (stable) | `null` | `"stable"` |
| `> 686.01` | Exceeds measurement window | `686.01` | `"exceeds_window"` |
| `< 1.0` | Below quantification | `1.0` | `"below_loq"` |

---

## Quality Control Flags

Automatically set flags based on results:

| Flag | Condition | Assay Types |
|------|-----------|-------------|
| `"stable"` | t1/2 = inf or > 120 min | LM, BS, GSH |
| `"unstable"` | t1/2 < 10 min | LM, BS, GSH |
| `"below_lod"` | Any value is BLOD | All |
| `"low_recovery"` | Recovery < 70% | Caco-2 |
| `"efflux_substrate"` | Efflux ratio > 2 | Caco-2 |
| `"high_efflux"` | Efflux ratio > 10 | Caco-2 |
| `"gsh_reactive"` | t1/2(+GSH) << t1/2(-GSH) | GSH |
| `"high_clearance"` | CLint > 100 uL/min/mg | LM |
| `"poor_permeability"` | Papp A-B < 1 | Caco-2 |

---

## Error Handling

### Parse Errors
- Missing required columns -> skip file, report error
- Malformed values -> skip row, log warning
- Empty compound ID -> skip row

### Compound Matching Errors
- No match found -> log warning, optionally create placeholder
- Multiple matches -> require manual resolution
- Control compound -> skip (don't import)

### Database Errors
- Duplicate import -> check for existing Assay with same file, warn
- Constraint violation -> rollback, report

---

## Testing Strategy

1. **Unit tests**: Parser functions with sample data
2. **Integration tests**: Full import cycle with test database
3. **Fixture files**: Anonymized versions of real NCU files

Test files location: `apps/compounds/assays/importers/tests/fixtures/`

---

## Future Extensions

### Additional Vendors
The architecture supports adding new vendors:
1. Create new parser module under `importers/{vendor}/`
2. Register parsers in `registry.py`
3. Add Protocol records for vendor-specific assay types

### Additional Assay Types
Common ADME assays that may be added:
- Plasma Protein Binding (PPB)
- CYP Inhibition (1A2, 2C9, 2C19, 2D6, 3A4)
- CYP Induction
- hERG Inhibition
- Solubility (kinetic, thermodynamic)
- LogD / LogP

### Bulk Import
Support for importing multiple files:
```bash
python manage.py import_adme /path/to/adme_files/*.xlsx --batch
```

### API Endpoint
REST endpoint for programmatic import:
```
POST /api/assays/import/adme/
Content-Type: multipart/form-data
file: <xlsx file>
```

---

## Appendix: Sample Data

### Liver Microsome (Table 1)
```
Compound ID          | Species | t1/2 (min) | CLint      | Scaled CLint | CLH    | ER
---------------------|---------|------------|------------|--------------|--------|------
Verapamil            | Human   | 5.45       | 255.19     | 262.34       | 19.44  | 0.93
                     | Mouse   | 2.39       | 579.53     | 2396.94      | 86.74  | 0.96
PH-NTU-03-L1-559-0   | Human   | 53.46      | 25.95      | 26.67        | 11.75  | 0.56
                     | Mouse   | 68.82      | 20.14      | 83.31        | 43.26  | 0.48
```

### Caco-2 Permeability (Table 1)
```
Compound ID          | Papp A-B | Papp B-A | Efflux Ratio | Recovery A-B | Recovery B-A
---------------------|----------|----------|--------------|--------------|-------------
Metoprolol           | 24.34    | 23.31    | 0.96         | 112.73%      | 100.44%
Digoxin              | 0.38     | 14.34    | 37.38        | 91.15%       | 97.72%
PH-NTU-03-L1-559-0   | 0.56     | 4.36     | 7.79         | 36.81%       | 82.01%
```

---

## References

- Django Assays Models: `apps/compounds/assays/models.py`
- AnalysisResult JSON schema: See `AnalysisResult.results` field
- Existing analysis flow: `apps/compounds/assays/analysis.py`
