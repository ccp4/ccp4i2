# Assay Data Architecture

## Overview

The assay module handles dose-response curve analysis for compound screening. Data flows from plate reader output through a configurable fitting pipeline to produce KPI values (IC50, Hill slope, etc.).

## Key Components

### 1. Plate Layout Configuration

Protocols now include a `plate_layout` JSON field that captures the physical arrangement of a plate:

```json
{
  "plate_format": 384,
  "controls": {
    "placement": "edge_columns",
    "max": {"columns": [1, 2], "rows": ["A", "B", "..."]},
    "min": {"columns": [23, 24], "rows": ["A", "B", "..."]}
  },
  "sample_region": {
    "start_column": 3,
    "end_column": 22,
    "start_row": "A",
    "end_row": "P"
  },
  "dilution": {
    "direction": "horizontal",
    "num_concentrations": 10
  },
  "replicate": {
    "count": 2,
    "pattern": "adjacent_rows"
  },
  "compound_source": {
    "type": "row_order"
  },
  "spreadsheet_origin": {
    "column": "A",
    "row": 1
  }
}
```

**Control placement types**:
- `edge_columns` - Controls in dedicated columns at plate edges
- `edge_rows` - Controls in dedicated rows at top/bottom
- `per_compound` - Embedded controls per compound (strip layout)

**Strip layout configuration** (for `per_compound` placement):

When controls are embedded per compound, each row contains repeating strips:
```json
{
  "controls": {"placement": "per_compound"},
  "strip_layout": {
    "strip_width": 12,
    "min_wells": 2,
    "data_wells": 8,
    "max_wells": 2,
    "strips_per_row": 2
  }
}
```
This defines a pattern: [min×2][data×8][max×2] repeated 2 times per row.

**Spreadsheet origin**: Defines where plate data starts in imported Excel files (cell A1 of the plate). The default is column A, row 1.

**Supported plate formats**: 24, 96, 384, 1536-well

**Replicate patterns**:
- `adjacent_rows` - A1,B1 are replicates of compound 1
- `adjacent_columns` - A1,A2 are replicates of compound 1
- `grouped_rows` - A1-A12, B1-B12 are replicates of compounds 1-12
- `interleaved_rows` - Alternating rows for same compound

**Compound source types**:
- `row_order` - Compounds assigned by row order in sample region
- `column_header` - Compound IDs from column headers in import file
- `row_header` - Compound IDs from row labels
- `plate_map_file` - Separate plate map defines compound positions
- `explicit_wells` - Per-well compound assignment in layout

### 2. Fitting Methods (Versioned Scripts)

The `FittingMethod` model stores Python fitting scripts as text in the database:

| Field | Purpose |
|-------|---------|
| `name` | Display name (e.g., "Four Parameter Logistic") |
| `slug` | URL-safe identifier (e.g., "four-parameter-logistic") |
| `version` | Semantic version for tracking changes |
| `script` | Python code implementing `fit()` function |
| `input_schema` | JSON Schema for input validation |
| `output_schema` | JSON Schema documenting output structure |
| `is_builtin` | True for system-provided methods |

**Why text storage (not files)?**
- Version history tracked in database
- Admin-editable through Django admin
- Migrations can seed/update scripts
- No file system dependencies for deployment

### 3. Script Interface

All fitting scripts implement a standard interface:

```python
def fit(input_data: dict) -> dict:
    """
    Args:
        input_data: {
            "concentrations": [10000, 3333, 1111, ...],  # high to low
            "responses": [95.2, 87.1, 62.3, ...],
            "controls": {"max": 100.0, "min": 2.3},
            "parameters": {"fix_hill": 1.0, ...}  # method-specific
        }

    Returns:
        {
            "ic50": 234.5,
            "hill_slope": 1.2,
            "top": 98.7,
            "bottom": 3.1,
            "r_squared": 0.994,
            "curve_points": [[x1, y1], ...],  # for plotting
            "flags": ["poor_fit", "ic50_extrapolated"],
            "kpi": "ic50",  # primary result field
            "fit_successful": True
        }
    """
```

### 4. Built-in Fitting Methods

#### Four Parameter Logistic (4PL)

Standard Hill-Langmuir sigmoidal curve:

```
y = bottom + (top - bottom) / (1 + (x / IC50)^hill)
```

**Parameters**:
- `fix_hill` - Constrain Hill coefficient (null = free)
- `fix_top` - Constrain top asymptote
- `fix_bottom` - Constrain bottom asymptote

**Quality flags**:
- `poor_fit` - R² < 0.8
- `ic50_extrapolated` - IC50 outside measured range
- `low_dynamic_range` - |top - bottom| < 10
- `unusual_hill_slope` - Hill < 0.3 or > 5
- `incomplete_top/bottom` - Curve doesn't reach asymptotes

## Frontend Upload Wizard

The `AssayUploadDrawer` component provides a 3-step wizard for adding assays:

### Step 1: Upload Excel
- Drag-and-drop or browse for Excel file (.xlsx, .xls, .csv)
- Parses file into a 2D cell grid using xlsx library
- Validates file format and shows row/column count

### Step 2: Review Extraction
- Applies protocol's plate layout to extract data series
- Uses `spreadsheet_origin` to locate plate data in the spreadsheet
- For each row in the sample region:
  - Extracts data values from sample columns
  - Extracts min/max control values from control regions
  - Identifies compound names (from row headers if configured)
- Shows validation status (missing data, missing controls)
- Displays extraction summary with valid/invalid series counts

### Step 3: Create Assay
- Configure target, lab book number, page number, comments
- **Run curve fitting analysis** toggle (enabled by default)
  - When enabled, automatically runs the protocol's fitting method
  - Uses the analysis method associated with the protocol
- Creates assay record with uploaded file
- Creates data series records for each extracted row
- Optionally triggers analysis pipeline

## Data Flow

```
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────┐
│  Plate Reader   │────▶│  Import/Parse    │────▶│  DataSeries     │
│  Output (Excel) │     │  (plate_layout)  │     │  (raw values)   │
└─────────────────┘     └──────────────────┘     └─────────────────┘
                                                         │
                                                         ▼
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────┐
│  Analysis       │◀────│  FittingMethod   │◀────│  Extract for    │
│  (ic50, hill)   │     │  .fit()          │     │  fitting        │
└─────────────────┘     └──────────────────┘     └─────────────────┘
```

### Upload API Request

The upload wizard sends a multipart/form-data POST to `/api/proxy/compounds/assays/`:

```
POST /api/proxy/compounds/assays/
Content-Type: multipart/form-data

Fields:
- protocol: Protocol ID
- data_file: Excel file
- target: Target ID (optional)
- labbook_number: Lab book number (optional)
- page_number: Page number (optional)
- comments: Comments (optional)
- extracted_series: JSON array of extracted data series
- run_analysis: "true" or "false" to trigger curve fitting
```

## File Structure

```
apps/compounds/assays/
├── models.py              # FittingMethod, Protocol (with plate_layout)
├── plate_layout.py        # Layout validation and utilities
├── fitting_scripts/
│   ├── __init__.py
│   └── four_parameter_logistic.py   # 4PL implementation
└── migrations/
    ├── 0004_add_fitting_method_and_plate_layout.py
    └── 0005_seed_builtin_fitting_methods.py
```

## Security Model

Fitting scripts are **admin-managed only**:
- No user upload of arbitrary code
- Scripts reviewed before deployment
- Built-in methods seeded via migrations
- Future: Container sandboxing for additional isolation

## Legacy Compatibility

The `Protocol.analysis_method` field is preserved for existing data:
- New protocols should use `fitting_method` FK
- `get_effective_fitting_method()` provides fallback logic
- Old analysis methods map to equivalent FittingMethod records
