# Assay Data Analysis - Implementation Summary

## Overview

The assay analysis system performs 4-parameter logistic (4PL) curve fitting on dose-response data to extract IC50/EC50 values and other pharmacological parameters.

## Data Flow

```
Excel File → Assay Creation → DataSeries extraction → 4PL Fitting → AnalysisResult
```

## Data Format

### extracted_data Structure

Each `DataSeries.extracted_data` array follows this format:

```
[high_signal_control, response_1, response_2, ..., response_N, low_signal_control]
```

| Position | Content | Signal Level | Curve Position |
|----------|---------|--------------|----------------|
| First (index 0) | Uninhibited control | HIGH (~8000-9000) | Top asymptote |
| Middle (indices 1 to N) | Dose-response data | Varies | Transition region |
| Last (index -1) | Fully inhibited control | LOW (~2000-3000) | Bottom asymptote |

### Concentration Series

Concentrations come from the linked `DilutionSeries` and are stored separately (not in extracted_data). They correspond to the middle elements only.

Example with 10 concentrations:
- `extracted_data` has 12 elements: [control_high, 10 responses, control_low]
- `dilution_series.concentrations` has 10 elements: [0.02, 0.06, 0.2, ..., 600] µM

## 4-Parameter Logistic Model

### Hill-Langmuir Equation

```
y = bottom + (top - bottom) / (1 + (x / IC50)^hill)
```

Where:
- `top` = maximum response (high signal, uninhibited)
- `bottom` = minimum response (low signal, fully inhibited)
- `IC50` = concentration at 50% response
- `hill` = Hill coefficient (slope steepness)

### Implementation

**File:** `apps/compounds/assays/fitting_scripts/four_parameter_logistic.py`

Key features:
- Uses scipy.optimize.curve_fit with bounds constraints
- Dynamic bounds based on actual data range
- Supports fixed parameters (fix_hill, fix_top, fix_bottom)
- Quality flags for incomplete curves, poor fits, extrapolated IC50

## Key Files

| File | Purpose |
|------|---------|
| `apps/compounds/assays/models.py` | DataSeries, AnalysisResult models |
| `apps/compounds/assays/analysis.py` | Analysis orchestration, control extraction |
| `apps/compounds/assays/fitting_scripts/four_parameter_logistic.py` | 4PL curve fitting algorithm |
| `apps/compounds/assays/fitting_scripts/tight_binding_wang.py` | Tight-binding Wang equation fitting |
| `apps/compounds/assays/views.py` | API endpoints including analyse_all action |
| `apps/compounds/frontend/components/DoseResponseChart.tsx` | Chart visualization |
| `apps/compounds/frontend/components/TightBindingParametersForm.tsx` | Tight-binding parameter configuration UI |

## API Endpoints

### Trigger Analysis

```bash
POST /compounds/assays/{id}/analyse_all/
```

Re-runs analysis on all data series in an assay.

### Get Results

```bash
GET /compounds/assays/{id}/
```

Returns assay with nested data_series, each containing analysis results.

## Analysis Results Schema

```json
{
  "status": "valid|invalid",
  "results": {
    "EC50": 195.5,
    "Hill": 0.57,
    "minVal": 6754,
    "maxVal": 9158,
    "r_squared": 0.99,
    "KPI": "ic50",
    "flags": ["incomplete_top", "incomplete_bottom"],
    "curve_points": [[x1, y1], ...],
    "fit_successful": true,
    "error": null
  }
}
```

### Quality Flags

| Flag | Meaning |
|------|---------|
| `poor_fit` | R² < 0.8 |
| `ic50_extrapolated` | IC50 outside tested concentration range |
| `low_dynamic_range` | |top - bottom| < 10 |
| `unusual_hill_slope` | Hill < 0.3 or Hill > 5 |
| `incomplete_top` | Curve doesn't reach top asymptote |
| `incomplete_bottom` | Curve doesn't reach bottom asymptote |
| `fit_error` | Fitting algorithm failed |

## Recent Fixes

### 1. Control Value Mapping (2026-01-07)

**Problem:** Controls were swapped when passed to the fitter, causing "x0 is infeasible" errors.

**Fix:** In `analysis.py`, corrected the mapping:
```python
# First element is high signal control (top of curve)
controls['max'] = raw_data[0]
# Last element is low signal control (bottom of curve)
controls['min'] = raw_data[-1]
```

### 2. Legacy Data Min/Max Inversion (2026-01-07)

**Problem:** Legacy analysis results used "max" to mean "max inhibition" (low signal), while new code uses "max" to mean "max signal" (high signal).

**Fix:** In `DoseResponseChart.tsx`, auto-detect and swap inverted values:
```typescript
if (normalizedMin > normalizedMax) {
  [normalizedMin, normalizedMax] = [normalizedMax, normalizedMin];
}
```

### 3. Hill-Langmuir Equation Direction (prior)

**Problem:** Original equation `(ec50/x)^hill` produced inverted curves.

**Fix:** Changed to `(x/ec50)^hill` for standard dose-response behavior.

## Database Relationships

```
Protocol (1) ──── (*) Assay
    │                  │
    │                  └──── (*) DataSeries ──── (1) AnalysisResult
    │                              │
    └──── (*) DilutionSeries ──────┘
```

- **Assay → DataSeries:** CASCADE (deleting assay deletes its data series)
- **DataSeries → AnalysisResult:** SET_NULL (requires explicit cleanup)
- **DataSeries → DilutionSeries:** PROTECT (dilution series are protocol-level)

## Frontend Components

### DoseResponseChart

Full-size interactive chart with:
- Log-scale X-axis (concentration)
- Scatter points for raw data
- Fitted curve line
- EC50 marker
- Axis labels with units

### DoseResponseThumb

Compact thumbnail version for table displays (120x120px default).

## Testing

```bash
# Re-analyse an assay
curl -X POST http://localhost:8000/compounds/assays/{id}/analyse_all/

# Check results
curl http://localhost:8000/compounds/assays/{id}/ | python3 -m json.tool
```

## Tight-Binding Competition Analysis (Wang Equation)

### When to Use

Use tight-binding analysis when:
- Inhibitor concentrations approach or exceed protein concentration
- Standard Cheng-Prusoff IC50-to-Ki conversion is inappropriate
- You need true Ki values for SAR analysis of potent compounds

**Rule of thumb:** If IC50 < 10 * [P]_total, use tight-binding analysis.

### Required Parameters

Configure these in the protocol's `fitting_parameters`:

| Parameter | Field | Description |
|-----------|-------|-------------|
| [P]t | `protein_conc` | Total protein concentration (nM) |
| [L]t | `ligand_conc` | Total labeled ligand/competitor concentration (nM) |
| Kd_L | `ligand_kd` | Dissociation constant of labeled ligand (nM) |

### Wang Equation

The Wang equation solves the competitive binding equilibrium exactly:

```
P + L ⇌ PL     (Kd_L)
P + I ⇌ PI     (Ki)
```

With mass balance:
```
[P]t = [P] + [PL] + [PI]
[L]t = [L] + [PL]
[I]t = [I] + [PI]
```

This leads to a cubic equation in [L]_free which is solved numerically.

Reference: Wang, Z-X. (1995) FEBS Letters 360, 111-114.

### Output

| Field | Description |
|-------|-------------|
| `Ki` | True inhibition constant (primary KPI) |
| `IC50_apparent` | Observed IC50 for reference |
| `maxVal` | Uninhibited response (top asymptote) |
| `minVal` | Fully inhibited response (bottom asymptote) |
| `r_squared` | Goodness of fit |

### Quality Flags

| Flag | Meaning |
|------|---------|
| `missing_tight_binding_params` | Required parameters not configured |
| `ki_extrapolated` | Ki outside tested concentration range |
| `weak_binding_use_standard_analysis` | Ki >> [P]t, standard analysis may be appropriate |

### Example Configuration

```json
{
  "fitting_method": "<UUID of tight-binding-wang>",
  "fitting_parameters": {
    "protein_conc": 50,
    "ligand_conc": 10,
    "ligand_kd": 5
  }
}
```

## Known Limitations

1. **Inhibition assays only:** Current implementation assumes decreasing response with increasing concentration
2. **No weighting:** All data points weighted equally in fit
3. **Single replicate:** No averaging of technical replicates before fitting
4. **Fixed KPI:** Reports IC50 (standard 4PL) or Ki (tight-binding) based on fitting method
