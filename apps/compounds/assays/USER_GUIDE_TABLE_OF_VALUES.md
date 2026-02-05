# Table of Values Import - User Guide

## Overview

The **Table of Values** import allows you to upload pre-analyzed assay data from external analysis tools. Unlike raw dose-response data imports, Table of Values data comes with KPI values (like EC50, IC50, etc.) already calculated.

## When to Use Table of Values Import

Use this import type when:
- Your data has already been analyzed in another system (e.g., GraphPad Prism, TIBCO Spotfire)
- You have summary statistics rather than raw plate reader data
- You want to import results from external collaborators

**Important**: Table of Values data cannot be re-analyzed within this system because it lacks the raw concentration/response data needed for curve fitting.

## Spreadsheet Format

Your spreadsheet must include these columns:

| Column Type | Description | Example |
|-------------|-------------|---------|
| **Compound** | Compound identifiers (batch notation supported) | `NCL-00031221`, `NCL-00031221/1` |
| **KPI** | Contains the name of the KPI column (same value in all rows) | `EC50 (uM)` |
| **Data columns** | Numeric values for each metric | `EC50`, `Hill`, `R2`, etc. |
| **Image File** (optional) | Filenames for associated plot images | `compound1_plot.png` |

### Example Spreadsheet

| Compound | KPI | EC50 (uM) | Hill | R2 |
|----------|-----|-----------|------|-----|
| NCL-00001 | EC50 (uM) | 0.45 | 1.2 | 0.98 |
| NCL-00002 | EC50 (uM) | 2.31 | 0.9 | 0.95 |
| NCL-00003 | EC50 (uM) | 15.7 | 1.1 | 0.92 |

## Units

### Supported Units

The system recognizes these unit formats:

| Category | Valid Units | Notes |
|----------|-------------|-------|
| **Molar concentration** | `nM`, `uM`, `mM`, `pM`, `M` | Use `uM` not `μM` (ASCII preferred) |
| **Percentage** | `%` | Not `pc`, `pct`, or `percent` |
| **Time** | `min`, `s`, `h` | Minutes, seconds, hours |
| **Clearance rates** | `uL/min/mg`, `mL/min/kg` | For ADME data |
| **Permeability** | `cm/s`, `1e-6 cm/s` | For Caco-2 Papp |

### How Units Are Detected

1. **From the KPI column name**: If your KPI column is named `EC50 (uM)` or `IC50 [nM]`, the unit in parentheses/brackets is automatically extracted
2. **Manual override**: You can specify a unit explicitly during import, which takes precedence

### Common Unit Errors

| Wrong | Correct | Notes |
|-------|---------|-------|
| `pc` | `%` | "pc" for percent is not recognized |
| (empty) | Specify a unit | Empty unit may cause issues |
| `percent` | `%` | Use the symbol |
| `μM` | `uM` | ASCII 'u' preferred over Unicode mu |
| `nM ` | `nM` | No trailing spaces |

If you enter an invalid unit, you'll see an error message with a suggestion for the correct format.

## Import Process

### Step 1: Upload File

1. Navigate to **Assays** > **Import Table of Values**
2. Upload your Excel (.xlsx, .xls) or CSV file
3. The system will display a preview of your data

### Step 2: Configure Columns

1. **Protocol**: Select a Table of Values protocol (or this will be pre-selected if you started from a protocol page)
2. **Target** (optional): Select the biological target
3. **Compound Column**: Select the column containing compound identifiers
4. **KPI Column**: Select the column that indicates which metric is the primary KPI
5. **KPI Unit**: Override the detected unit if needed, or leave blank to use the detected unit
6. **Image File Column** (optional): Select if your data includes plot image filenames

### Step 3: Review and Import

1. Review the import summary showing row count and detected settings
2. Click **Create Assay & Import Data**
3. You'll be redirected to the new assay detail page

## After Import

### What Works

- View compound data and KPI values
- Set validation status (Valid/Invalid) for each data series
- Upload plot images via batch upload
- View uploaded plot images in data series detail
- Export data

### What Doesn't Work for Table of Values

Because Table of Values data is pre-analyzed, these features are **not available**:

- **Re-analyse**: The Re-analyse button is hidden because there's no raw data to fit
- **Change Dilutions**: Not applicable since data isn't based on a dilution series

If you need to update the values, you should re-import the data with a new spreadsheet.

## Troubleshooting

### "Invalid unit" Error

**Problem**: You see an error like "Invalid unit 'pc'. Did you mean '%'?"

**Solution**: Use the suggested unit format. Common mistakes:
- `pc` should be `%`
- `percent` should be `%`
- `um` should be `uM` (capital M for Molar)

### Import Completes But Data Looks Wrong

**Problem**: Data imports but values seem incorrect or missing

**Solutions**:
1. Check that your KPI column contains the same value in every row
2. Verify that the KPI value matches a column header exactly
3. Ensure numeric columns contain only numbers (no units in the cell)

### Compound Not Matched

**Problem**: Compound shows as "unmatched" in the assay

**Solution**:
- Check that the compound ID format matches your registry format
- Verify the compound exists in the registry
- Batch notation format is `COMPOUND-ID/BATCH-NUMBER` (e.g., `NCL-00001/1`)

### Cannot Re-analyse Data

**This is expected behavior**. Table of Values imports contain pre-analyzed results without raw data, so they cannot be re-analyzed. If you need different analysis, re-import from your source system.

## Best Practices

1. **Include units in KPI column names**: `EC50 (uM)` is better than just `EC50`
2. **Use consistent compound naming**: Match your registry format exactly
3. **Keep one KPI per import**: Each import should have a single KPI type in the KPI column
4. **Validate before import**: Review the preview to catch issues early
5. **Include plot images**: If you have visualization from your analysis tool, include the filenames and upload them after import

## Technical Notes

- Maximum file size: 50 MB
- Supported formats: .xlsx, .xls, .csv
- Compound batch notation: `COMPOUND-ID/BATCH` (e.g., `NCL-00001/1`)
- Empty cells are treated as missing values
- Duplicate compound names create separate data series entries
