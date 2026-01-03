# TODO: splitHklout() Type-Based Column Matching

## Current Implementation (v1)
- Auto-infers output column names from `CONTENT_SIGNATURE_LIST[0]`
- Works for simple cases like FreeR (single column)
- Does **NOT** validate that input column types match the signature

## Future Enhancement: Type-Based Signature Matching

### Problem
For multi-column file types with multiple signatures (e.g., IPAIR), we need to:
1. Read the **source MTZ** to get input column types
2. Match those types against **all signatures** in `CONTENT_SIGNATURE_LIST`
3. Use the **matching signature** for relabeling (not just `[0]`)

### Example: IPAIR Files
```python
# Hypothetical IPAIR signatures
CONTENT_SIGNATURE_LIST = [
    ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)'],  # Standard names
    ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'],  # Alternative names
]
```

**Current behavior**: Always uses first signature `['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)']`

**Desired behavior**:
- If input columns are `['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus']` with types `[K, M, K, M]`
- Match types `[K, M, K, M]` against signature definitions
- Find matching signature based on column types
- Use that signature for relabeling

### Implementation Strategy

1. **Extend CONTENT_SIGNATURE_LIST format**:
   ```python
   CONTENT_SIGNATURE_LIST = [
       {
           'names': ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)'],
           'types': ['K', 'M', 'K', 'M'],  # NEW: type specifications
       },
       {
           'names': ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'],
           'types': ['K', 'M', 'K', 'M'],
       },
   ]
   ```

2. **Modify `splitHklout()` logic**:
   ```python
   # Read source MTZ to get input column types
   import gemmi
   mtz = gemmi.read_mtz_file(inFile)
   input_cols = [col for col in mtz.columns if col.label in input_col_names]
   input_types = [col.type for col in input_cols]

   # Match against signatures
   for sig in file_obj.CONTENT_SIGNATURE_LIST:
       if sig['types'] == input_types:
           output_col_string = ','.join(sig['names'][0])  # Use first matching
           break
   ```

3. **Maintain backward compatibility**:
   - Keep supporting old format: `[['FREER']]` (list of lists)
   - New format: `[{'names': ['FREER'], 'types': ['I']}]` (list of dicts)
   - Auto-detect format and handle both

### Benefits
- Correct relabeling for complex multi-column types
- More robust matching based on MTZ column types, not just names
- Handles cases where programs generate different column names but same types

### Migration Path
1. ✅ **Phase 1 (Current)**: Simple name-based relabeling using `[0]`
2. **Phase 2**: Add type-based matching for complex file types
3. **Phase 3**: Migrate all `CONTENT_SIGNATURE_LIST` definitions to new format

### Files to Update (Phase 2)
- `core/CCP4PluginScript.py`: `splitHklout()` method
- `core/CCP4XtalData.py`: Extend `CONTENT_SIGNATURE_LIST` format for complex types
- Add tests for type-based matching in `tests/test_splitHklout.py`

### Related Issues
- See conversation about FreeR column standardization
- User feedback: "The relabelling will depend on the types of the column to be extracted"
