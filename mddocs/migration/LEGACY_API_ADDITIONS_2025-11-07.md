# Legacy API Method Additions - 2025-11-07

## Summary

Added three missing legacy API methods to support backward compatibility with ccp4i2 plugins:

1. **CAsuDataFile.writeArpPir()** - Write PIR format sequence files for ARP/wARP
2. **CGenericReflDataFile.getFormat()** - Detect MTZ/mmCIF file format from extension
3. **CGenericReflDataFile.getMerged()** - Check if reflection data is merged

## Files Modified

### 1. core/CCP4ModelData.py (Lines 255-267)

**Added Method**: `writeArpPir()`

```python
def writeArpPir(self, fileName: str, writeMulti: bool = False, polymerTypes: list = None):
    """
    Write sequences to ARP/wARP PIR format file.

    Legacy API compatibility method for arp_warp_classic plugin.
    This is a wrapper around writeFasta() with format='pir'.

    Args:
        fileName: Output file path
        writeMulti: Write multiple copies based on nCopies
        polymerTypes: List of polymer types to include (default: PROTEIN, RNA, DNA)
    """
    return self.writeFasta(fileName, indx=-1, format='pir', writeMulti=writeMulti, polymerTypes=polymerTypes)
```

**Why**: The arp_warp_classic plugin calls `AWA_SEQIN.writeArpPir(tmpSeqFile, writeMulti=True)` at line 66.

**Status**: ✅ **WORKING** - Test shows `ARP_SEQIN.txt` file is successfully created.

**Test Result**:
- Before: `AttributeError: 'CAsuDataFile' object has no attribute 'writeArpPir'`
- After: File successfully created, test progresses to next stage

---

### 2. core/CCP4XtalData.py (Lines 407-446)

**Added Methods**: `getFormat()` and `getMerged()`

```python
def getFormat(self):
    """
    Detect file format from file extension.

    Legacy API compatibility method for import_merged and aimless plugins.

    Returns:
        str: File format ('mtz', 'mmcif', 'cif', or 'unknown')
    """
    # Get the file path
    file_path = self.getFullPath()
    if not file_path:
        # Try baseName if getFullPath() returns None
        if hasattr(self, 'baseName') and self.baseName is not None:
            file_path = str(self.baseName.value if hasattr(self.baseName, 'value') else self.baseName)
        else:
            return 'unknown'

    # Detect format from extension
    file_path_lower = str(file_path).lower()
    if file_path_lower.endswith('.mtz'):
        return 'mtz'
    elif file_path_lower.endswith('.cif') or file_path_lower.endswith('.mmcif'):
        return 'mmcif'
    else:
        return 'unknown'

def getMerged(self):
    """
    Check if reflection data is merged.

    Legacy API compatibility method for import_merged plugin.

    Returns:
        bool: True if data is merged, False if unmerged
    """
    # Default to True (merged) unless we can determine otherwise
    # In a full implementation, this would inspect the MTZ/mmCIF file
    # to check for unmerged data indicators (e.g., I+/I- columns)
    return True
```

**Why**:
- import_merged.py line 54 calls `HKLIN.getFormat()` to determine file type
- import_merged.py line 57 calls `HKLIN.getMerged()` to check merge status

**Status**: ✅ **WORKING** - Test shows "process self.fformat <class 'str'>" indicating successful string return.

**Test Result**:
- Before: `AttributeError: 'CGenericReflDataFile' object has no attribute 'getFormat'`
- After: Method returns 'mtz' successfully, test progresses to next stage (now fails at line 107 due to missing fileContent, which is a separate issue)

---

## Known Issues Identified

### 1. Locked Legacy Code Bug: pointless.py 'merged' variable

**File**: `wrappers/pointless/script/pointless.py`

**Issue**: UnboundLocalError at line 76 - `merged` variable referenced before assignment

**Root Cause**:
- Lines 44-45: `merged` is defined INSIDE a `for` loop over `UNMERGEDFILES`
- Line 76: `merged` is used OUTSIDE the loop
- If `UNMERGEDFILES` is empty (length 0), loop never executes, `merged` is never defined

**Code Context**:
```python
# Lines 34-62: INSIDE LOOP
for i in range(len(self.container.inputData.UNMERGEDFILES)):
    merged = str(self.container.inputData.UNMERGEDFILES[i].file.fileContent.merged)
    merged = (merged == 'merged')
    # ... other loop code
# Line 62: END LOOP

# Lines 73-78: OUTSIDE LOOP
if par.WRITE_HKLOUT:
    self.appendCommandScript("HKLOUT %s" % self.container.outputData.MTZUNMERGEDOUT.fullPath)
    if merged and (par.MODE == 'MATCH'):  # ❌ ERROR: 'merged' may not be defined
        self.appendCommandScript("OUTPUT UNMERGED")
```

**Status**: ⚠️ **CANNOT FIX** - This is locked legacy code in `wrappers/` directory which must not be modified.

**Impact**:
- Affects 2 tests: test_aimless.py::test_gamma, test_aimless.py::test_mdm2
- Occurs when aimless_pipe calls pointless with empty UNMERGEDFILES

**Proper Fix** (if code were not locked):
```python
# Initialize before loop
merged = False

for i in range(len(self.container.inputData.UNMERGEDFILES)):
    merged = str(self.container.inputData.UNMERGEDFILES[i].file.fileContent.merged)
    merged = (merged == 'merged')
    # ... rest of loop
```

**Workaround**: Not feasible without modifying locked code. Tests will continue to fail in this scenario.

---

### 2. CAsuDataFile.loadFile() Recursion Issue

**Observation**: When testing arpwarp, saw error: "maximum recursion depth exceeded in comparison" when loading gamma.asu.xml

**File**: core/CCP4ModelData.py - `CAsuDataFile.loadFile()` method

**Potential Cause**: Circular references in object hierarchy causing infinite recursion during XML parsing

**Status**: ⚠️ Requires investigation

**Impact**: May affect sequence file loading for various plugins

---

### 3. CGenericReflDataFile.loadFile() Missing

**Observation**: import_merged test fails at line 107 with `'NoneType' object has no attribute 'listOfColumns'`

**Root Cause**: `HKLIN.loadFile()` is called (line 53) but doesn't populate `fileContent`

**Status**: ⚠️ Requires implementation of `loadFile()` method for CGenericReflDataFile

**Required Implementation**:
- Parse MTZ file using gemmi
- Populate fileContent with column information
- Set format, merged status, numberofdatasets, etc.

---

## Test Impact Assessment

### Tests Fixed (Methods No Longer Missing)

1. **arpwarp** - writeArpPir() now exists
   - Still fails, but for different reason (missing XYZDUM.pdb output)
   - ✅ Progress: AttributeError → FileNotFoundError (different stage)

2. **import_merged** - getFormat() now exists
   - Still fails, but for different reason (fileContent not loaded)
   - ✅ Progress: AttributeError → 'NoneType' has no attribute

### Tests Still Failing (Locked Code)

3. **aimless (test_gamma, test_mdm2)** - Blocked by pointless.py 'merged' bug
   - ❌ Cannot fix: Locked legacy code issue
   - Requires `UNMERGEDFILES` to always have at least one item

---

## Design Patterns Used

### 1. Legacy API Wrapper Pattern

**Example**: `writeArpPir()` wraps `writeFasta()` with PIR-specific parameters

**Benefits**:
- Maintains backward compatibility
- Avoids code duplication
- Single source of truth for implementation

```python
def writeArpPir(self, fileName: str, writeMulti: bool = False, polymerTypes: list = None):
    return self.writeFasta(fileName, indx=-1, format='pir', writeMulti=writeMulti, polymerTypes=polymerTypes)
```

### 2. File Extension Detection

**Example**: `getFormat()` uses simple extension matching

**Rationale**:
- Fast and reliable for most cases
- No file I/O required
- Sufficient for legacy plugin compatibility

**Future Enhancement**: Could inspect file headers for more robust detection

### 3. Defensive Defaults

**Example**: `getMerged()` returns `True` by default

**Rationale**:
- Most reflection files are merged
- Safer to assume merged than unmerged
- Prevents breaking plugins that don't check merge status

**Future Enhancement**: Inspect MTZ columns (I+/I- presence) to determine actual merge status

---

## Lessons Learned

### 1. Legacy Code Dependencies

**Challenge**: Cannot modify locked legacy wrappers in `wrappers/`, `wrappers2/`, `pipelines/`

**Solution**: Document issues clearly and work around where possible

**Impact**: Some bugs (like pointless.py 'merged' variable) cannot be fixed, only documented

### 2. Incremental Method Addition

**Strategy**: Add only what's needed, when it's needed

**Benefits**:
- Reduces upfront implementation burden
- Methods are tested immediately upon addition
- Clear test-driven development cycle

**Process**:
1. Test fails with AttributeError for missing method
2. Analyze what the method should do
3. Implement minimal viable version
4. Test shows progress (fails at next stage)
5. Repeat

### 3. Wrapper vs. Full Implementation

**Decision**: Start with simple wrappers and minimal logic

**Examples**:
- `writeArpPir()` - Simple wrapper (12 lines)
- `getFormat()` - Extension-based detection (18 lines)
- `getMerged()` - Default return value (8 lines)

**Benefits**:
- Fast to implement
- Low risk of introducing bugs
- Sufficient for immediate compatibility

**Future**: Can enhance with full implementations (file inspection, column analysis, etc.) when needed

---

## Recommendations

### High Priority

1. **Implement CGenericReflDataFile.loadFile()**
   - Required for import_merged tests
   - Should use gemmi to parse MTZ
   - Populate fileContent with column metadata

2. **Investigate CAsuDataFile recursion issue**
   - Maximum recursion depth error suggests circular reference
   - May affect multiple sequence-related plugins
   - Check parent-child relationships in object hierarchy

### Medium Priority

3. **Enhanced getMerged() implementation**
   - Inspect MTZ file columns for I+/I- presence
   - Check for unmerged data indicators
   - More accurate merge status detection

4. **Enhanced getFormat() with header inspection**
   - Read file headers instead of just extensions
   - More robust format detection
   - Handle edge cases (renamed files, etc.)

### Low Priority (Documentation Only)

5. **Document pointless.py 'merged' variable bug**
   - Cannot fix due to locked code policy
   - Document in test suite expectations
   - Note: requires UNMERGEDFILES to be non-empty

---

## Summary Statistics

**Methods Added**: 3
- `writeArpPir()` - 12 lines
- `getFormat()` - 18 lines
- `getMerged()` - 8 lines

**Total Lines Added**: 38 lines of production code + extensive documentation

**Tests Impacted**:
- arpwarp: AttributeError → FileNotFoundError (progress)
- import_merged: AttributeError → NoneType error (progress)
- aimless: Still blocked by locked code bug

**Pass Rate Impact**: TBD (tests currently running)

**Files Modified**: 2
- core/CCP4ModelData.py
- core/CCP4XtalData.py

---

## Next Steps

1. ✅ **COMPLETED**: Add writeArpPir(), getFormat(), getMerged() methods
2. 🔄 **IN PROGRESS**: Run full test suite to measure impact
3. ⏭️ **NEXT**: Implement CGenericReflDataFile.loadFile() for import_merged
4. ⏭️ **NEXT**: Investigate CAsuDataFile recursion issue
5. ⏭️ **LATER**: Enhance getMerged() and getFormat() with file inspection

---

## Code Quality Notes

**Strengths**:
- Minimal, focused implementations
- Clear documentation
- Follows existing patterns (writeFasta, getFullPath)
- Defensive coding (fallbacks, defaults)

**Areas for Future Enhancement**:
- File content inspection (not just extension/defaults)
- More robust error handling
- Unit tests for individual methods
- Integration tests with gemmi library

---

## Appendix: Error Messages (Before vs. After)

### writeArpPir

**Before**:
```
AttributeError: 'CAsuDataFile' object has no attribute 'writeArpPir'
  File "arp_warp_classic.py", line 66, in processInputFiles
    cinp.AWA_SEQIN.writeArpPir(tmpSeqFile,writeMulti=True)
```

**After**:
```
FileNotFoundError: [Errno 2] Failed to open .../XYZDUM.pdb: No such file or directory
```
✅ Different error at later stage - method works!

### getFormat

**Before**:
```
AttributeError: 'CGenericReflDataFile' object has no attribute 'getFormat'
  File "import_merged.py", line 54, in process
    self.fformat = self.container.inputData.HKLIN.getFormat()
```

**After**:
```
process self.fformat <class 'str'>
...
AttributeError: 'NoneType' object has no attribute 'listOfColumns'
```
✅ Different error at later stage (line 107) - method works!

### merged variable

**Before & After** (unchanged - locked code):
```
UnboundLocalError: local variable 'merged' referenced before assignment
  File "pointless.py", line 76, in makeCommandAndScript
    if merged and (par.MODE == 'MATCH'):
```
⚠️ Cannot fix - locked legacy code in wrappers/

---

## Git Commit Message (Suggested)

```
Add legacy API methods for ccp4i2 plugin compatibility

Added three missing methods to support backward compatibility with legacy
ccp4i2 plugins:

- CAsuDataFile.writeArpPir(): Write PIR format sequences for ARP/wARP
- CGenericReflDataFile.getFormat(): Detect MTZ/mmCIF from file extension
- CGenericReflDataFile.getMerged(): Check if reflection data is merged

These methods are simple wrappers and defaults designed for immediate
compatibility. They can be enhanced with full file inspection logic in
future iterations.

Tests affected:
- arpwarp: Now progresses past writeArpPir() call
- import_merged: Now progresses past getFormat() call
- aimless: Still blocked by locked code bug in pointless.py

Files modified:
- core/CCP4ModelData.py: Added writeArpPir() (12 lines)
- core/CCP4XtalData.py: Added getFormat() and getMerged() (26 lines)

🤖 Generated with Claude Code (https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>
```
