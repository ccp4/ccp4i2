# Stub/Implementation Inheritance Pattern

**Date**: 2025-10-27
**Status**: ✅ DOCUMENTED

## The Problem

In the stub/implementation pattern, **stub classes inherit from stubs**, but **full-fat implementations need to inherit methods from their full-fat parents**, not just from stubs.

### Example Hierarchy

```
Stubs (generated):
  CMiniMtzDataFileStub(CDataFile)
    ├─ CObsDataFileStub(CMiniMtzDataFileStub)
    ├─ CPhsDataFileStub(CMiniMtzDataFileStub)
    └─ CFreeRDataFileStub(CMiniMtzDataFileStub)

Implementations (manual):
  CMiniMtzDataFile(CMiniMtzDataFileStub)
    └─ has methods: foo(), bar(), etc.

  CObsDataFile(CObsDataFileStub)  ← ❌ Can't see CMiniMtzDataFile.foo()!
```

**Problem**: `CObsDataFile` inherits from `CObsDataFileStub`, which inherits from `CMiniMtzDataFileStub` (the stub). It **never sees** the full-fat `CMiniMtzDataFile` and its methods!

## The Solution: Multiple Inheritance

Full-fat implementation classes should inherit from **both** their stub AND their full-fat parent:

```python
class CObsDataFile(CObsDataFileStub, CMiniMtzDataFile):
    """
    Inherits from:
    - CObsDataFileStub: Metadata and structure
    - CMiniMtzDataFile: Shared full-fat methods
    """
    pass
```

### Method Resolution Order (MRO)

Python uses C3 linearization to resolve the diamond inheritance:

```
CObsDataFile.__mro__ =
  0: CObsDataFile              # Self
  1: CObsDataFileStub          # Primary parent (stub)
  2: CMiniMtzDataFile          # Secondary parent (full-fat) ✓
  3: CMiniMtzDataFileStub      # Base stub
  4: CMtzDataFileStub
  5: CDataFile
  6: CData
  7: HierarchicalObject
```

Now `CObsDataFile` **can access** all methods from `CMiniMtzDataFile`!

## Implementation in ccp4i2

### File: `core/CCP4XtalData.py`

```python
# ============================================================================
# Base Class (must be defined first)
# ============================================================================

class CMiniMtzDataFile(CMiniMtzDataFileStub):
    """Base class for all MTZ mini-file types."""

    def shared_method(self):
        """This method is available to all descendants."""
        pass

# ============================================================================
# Descendant Classes (defined after parent)
# ============================================================================

class CObsDataFile(CObsDataFileStub, CMiniMtzDataFile):
    """
    Observed diffraction data.

    Inherits from:
    - CObsDataFileStub: Metadata and structure
    - CMiniMtzDataFile: Shared full-fat methods
    """

    def as_IPAIR(self, work_directory=None):
        # Can call shared_method() from CMiniMtzDataFile!
        self.shared_method()
        ...

class CPhsDataFile(CPhsDataFileStub, CMiniMtzDataFile):
    """
    Phase probability data.

    Inherits from:
    - CPhsDataFileStub: Metadata and structure
    - CMiniMtzDataFile: Shared full-fat methods
    """

    def as_HL(self, work_directory=None):
        # Can also call shared_method()!
        self.shared_method()
        ...
```

## Critical Constraint: Definition Order

**Classes can only inherit from classes defined earlier in the file.**

### Problem: Classes Defined Before Parent

```python
# ❌ THIS DOESN'T WORK

class CFreeRDataFile(CFreeRDataFileStub, CMiniMtzDataFile):
    # ERROR: CMiniMtzDataFile not defined yet!
    pass

# ... 200 lines later ...

class CMiniMtzDataFile(CMiniMtzDataFileStub):
    pass
```

**Error**: `NameError: name 'CMiniMtzDataFile' is not defined`

### Solution: Single Inheritance for Early Classes

Classes defined **before** their full-fat parent must use single inheritance:

```python
class CFreeRDataFile(CFreeRDataFileStub):
    """
    Free R flags.

    NOTE: Cannot inherit from CMiniMtzDataFile due to definition order.
    CMiniMtzDataFile is defined later in the file.
    """
    pass

class CMapCoeffsDataFile(CMapCoeffsDataFileStub):
    """
    Map coefficients.

    NOTE: Cannot inherit from CMiniMtzDataFile due to definition order.
    CMiniMtzDataFile is defined later in the file.
    """
    pass

# ... Later in file ...

class CMiniMtzDataFile(CMiniMtzDataFileStub):
    """Now defined - later classes can inherit from this."""
    pass

# ... Even later ...

class CObsDataFile(CObsDataFileStub, CMiniMtzDataFile):
    """✓ This works - CMiniMtzDataFile is already defined."""
    pass
```

## Current Status in CCP4XtalData.py

| Class | Line | Inherits from CMiniMtzDataFile? | Reason |
|-------|------|--------------------------------|--------|
| CFreeRDataFile | 329 | ❌ No | Defined before parent (line 532) |
| CMapCoeffsDataFile | 452 | ❌ No | Defined before parent (line 532) |
| **CMiniMtzDataFile** | **532** | **N/A (parent)** | **Base class** |
| CObsDataFile | 653 | ✅ Yes | Defined after parent |
| CPhsDataFile | 799 | ✅ Yes | Defined after parent |

## When to Use This Pattern

### ✅ Use Multiple Inheritance When:

1. **Sharing implementation methods** across sibling classes
   ```python
   class CMiniMtzDataFile:
       def validate_mtz(self): ...
       def get_column_info(self): ...

   # All descendants get these methods!
   class CObsDataFile(CObsDataFileStub, CMiniMtzDataFile): pass
   class CPhsDataFile(CPhsDataFileStub, CMiniMtzDataFile): pass
   ```

2. **Class is defined AFTER its full-fat parent** in the file

3. **You want to avoid code duplication** across related classes

### ❌ Don't Use Multiple Inheritance When:

1. **Class is defined BEFORE its full-fat parent**
   - Use single inheritance: `class X(XStub)`
   - Add note explaining why

2. **No shared methods exist** in the parent class
   - Current CMiniMtzDataFile is empty (`pass`)
   - Multiple inheritance not needed (yet)

3. **Methods are already in CDataFile** base class
   - Example: `_get_conversion_output_path()` is in CDataFile
   - All CDataFile descendants get it automatically

## Advantages

### 1. Code Reuse
```python
class CMiniMtzDataFile:
    def validate_columns(self, required_cols):
        """Shared validation logic."""
        mtz = gemmi.read_mtz_file(self.getFullPath())
        actual_cols = {col.label for col in mtz.columns}
        missing = set(required_cols) - actual_cols
        if missing:
            raise ValueError(f"Missing columns: {missing}")
        return True

# All children can use this!
class CObsDataFile(CObsDataFileStub, CMiniMtzDataFile):
    def as_IPAIR(self, work_directory=None):
        self.validate_columns(['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'])
        # ... conversion logic ...
```

### 2. Centralized Updates

Update the parent once → all children get the update:

```python
# Add new method to parent
class CMiniMtzDataFile:
    def get_space_group(self):
        mtz = gemmi.read_mtz_file(self.getFullPath())
        return mtz.spacegroup.hm

# All descendants automatically have it!
obs = CObsDataFile()
print(obs.get_space_group())  # Works!

phs = CPhsDataFile()
print(phs.get_space_group())  # Works!
```

### 3. Maintains Stub Metadata

The stub is still the **primary parent**, so all metadata is preserved:

```python
class CObsDataFile(CObsDataFileStub, CMiniMtzDataFile):
    # Still gets all CONTENT_FLAGS, SUBTYPES, etc. from stub!
    pass

print(CObsDataFile.CONTENT_FLAG_FMEAN)  # 4 (from stub)
print(CObsDataFile.SUBTYPE_OBSERVED)    # 1 (from stub)
```

## Diamond Inheritance Diagram

```
                    CDataFile
                        ↑
                        |
                CMiniMtzDataFileStub
                    ↑       ↑
                    |       |
    CObsDataFileStub        CMiniMtzDataFile
                ↑               ↑
                |               |
                +---------------+
                        |
                  CObsDataFile
```

Python's MRO ensures no method is called twice and no conflicts occur.

## Verification

Check MRO at runtime:

```python
from core.CCP4XtalData import CObsDataFile

print("CObsDataFile MRO:")
for i, cls in enumerate(CObsDataFile.__mro__):
    print(f"  {i}: {cls.__name__}")

# Output:
#   0: CObsDataFile
#   1: CObsDataFileStub        ← Primary parent (stub)
#   2: CMiniMtzDataFile        ← Full-fat parent (methods)
#   3: CMiniMtzDataFileStub    ← Base stub
#   4: CMtzDataFileStub
#   5: CDataFile
#   6: CData
#   7: HierarchicalObject
```

## Testing

All tests pass with multiple inheritance:

```bash
$ python -m pytest tests/test_cpluginscript_makehklin.py -v
==================== 19 passed in 0.13s ====================
```

## Alternative Solutions (Not Recommended)

### Alternative 1: Reorder File

Move `CMiniMtzDataFile` to the beginning of the file, before all descendants.

**Problems**:
- Breaks generator output order
- Stubs and implementations would be interleaved
- Harder to maintain

### Alternative 2: Duplicate Methods

Copy methods into each subclass.

**Problems**:
- Code duplication
- Maintenance nightmare
- Violates DRY principle

### Alternative 3: Mixin Class

Create a separate mixin for shared methods.

**Problems**:
- Extra complexity
- Multiple inheritance still needed
- Doesn't solve ordering issue

## Conclusion

**Use multiple inheritance for full-fat implementation classes when:**
1. The full-fat parent is defined earlier in the file
2. You want to share methods across sibling classes
3. You need both metadata (from stub) and methods (from parent)

**Pattern**:
```python
class ChildClass(ChildClassStub, ParentClass):
    """
    Inherits from:
    - ChildClassStub: Metadata and structure
    - ParentClass: Shared full-fat methods
    """
    pass
```

This pattern is **already in use** for:
- ✅ `CObsDataFile` (line 653)
- ✅ `CPhsDataFile` (line 799)

And **documented for future use** when `CMiniMtzDataFile` gets shared methods.
