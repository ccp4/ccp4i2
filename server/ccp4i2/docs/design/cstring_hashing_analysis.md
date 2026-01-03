# CString Hashing Analysis: Identity vs Value-Based

## Problem Statement

CString objects currently use identity-based hashing (`hash(id(self))`), which means they cannot be used directly as dictionary keys to match string keys:

```python
my_dict = {"PROTEIN": some_value}
cstring_obj = CString("PROTEIN")

my_dict[cstring_obj]  # KeyError! CString hashes differently than "PROTEIN"
my_dict[str(cstring_obj)]  # Works, but requires explicit conversion
```

This has caused bugs in several places:
- `crank2_basepipe.py` - KeyError when looking up pipeline steps
- `modelASUCheck.py` - KeyError when looking up polymer types

## Current Implementation

All CData fundamental types use identity-based hashing:

```python
# From core/base_object/fundamental_types.py

class CString(CData):
    def __hash__(self):
        """Make CString hashable by object identity for use in sets and as dict keys.

        Identity-based hashing is required because:
        1. CString objects are mutable (value can change)
        2. Multiple CString objects can have the same value
        3. Children tracking uses sets which require stable hashing
        """
        return hash(id(self))
```

## Analysis: Is Reason #3 Valid?

**No.** Investigation of `hierarchy_system.py` reveals that children tracking does NOT use CData's `__hash__`:

```python
class HierarchicalObject:
    def __init__(self, ...):
        self._children: Set[weakref.ReferenceType] = set()

    def _add_child(self, child):
        child_ref = weakref.ref(child)  # Creates a NEW weakref wrapper
        self._children.add(child_ref)   # Adds the WRAPPER to the set
```

The set stores `weakref.ReferenceType` objects, not CData objects directly. Each `weakref.ref(child)` creates a new wrapper with its own identity-based hash. CData's `__hash__` method is **never invoked** during hierarchy tracking.

## Remaining Constraints

### 1. Mutability (Valid Concern)

Python's hash/eq contract states:
- If `a == b`, then `hash(a) == hash(b)`
- An object's hash must remain stable while it's in a set/dict

A CString's value can change after it's been added to a dictionary, which would violate this contract if we use value-based hashing.

### 2. Object Identity vs Value Equality

With value-based hashing, two different CString objects with the same value would be treated as the same dictionary key:

```python
a = CString("PROTEIN")
b = CString("PROTEIN")
d = {a: 1}
d[b]  # Would return 1 with value-based hashing
```

This may or may not be desired depending on use case.

## Possible Solutions

### Option A: Value-Based Hashing (Recommended)

```python
def __hash__(self):
    return hash(self.value)

def __eq__(self, other):
    if isinstance(other, str):
        return self.value == other
    if isinstance(other, CString):
        return self.value == other.value
    return NotImplemented
```

**Pros:**
- `dict[cstring_obj]` matches `dict["string"]` - natural Python behavior
- No need for explicit `str()` conversions throughout codebase
- Matches how Python's built-in str works

**Cons:**
- Mutating a CString while it's a dict key causes undefined behavior (same as any mutable-value hashable)
- Two CString objects with same value are "equal" for dict purposes

**Mitigation:** Document that CString values should not be mutated while the object is used as a dict key. This is the same contract as Python's str (which is immutable) but requires developer discipline.

### Option B: Immutable CString Variant

Create a separate `FrozenCString` class that's immutable and uses value-based hashing, similar to `frozenset` vs `set`.

**Pros:**
- Clear contract - FrozenCString is safe for dict keys
- Original CString remains mutable

**Cons:**
- Adds complexity
- Requires changes to code that wants dict-key behavior

### Option C: Status Quo with str() Conversion

Keep identity-based hashing and require `str()` conversion when using CString as dict keys.

```python
# Current workaround pattern
dict_key = str(cstring_obj)  # Explicit conversion
```

**Pros:**
- No changes to CData system
- Explicit about what's happening

**Cons:**
- Error-prone - easy to forget `str()` conversion
- Already caused multiple bugs
- Feels unnatural for Python developers

### Option D: Custom __eq__ Only

Keep identity-based `__hash__` but add value-based `__eq__` for string comparison:

```python
def __eq__(self, other):
    if isinstance(other, str):
        return self.value == other
    return self is other  # Identity for CString-CString
```

**Pros:**
- `cstring_obj == "PROTEIN"` returns True

**Cons:**
- Does NOT fix dict lookup (dict uses hash first, then eq)
- Violates hash/eq contract: `a == b` but `hash(a) != hash(b)`

## Recommendation

**Implement Option A (Value-Based Hashing)** for the following reasons:

1. The main stated reason for identity hashing (hierarchy tracking) is invalid
2. The current behavior has caused real bugs
3. Value-based hashing matches Python developer expectations
4. The mutability concern applies to any object used as a dict key

### Implementation Plan

1. Update `CString.__hash__()` to return `hash(self.value)`
2. Update `CString.__eq__()` to compare values with both str and CString
3. Consider same changes for other CData types where value comparison makes sense
4. Add unit tests for dict key behavior
5. Document the contract: don't mutate CString values while used as dict keys

### Files to Modify

- `core/base_object/fundamental_types.py` - CString class (lines 881-895)
- Possibly CInt, CFloat, CBoolean if similar behavior is desired

### Test Cases to Add

```python
def test_cstring_as_dict_key():
    d = {"PROTEIN": 1, "DNA": 2}
    cs = CString("PROTEIN")
    assert d[cs] == 1  # Should work with value-based hashing

def test_cstring_equality():
    cs1 = CString("test")
    cs2 = CString("test")
    assert cs1 == cs2  # Value equality
    assert cs1 == "test"  # String equality
    assert "test" == cs1  # Reverse comparison

def test_cstring_in_set():
    s = {"PROTEIN", "DNA"}
    cs = CString("PROTEIN")
    assert cs in s  # Should find by value
```

## Related Files

- `pipelines/crank2/script/crank2_basepipe.py` - Fixed with str() workaround
- `wrappers/modelASUCheck/script/modelASUCheck.py` - Fixed with str() workaround
- These workarounds can be removed after implementing value-based hashing

## Decision History

| Date | Decision | Rationale |
|------|----------|-----------|
| 2025-12-04 | Analysis document created | Investigating CString dict key issues |
| TBD | Implementation decision | Pending review |
