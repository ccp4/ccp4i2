# CData Code Generation Migration Strategy

## Problem Statement

The current code generation approach has several critical issues:

1. **Incomplete Generated Code**: Generated classes are "stubs" that require manual patching
2. **Missing Imports**: Type aliases like `CUUID`, `CFilePath`, `CProjectId` are not imported
3. **Class Definition Order**: Classes are not properly topologically sorted (dependencies must come before dependents)
4. **Messy Repository**: Multiple output directories (`core/`, `core/cdata_stubs/`, `generated_cdata_classes_by_file/`)
5. **False Starts**: Accumulation of experimental generator scripts and partial implementations
6. **Moving Target**: Need to regenerate from evolving CCP4i2 codebase without manual rework

## Current State Analysis

### Test Failure Example
```
core/CCP4Data.py:62: in <module>
    class CFollowFromJob(CUUID):
                         ^^^^^
E   NameError: name 'CUUID' is not defined
```

**Root cause**: `CUUID` is a type alias defined in `fundamental_types.py:983` but not imported.

### Generator Issues

**Current generators**:
- `generate_new_files.py` - Outputs to `core/cdata_stubs/` with incomplete imports
- `cdata_class_generator.py` - Alternative approach, outputs to `generated_cdata_classes_by_file/`
- Both produce decorator-only classes with `pass` bodies

**What's missing**:
1. Type alias imports (`CUUID`, `CFilePath`, `CProjectId`, etc.)
2. Cross-file class references (e.g., `CAsuContent` references `CAsuContentSeqList`)
3. Proper topological sorting across files
4. Complete `__init__` methods with attribute initialization

## Proposed Solution: Single-Source Production Generator

### Goal
Generate **complete, production-ready, directly-usable** Python classes that:
- Pass all tests without manual modification
- Handle the moving CCP4i2 codebase
- Maintain clean repository structure
- Provide clear audit trail of what's generated vs. manually extended

### Architecture

```
ccp4i2/
├── migration/
│   └── CData/
│       ├── cdata.json                    # Source metadata (generated from CCP4i2)
│       ├── production_generator.py       # NEW: Single production generator
│       ├── type_resolver.py              # NEW: Resolves all type references
│       └── class_graph.py                # NEW: Full dependency graph
├── core/
│   ├── base_object/                      # Hand-written base architecture (unchanged)
│   ├── generated/                        # NEW: All generated code here
│   │   ├── __init__.py                   # Exports all classes
│   │   ├── CCP4Data.py                   # Generated (complete, importable)
│   │   ├── CCP4File.py                   # Generated (complete, importable)
│   │   ├── CCP4ModelData.py              # Generated (complete, importable)
│   │   └── ...                           # All other generated files
│   └── extensions/                       # NEW: Manual extensions (if needed)
│       ├── __init__.py
│       └── README.md                     # "Add custom methods here"
└── tests/                                # Tests import from core.generated.*
```

### Key Design Decisions

#### 1. Complete Type Resolution

**Problem**: Generated code references types that aren't imported

**Solution**: Build a complete type registry before code generation

```python
class TypeResolver:
    def __init__(self, cdata_json):
        self.fundamental_types = {
            'CInt', 'CFloat', 'CString', 'CBoolean', 'CList',
            'CUUID', 'CFilePath', 'CProjectId', 'COneWord', ...
        }
        self.custom_classes = {...}  # From cdata.json

    def resolve_import(self, class_name: str, target_file: str) -> str:
        """
        Returns proper import statement for using class_name in target_file.

        Examples:
        - resolve_import("CUUID", "CCP4Data.py")
          → "from core.base_object.fundamental_types import CUUID"
        - resolve_import("CAsuContentSeqList", "CCP4ModelData.py")
          → "" (same file, no import needed)
        - resolve_import("CDataFileContent", "CCP4ModelData.py")
          → "from core.base_object.base_classes import CDataFileContent"
        """
```

#### 2. Global Topological Sort

**Problem**: Current sort is per-file; doesn't handle cross-file dependencies

**Solution**: Build complete dependency graph across all classes

```python
class ClassDependencyGraph:
    def __init__(self, cdata_json):
        self.graph = {}  # class_name → {dependencies}
        self._build_graph()

    def _build_graph(self):
        """
        Extract dependencies from:
        1. immediate_parent (inheritance)
        2. CONTENTS (attribute types)
        3. Type aliases used in qualifiers
        """

    def get_sorted_classes_by_file(self) -> Dict[str, List[Tuple[str, dict]]]:
        """
        Returns: {
            "CCP4Data.py": [(CBaseData, data), (CRange, data), ...],
            "CCP4ModelData.py": [(CAsuContent, data), ...]
        }
        Where classes within each file are topologically sorted,
        AND classes using types from other files come after those files.
        """
```

#### 3. Complete Class Generation

**Problem**: Generated classes only have decorators and `pass`

**Solution**: Generate complete class bodies

```python
def render_complete_class(name: str, data: dict, type_resolver: TypeResolver) -> str:
    """
    Generate complete class with:
    - All necessary imports
    - Complete decorator with all metadata
    - Type-annotated attributes
    - Full __init__ with proper super().__init__()
    - Docstring
    """
    lines = []

    # Resolve all imports needed for this class
    imports = type_resolver.get_imports_for_class(name, data)

    # Decorator
    lines.append(render_decorator(data))

    # Class definition
    base = data.get('immediate_parent', 'CData')
    lines.append(f"class {name}({base}):")

    # Docstring
    if doc := data.get('docstring'):
        lines.append(f'    """{doc}"""')

    # Attribute declarations with types
    for attr, info in data.get('CONTENTS', {}).items():
        attr_type = type_resolver.resolve_type(info.get('class'))
        lines.append(f"    {attr}: {attr_type}")

    # __init__ method
    lines.append(f"    def __init__(self, parent=None, name=None, **kwargs):")
    lines.append(f"        super().__init__(parent=parent, name=name, **kwargs)")

    # Initialize attributes to None (will be created by decorator)
    for attr in data.get('CONTENTS', {}).keys():
        lines.append(f"        # self.{attr} initialized by @cdata_class decorator")

    return "\n".join(lines)
```

#### 4. Single Generation Command

**Problem**: Multiple generators with different outputs

**Solution**: One authoritative generator with clear options

```bash
# Generate all production code
python migration/CData/production_generator.py

# Options:
python migration/CData/production_generator.py \
    --input migration/CData/cdata.json \
    --output core/generated/ \
    --verify  # Run basic import checks after generation
```

#### 5. Clean Separation: Generated vs. Extended

**Problem**: Hard to distinguish generated code from manual additions

**Solution**: Clear directory structure

```python
# In core/generated/__init__.py (GENERATED - DO NOT EDIT)
from .CCP4Data import *
from .CCP4File import *
from .CCP4ModelData import *
# ... all generated exports

# In core/extensions/__init__.py (Manual extensions welcome)
# Example: Add custom methods to generated classes
from core.generated import CPdbDataFile

class CPdbDataFileExtended(CPdbDataFile):
    """Extended version with custom save logic."""

    def custom_save_method(self):
        """Add your custom logic here."""
        pass
```

### Implementation Plan

#### Phase 1: Build Type Resolution System

**Files to create**:
- `migration/CData/type_resolver.py`
- `migration/CData/class_graph.py`

**Tasks**:
1. Extract all type aliases from `fundamental_types.py` (CUUID, CFilePath, etc.)
2. Parse `cdata.json` to build complete class registry
3. Implement dependency graph construction
4. Implement global topological sort

**Test**: Can resolve all type references in cdata.json

#### Phase 2: Build Production Generator

**Files to create**:
- `migration/CData/production_generator.py`

**Tasks**:
1. Implement complete class rendering with proper imports
2. Handle edge cases:
   - Circular dependencies (flag as errors)
   - Unknown types (flag as errors)
   - Base classes from `base_object/`
3. Generate file headers with warnings "AUTO-GENERATED - DO NOT EDIT"
4. Run autopep8 formatting

**Test**: Generated code has no import errors

#### Phase 3: Repository Restructure

**Tasks**:
1. Create `core/generated/` directory
2. Run production generator to populate `core/generated/`
3. Update all tests to import from `core.generated.*`
4. Run full test suite
5. Fix any remaining issues
6. Delete old directories:
   - `core/cdata_stubs/`
   - `generated_cdata_classes_by_file/`
   - Old manual implementation files in `core/CCP4*.py`

#### Phase 4: Documentation and Handoff

**Tasks**:
1. Update CLAUDE.md with new structure
2. Create `core/extensions/README.md` with examples
3. Document regeneration workflow
4. Add CI check: "git diff core/generated/ should be empty"

### Migration Checklist

**Before starting**:
- [ ] Commit current state to git
- [ ] Run tests and document current failure state
- [ ] Back up `migration/CData/cdata.json`

**Phase 1**:
- [ ] Create `type_resolver.py` with type alias extraction
- [ ] Create `class_graph.py` with dependency graph
- [ ] Write unit tests for topological sort
- [ ] Verify all types in cdata.json can be resolved

**Phase 2**:
- [ ] Create `production_generator.py`
- [ ] Generate to temporary directory and test imports
- [ ] Add `--verify` mode that tests all imports
- [ ] Compare output with current manual files

**Phase 3**:
- [ ] Create `core/generated/` and run generator
- [ ] Update test imports one file at a time
- [ ] Run test suite after each update
- [ ] Document any needed manual extensions
- [ ] Clean up old generated code

**Phase 4**:
- [ ] Update all documentation
- [ ] Add regeneration script to CI/CD
- [ ] Create "extension pattern" examples
- [ ] Final test suite run

## Expected Outcomes

### After Migration

1. **Clean regeneration**: `python migration/CData/production_generator.py` produces working code
2. **All tests pass**: No manual patching needed
3. **Clear structure**:
   - `core/base_object/` - Hand-written base classes
   - `core/generated/` - 100% generated, never manually edited
   - `core/extensions/` - Optional manual extensions
4. **Audit trail**: Git history shows clear separation
5. **Moving target**: Can regenerate from updated CCP4i2 without manual rework

### Code Quality Improvements

- No more `NameError: name 'CUUID' is not defined`
- No more circular import issues
- No more class-defined-before-its-base-class errors
- Clear generated code markers
- Tests import from consistent location

### Developer Experience

```bash
# Update from CCP4i2
python migration/CData/scan_ccp4i2.py --output migration/CData/cdata.json

# Regenerate all classes
python migration/CData/production_generator.py

# Tests pass
pytest tests/

# Commit
git add core/generated/
git commit -m "Regenerate from CCP4i2 vX.Y.Z"
```

Clean, repeatable, maintainable.

## Appendix: Type Alias Registry

**From `fundamental_types.py`** (line 983+):
```python
CUUID = CString
CFilePath = CString
CProjectId = CString
COneWord = CString
# ... extract all aliases
```

**Needed in generator**:
```python
TYPE_ALIASES = {
    'CUUID': ('CString', 'core.base_object.fundamental_types'),
    'CFilePath': ('CString', 'core.base_object.fundamental_types'),
    'CProjectId': ('CString', 'core.base_object.fundamental_types'),
    'COneWord': ('CString', 'core.base_object.fundamental_types'),
    # ... complete mapping
}
```

## Appendix: Minimal Import Set

**Every generated file needs**:
```python
"""Auto-generated from CCP4i2 metadata. DO NOT EDIT.

To extend this class, create a subclass in core/extensions/
"""

from __future__ import annotations
from typing import TYPE_CHECKING

# Base classes (always needed)
from core.base_object.base_classes import CData, CContainer, CDataFile, CDataFileContent
from core.base_object.class_metadata import cdata_class, attribute, AttributeType
from core.base_object.fundamental_types import CInt, CFloat, CString, CBoolean, CList

# Type aliases (import only if used)
from core.base_object.fundamental_types import CUUID, CFilePath, CProjectId, COneWord
# ... others as needed

# Cross-file imports (only if needed)
if TYPE_CHECKING:
    from .CCP4ModelData import CAsuContent  # Avoid circular imports
```

## Questions for Review

1. **Manual extensions**: Should we support them, or require pure generation?
   - **Recommendation**: Support extensions via subclassing in `core/extensions/`

2. **File structure**: One file per original CCP4i2 file, or split by size?
   - **Recommendation**: Maintain one-to-one mapping with CCP4i2 files for traceability

3. **Circular dependencies**: How to handle if found in cdata.json?
   - **Recommendation**: Use `TYPE_CHECKING` imports and flag as warnings

4. **Legacy compatibility**: Keep old files during migration?
   - **Recommendation**: Keep in separate branch, delete from main after tests pass
