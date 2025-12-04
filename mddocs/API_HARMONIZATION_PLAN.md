# CData API Harmonization Plan

## Overview

This document outlines the plan to harmonize the "new" CData implementation (modern hierarchy/metadata/event system) with the "old" CData API documented in `rstdocs/source/developers/data_classes.rst`.

## Current Status

### New Implementation (What We Have)

- **Location**: `core/base_object/`
- **Architecture**:
  - `HierarchicalObject` → `CData` → fundamental types & generated classes
  - Signal/event system (replacing Qt signals)
  - Modern metadata system
  - Value state tracking (NOT_SET, DEFAULT, EXPLICITLY_SET)
  - Property-based value validation

### Old API (What's Expected)

- **Qt-based**: Originally subclassed `QtCore.QObject`
- **XML-based**: Heavy use of eTree for serialization
- **Container-focused**: CContainer with dynamic contents from DEF files
- **Class attributes**: CONTENTS, QUALIFIERS, PYTHONTYPE, etc.

---

## API Mapping Analysis

### Phase 1: Core Data Methods ✅ (Already Implemented)

| Old API Method                 | New Implementation                                 | Status         | Notes                          |
| ------------------------------ | -------------------------------------------------- | -------------- | ------------------------------ |
| `__init__(parent, qualifiers)` | `__init__(parent, name, **kwargs)`                 | ✅ Implemented | Modern uses kwargs             |
| `parent()`                     | `get_parent()` via HierarchicalObject              | ✅ Implemented |                                |
| `objectName()`                 | `name` property                                    | ✅ Implemented | Direct attribute access        |
| `objectPath()`                 | **❌ Missing**                                     | 🔨 Need to add | Path with underscore separator |
| `set(data)`                    | `set()` in CDataFile, value setter in fundamentals | ⚠️ Partial     | Need standardization           |
| `get(name=None)`               | `value` property / direct access                   | ⚠️ Partial     | Need dict-style get()          |
| `unSet()`                      | **❌ Missing**                                     | 🔨 Need to add | Reset to NOT_SET state         |
| `isSet()`                      | `isSet()`                                          | ✅ Implemented | In fundamental types           |
| `validity()`                   | Validation in setters                              | ⚠️ Partial     | Need CErrorReport return       |
| `fix()`                        | **❌ Missing**                                     | 📝 Optional    | Auto-fix invalid data          |
| `setDefault()`                 | **❌ Missing**                                     | 🔨 Need to add | Set to default value           |

### Phase 2: Qualifier Methods ⚠️ (Partially Implemented)

| Old API Method                      | New Implementation         | Status      | Notes                 |
| ----------------------------------- | -------------------------- | ----------- | --------------------- |
| `qualifiers(name, default, custom)` | `get_qualifier()`          | ⚠️ Partial  | Need dict return mode |
| `qualifiersDefinition(name)`        | Class metadata             | ⚠️ Partial  | Need accessor method  |
| `qualifiersOrder()`                 | **❌ Missing**             | 📝 Optional | For GUI ordering      |
| `pythonType()`                      | `PYTHONTYPE` in decorators | ⚠️ Partial  | Need accessor method  |

### Phase 3: XML/Serialization Methods ❌ (Not Implemented)

| Old API Method                | New Implementation | Status               | Notes                      |
| ----------------------------- | ------------------ | -------------------- | -------------------------- |
| `getEtree()`                  | **❌ Missing**     | 🔨 **High Priority** | Export to eTree/XML        |
| `setEtree(element)`           | **❌ Missing**     | 🔨 **High Priority** | Import from eTree/XML      |
| `getQualifiersEtree()`        | **❌ Missing**     | 🔨 **High Priority** | Export qualifiers to XML   |
| `setQualifiersEtree(element)` | **❌ Missing**     | 🔨 **High Priority** | Import qualifiers from XML |

### Phase 4: Signal/Event Methods ✅ (Modernized)

| Old API Method      | New Implementation                  | Status        | Notes                      |
| ------------------- | ----------------------------------- | ------------- | -------------------------- |
| `emitDataChanged()` | Signal system in HierarchicalObject | ✅ Modernized | No Qt dependency           |
| Qt signals/slots    | Modern signal system                | ✅ Modernized | Pure Python implementation |

### Phase 5: Container Methods ⚠️ (Partially Implemented)

| Old API Method                         | New Implementation         | Status           | Notes                  |
| -------------------------------------- | -------------------------- | ---------------- | ---------------------- |
| `loadContentsFromXml(fileName)`        | **❌ Missing**             | 🔨 High Priority | Load DEF file          |
| `loadDataFromXml(fileName)`            | **❌ Missing**             | 🔨 High Priority | Load PARAMS file       |
| `saveContentsToXml(fileName)`          | **❌ Missing**             | 🔨 High Priority | Save DEF file          |
| `saveDataToXml(fileName)`              | **❌ Missing**             | 🔨 High Priority | Save PARAMS file       |
| `loadContentsFromEtree(element)`       | **❌ Missing**             | 🔨 High Priority | We have DEF XML parser |
| `loadDataFromEtree(element)`           | **❌ Missing**             | 🔨 High Priority | Need params import     |
| `saveContentsToEtree()`                | **❌ Missing**             | 🔨 High Priority |                        |
| `saveDataToEtree()`                    | **❌ Missing**             | 🔨 High Priority |                        |
| `addContent(name, cls, qualifiers)`    | `add_item()` in CContainer | ⚠️ Partial       | Need standardization   |
| `addObject(name, object, afterObject)` | **❌ Missing**             | 🔨 Need to add   |                        |
| `replaceObject(name, object)`          | **❌ Missing**             | 📝 Optional      |                        |
| `deleteObject(name)`                   | `remove_item()`            | ⚠️ Partial       |                        |
| `renameObject(oldName, newName)`       | **❌ Missing**             | 📝 Optional      |                        |
| `clear()`                              | **❌ Missing**             | 🔨 Need to add   | Remove all contents    |
| `dataOrder()`                          | **❌ Missing**             | 🔨 Need to add   | List all object names  |
| `addHeader()`                          | **❌ Missing**             | 🔨 Need to add   | For XML export         |
| `parseCommandLine()`                   | **❌ Missing**             | 📝 Optional      | For CLI usage          |

### Phase 6: Class Attributes ⚠️ (Partially Implemented)

| Old API Attribute       | New Implementation                   | Status         | Notes                |
| ----------------------- | ------------------------------------ | -------------- | -------------------- |
| `CONTENTS`              | `attributes_definition` in decorator | ✅ Modernized  | Via @cdata_class     |
| `PYTHONTYPE`            | `gui_type_hint` in decorator         | ✅ Modernized  |                      |
| `QUALIFIERS`            | `qualifiers_definition` in decorator | ✅ Modernized  |                      |
| `QUALIFIERS_ORDER`      | **❌ Missing**                       | 📝 Optional    | For GUI              |
| `QUALIFIERS_DEFINITION` | `qualifiers_definition`              | ✅ Implemented |                      |
| `ERROR_CODES`           | **❌ Missing**                       | 📝 Optional    | Task-specific errors |
| `PROPERTIES`            | Python `@property`                   | ✅ Modernized  | Native Python        |

---

## Implementation Plan

### Priority 1: Core Compatibility (Essential for Legacy Code)

#### Task 1.1: Add Missing Core Methods to CData

**File**: `core/base_object/base_classes.py`

```python
class CData(HierarchicalObject):
    def objectName(self):
        """Return object name (legacy API compatibility)."""
        return self.name

    def objectPath(self, separator='_'):
        """Return path with parent names (legacy API)."""
        path_parts = []
        obj = self
        while obj is not None:
            if obj.name:
                path_parts.insert(0, obj.name)
            obj = obj.get_parent()
        return separator.join(path_parts)

    def unSet(self):
        """Unset the data (reset to NOT_SET state)."""
        if hasattr(self, '_value_states'):
            for key in self._value_states:
                self._value_states[key] = ValueState.NOT_SET
        if hasattr(self, 'value'):
            # Reset to type default
            if isinstance(self, (CInt, CFloat)):
                super().__setattr__('_value', 0 if isinstance(self, CInt) else 0.0)
            elif isinstance(self, CBoolean):
                super().__setattr__('_value', False)
            elif isinstance(self, CString):
                super().__setattr__('_value', '')

    def setDefault(self):
        """Set data to default value from qualifiers."""
        default = self.get_qualifier('default')
        if default is not None and hasattr(self, 'value'):
            self.value = default
            if hasattr(self, '_value_states'):
                self._value_states['value'] = ValueState.DEFAULT

    def get(self, name=None):
        """Get data (legacy API compatibility).

        For simple types: returns value
        For complex types: returns dict of all values or named value
        """
        if name is not None:
            # Return specific named attribute
            return getattr(self, name, None)

        # For simple value types
        if hasattr(self, 'value'):
            return self.value

        # For containers: return dict of contents
        if isinstance(self, CContainer):
            result = {}
            for child_name in dir(self):
                if not child_name.startswith('_'):
                    child = getattr(self, child_name, None)
                    if isinstance(child, CData):
                        result[child_name] = child
            return result

        return None

    def pythonType(self):
        """Return Python type for simple classes (legacy API)."""
        if hasattr(self.__class__, 'PYTHONTYPE'):
            return self.__class__.PYTHONTYPE
        if hasattr(self, 'value'):
            return type(self.value)
        return None

    def qualifiers(self, name=None, default=True, custom=True):
        """Return qualifiers dict (legacy API).

        Args:
            name: If set, return value of that qualifier only
            default: If False, exclude default qualifiers
            custom: If False, exclude custom qualifiers
        """
        if name is not None:
            return self.get_qualifier(name)

        # Return all qualifiers
        result = {}
        if hasattr(self, '_qualifiers'):
            result.update(self._qualifiers)
        return result

    def qualifiersDefinition(self, name=None):
        """Return qualifier definitions (legacy API)."""
        if hasattr(self.__class__, 'qualifiers_definition'):
            definitions = self.__class__.qualifiers_definition
            if name is not None:
                return definitions.get(name)
            return definitions
        return {} if name is None else None
```

#### Task 1.2: Standardize set() Method

**File**: `core/base_object/fundamental_types.py`

Update all fundamental types to have consistent `set()` method:

```python
def set(self, value):
    """Set value (legacy API compatibility).

    Returns CErrorReport for compatibility (or None if no errors).
    """
    try:
        self.value = value
        return None  # Success (or create minimal CErrorReport)
    except ValueError as e:
        # Return error report
        return {'error': str(e)}  # Simplified - full CErrorReport later
```

#### Task 1.3: Add dataOrder() to CContainer

**File**: `core/base_object/base_classes.py`

```python
class CContainer(CData):
    def dataOrder(self):
        """Return list of all data object names (legacy API)."""
        result = []
        for attr_name in dir(self):
            if not attr_name.startswith('_'):
                attr = getattr(self, attr_name, None)
                if isinstance(attr, CData):
                    result.append(attr_name)
        return result

    def clear(self):
        """Remove all content from container (legacy API)."""
        for child_name in self.dataOrder():
            delattr(self, child_name)

    def deleteObject(self, name):
        """Delete object with given name (legacy API)."""
        if hasattr(self, name):
            delattr(self, name)

    def addObject(self, name, obj, afterObject=None):
        """Add existing object to container (legacy API)."""
        setattr(self, name, obj)
        if obj.get_parent() != self:
            obj.set_parent(self)
```

---

### Priority 2: XML Serialization (Critical for CCP4i2 Integration)

#### Task 2.1: Implement getEtree() / setEtree()

**File**: `core/base_object/xml_serialization.py` (NEW)

```python
import xml.etree.ElementTree as ET
from typing import Optional
from .base_classes import CData, CContainer, ValueState

class XMLSerializationMixin:
    """Mixin providing XML serialization for CData objects."""

    def getEtree(self) -> ET.Element:
        """Return eTree element representing this object (legacy API)."""
        element = ET.Element(self.name or self.__class__.__name__)

        # Add value for simple types
        if hasattr(self, 'value'):
            element.text = str(self.value)

        # Add children for containers
        if isinstance(self, CContainer):
            for child_name in self.dataOrder():
                child = getattr(self, child_name)
                if isinstance(child, CData):
                    child_elem = child.getEtree()
                    element.append(child_elem)

        return element

    def setEtree(self, element: ET.Element):
        """Parse eTree element and initialize object (legacy API)."""
        # Set name from element tag
        if element.tag and not self.name:
            self.name = element.tag

        # Set value for simple types
        if hasattr(self, 'value') and element.text:
            self.value = element.text

        # Set children for containers
        if isinstance(self, CContainer):
            for child_elem in element:
                # Would need class registry to instantiate children
                pass

    def getQualifiersEtree(self) -> ET.Element:
        """Return eTree element containing qualifiers (legacy API)."""
        qualifiers_elem = ET.Element('qualifiers')

        for key, value in self.qualifiers().items():
            qual_elem = ET.SubElement(qualifiers_elem, key)
            qual_elem.text = str(value)

        return qualifiers_elem

    def setQualifiersEtree(self, element: ET.Element):
        """Parse eTree element containing qualifiers (legacy API)."""
        for child in element:
            key = child.tag
            value = child.text
            self.set_qualifier(key, value)
```

#### Task 2.2: Integrate with DEF XML Parser

**File**: `core/task_manager/def_xml_handler.py`

Connect the parser to use `setEtree()` / `setQualifiersEtree()` methods.

---

### Priority 3: Container File I/O (For DEF/PARAMS Files)

#### Task 3.1: Add XML File Methods to CContainer

**File**: `core/base_object/base_classes.py`

```python
class CContainer(CData, XMLSerializationMixin):
    def loadContentsFromXml(self, fileName):
        """Load content definition from DEF file (legacy API)."""
        from ..task_manager.def_xml_handler import parse_def_xml_file
        result = parse_def_xml_file(fileName)
        # Copy contents from result to self
        return None  # Success (CErrorReport)

    def loadDataFromXml(self, fileName):
        """Load data values from PARAMS file (legacy API)."""
        tree = ET.parse(fileName)
        root = tree.getroot()
        self.loadDataFromEtree(root)
        return None

    def saveContentsToXml(self, fileName):
        """Save content definition to DEF file (legacy API)."""
        # Would generate DEF XML from current structure
        pass

    def saveDataToXml(self, fileName):
        """Save data values to PARAMS file (legacy API)."""
        element = self.saveDataToEtree()
        tree = ET.ElementTree(element)
        tree.write(fileName, encoding='utf-8', xml_declaration=True)
        return None

    def loadContentsFromEtree(self, element):
        """Load contents from eTree element."""
        self.setEtree(element)

    def loadDataFromEtree(self, element):
        """Load data from eTree element."""
        self.setEtree(element)

    def saveContentsToEtree(self):
        """Return eTree with content definitions."""
        return self.getEtree()

    def saveDataToEtree(self):
        """Return eTree with data values."""
        return self.getEtree()

    def addHeader(self, header):
        """Add header for XML export (legacy API)."""
        self.header = header
```

---

### Priority 4: Error Reporting (For Validation)

#### Task 4.1: Create CErrorReport Class

**File**: `core/base_object/error_reporting.py` (NEW)

```python
class CErrorReport:
    """Legacy-compatible error reporting."""

    def __init__(self, error_code=0, message=''):
        self.error_code = error_code
        self.message = message

    def isOK(self):
        return self.error_code == 0

    def __bool__(self):
        return not self.isOK()

    def __str__(self):
        return self.message if self.message else 'OK'
```

#### Task 4.2: Update set() to Return CErrorReport

Update all `set()` methods to return `CErrorReport` instead of None/raising exceptions.

---

## Implementation Schedule

### Week 1: Core Compatibility

- [ ] Task 1.1: Add missing core methods
- [ ] Task 1.2: Standardize set() method
- [ ] Task 1.3: Add CContainer methods
- [ ] Test: Verify basic API compatibility

### Week 2: XML Serialization

- [ ] Task 2.1: Implement getEtree/setEtree
- [ ] Task 2.2: Integrate with DEF parser
- [ ] Task 4.1: Create CErrorReport
- [ ] Test: Round-trip XML serialization

### Week 3: Container File I/O

- [ ] Task 3.1: Add XML file methods
- [ ] Task 4.2: Update error reporting
- [ ] Test: Load/save DEF and PARAMS files

### Week 4: Integration & Testing

- [ ] Integration testing with legacy code
- [ ] Documentation updates
- [ ] Migration guide for existing code

---

## Testing Strategy

### Unit Tests

- Test each legacy method for correct behavior
- Test backward compatibility with old calling patterns
- Test new + old API methods work together

### Integration Tests

- Load actual CCP4i2 DEF files
- Load/save PARAMS files
- Test with real task definitions

### Compatibility Tests

- Run against existing CCP4i2 code
- Verify no breaking changes
- Performance benchmarks

---

## Migration Guide (For Users)

### Calling Patterns That Will Work

**Getting data:**

```python
# Old way (still works)
value = obj.get()

# New way (also works)
value = obj.value
```

**Setting data:**

```python
# Old way (returns CErrorReport)
error = obj.set(42)

# New way (raises ValueError on error)
obj.value = 42
```

**Object hierarchy:**

```python
# Old way
parent_obj = obj.parent()
path = obj.objectPath()

# New way (also supported)
parent_obj = obj.get_parent()
path = obj.objectPath()  # Both work!
```

---

## Success Criteria

1. ✅ All methods from old API documentation are implemented
2. ✅ Existing CCP4i2 task code runs without modification
3. ✅ DEF and PARAMS files can be loaded/saved
4. ✅ No performance regression
5. ✅ Modern features (signals, validation) still work
6. ✅ Clear migration path for new code

---

## Notes

- Priority is **backward compatibility** while maintaining modern improvements
- Qt dependence has been eliminated (signals are pure Python)
- Modern property-based access is preferred but old methods work too
- XML serialization is critical for CCP4i2 integration
- Error reporting needs to be standardized across old/new patterns

## Status Legend

- ✅ Implemented
- ⚠️ Partially implemented
- ❌ Not implemented
- 🔨 High priority
- 📝 Optional/lower priority
