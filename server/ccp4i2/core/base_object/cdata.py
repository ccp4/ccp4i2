"""
CData - Base class for all CCP4i2 data objects.

This is the fundamental base class that provides:
- Hierarchical object relationships
- Value state tracking (NOT_SET, DEFAULT, EXPLICITLY_SET)
- Smart assignment (preserves CData objects, updates .value for value types)
- Metadata integration
- Qualifier system
"""

import sys
import os
from typing import Any, Dict
from enum import Enum, auto
import logging

from .hierarchy_system import HierarchicalObject

# Configure logger
logger = logging.getLogger(__name__)


class ValueState(Enum):
    """States that a value can be in."""
    NOT_SET = auto()  # Value has never been explicitly set
    DEFAULT = auto()  # Value is using a default from qualifiers
    EXPLICITLY_SET = auto()  # Value has been explicitly assigned


# Import decorator after ValueState definition
from .class_metadata import cdata_class


@cdata_class(gui_label="CData")
class CData(HierarchicalObject):
    """Base class for all CCP4i2 data objects with hierarchical relationships."""

    def __init__(self, parent=None, name=None, **kwargs):
        # Initialize hierarchical object first (it only takes parent and name)
        super().__init__(parent=parent, name=name)
        logger.debug("Initializing CData instance of type %s with name '%s'", self.__class__.__name__, name)
        # Initialize set state tracking
        self._value_states: Dict[str, ValueState] = {}
        self._default_values: Dict[str, Any] = {}
        # Flag to skip validation (used during .def.xml parsing)
        self._skip_validation: bool = False

        # Mark that hierarchy is initialized - now we can use custom setattr
        self._hierarchy_initialized = True

        # Load default values from qualifiers if available
        self._load_default_values()

        # NEW: Apply metadata-driven attribute creation
        self._apply_metadata_attributes()

        # --- Per-instance metadata copying and override ---
        # Copy class-level metadata to instance for override flexibility
        cls = self.__class__
        # Qualifiers - store in _qualifiers (private) to leave qualifiers() method free for legacy API
        if hasattr(cls, '_class_qualifiers'):
            class_qualifiers = getattr(cls, '_class_qualifiers')
            logger.debug("%s._class_qualifiers type: %s, value: %s", cls.__name__, type(class_qualifiers), class_qualifiers)

            # Handle case where qualifiers is a dict-like object
            if isinstance(class_qualifiers, dict):
                self._qualifiers = dict(class_qualifiers)
            elif hasattr(class_qualifiers, 'items') and callable(getattr(class_qualifiers, 'items', None)):
                try:
                    self._qualifiers = dict(class_qualifiers.items())
                except (AttributeError, TypeError) as e:
                    logger.error("Error calling .items() on _class_qualifiers for %s: %s (type: %s)", cls.__name__, e, type(class_qualifiers))
                    self._qualifiers = {}
            else:
                # Not a dict and doesn't have .items() - set to empty dict
                logger.warning("Class-level _class_qualifiers for %s is not dict-like: %s, setting to empty dict", cls.__name__, type(class_qualifiers))
                self._qualifiers = {}
        else:
            self._qualifiers = {}
        # Qualifiers order
        if hasattr(cls, 'qualifiers_order'):
            self.qualifiers_order = list(getattr(cls, 'qualifiers_order'))
        # Qualifiers definition
        if hasattr(cls, 'qualifiers_definition'):
            self.qualifiers_definition = dict(getattr(cls, 'qualifiers_definition'))
        # CONTENT_ORDER
        if hasattr(cls, 'CONTENT_ORDER'):
            self.CONTENT_ORDER = list(getattr(cls, 'CONTENT_ORDER'))
        # For CList: subitem
        if hasattr(cls, 'subItem'):
            self.subItem = getattr(cls, 'subItem')

        # Allow overrides via kwargs
        for meta_key in ['qualifiers', 'qualifiers_order', 'qualifiers_definition', 'CONTENT_ORDER', 'subitem']:
            if meta_key in kwargs:
                # Set qualifiers directly as a dict to avoid wrapping as CData
                if meta_key == 'qualifiers' and isinstance(kwargs[meta_key], dict):
                    self._qualifiers = kwargs.pop(meta_key)
                else:
                    setattr(self, meta_key, kwargs.pop(meta_key))

        # Apply provided attributes with hierarchy handling
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __getattribute__(self, name: str):
        """Override to check hierarchy first for CData children, then fall back to normal lookup.

        This enables storing CData children in hierarchy only (not in __dict__),
        providing O(1) access via _children_by_name cache while keeping a clean
        separation between hierarchy and instance storage.
        """
        # Use object.__getattribute__ for internal/private attributes to avoid recursion
        # Also skip for common attributes that are definitely in __dict__
        if name.startswith('_') or name in ('parent', 'children', 'signals'):
            return object.__getattribute__(self, name)

        # Check if object is initialized enough to have hierarchy
        try:
            children_by_name = object.__getattribute__(self, '_children_by_name')
        except AttributeError:
            # Object not fully initialized yet
            return object.__getattribute__(self, name)

        # O(1) lookup in hierarchy cache for CData children
        if name in children_by_name:
            child_ref = children_by_name[name]
            if child_ref is not None:
                child = child_ref()
                if child is not None:
                    return child
                else:
                    # Dead reference - clean it up
                    del children_by_name[name]

        # Fall back to normal attribute lookup (checks __dict__, class, __getattr__)
        return object.__getattribute__(self, name)

    def _load_default_values(self):
        """Load default values from qualifiers metadata."""
        # Try to get default values from metadata system
        try:
            from .metadata_system import MetadataRegistry
            metadata = MetadataRegistry.get_class_metadata(self.__class__.__name__)
            if metadata:
                for field_name, field_meta in metadata.fields.items():
                    if field_meta.default_value is not None:
                        self._default_values[field_name] = field_meta.default_value
                        self._value_states[field_name] = ValueState.NOT_SET
        except Exception:
            # If metadata not available, that's okay
            pass

    def _apply_metadata_attributes(self):
        """Apply metadata-driven attribute creation if metadata is available."""
        try:
            from .class_metadata import apply_metadata_to_instance

            apply_metadata_to_instance(self)
        except ImportError:
            # Metadata system not available, skip
            pass
        except Exception:
            # Any other error, skip silently to avoid breaking existing code
            pass

    def get_qualifier(self, key, default=None):
        """Get a qualifier value for this instance."""
        if hasattr(self, '_qualifiers') and self._qualifiers is not None:
            return self._qualifiers.get(key, default)
        return default

    def set_qualifier(self, key, value):
        """Set or override a qualifier value for this instance."""
        if not hasattr(self, '_qualifiers') or self._qualifiers is None:
            self._qualifiers = {}
        self._qualifiers[key] = value

    def qualifiers(self, key=None, default=None):
        """
        Legacy compatibility method for getting qualifiers.

        Provided for backward compatibility with legacy ccp4i2 plugin code
        that calls qualifiers('key') instead of get_qualifier('key').

        Args:
            key: Optional qualifier key to get. If None, returns all qualifiers dict.
            default: Default value if key not found

        Returns:
            If key is provided, returns the qualifier value (or default).
            If key is None, returns the _qualifiers dict.

        Example:
            # Legacy code pattern:
            label = obj.qualifiers('guiLabel')

            # All qualifiers:
            all_quals = obj.qualifiers()
        """
        if key is None:
            # Return all qualifiers
            return getattr(self, '_qualifiers', {})
        else:
            # Return specific qualifier
            return self.get_qualifier(key, default)

    def setQualifier(self, key, value):
        """
        Legacy compatibility alias for set_qualifier().

        Provided for backward compatibility with legacy ccp4i2 plugin code.
        New code should use set_qualifier() instead.
        """
        return self.set_qualifier(key, value)

    def setQualifiers(self, qualifiers_dict):
        """
        Set multiple qualifiers from a dictionary.

        Provided for backward compatibility with legacy ccp4i2 plugin code
        that calls setQualifiers() with a dictionary of qualifier key-value pairs.

        Args:
            qualifiers_dict: Dictionary of qualifier key-value pairs to set

        Example:
            obj.setQualifiers({'listMinLength': 0, 'listMaxLength': 10})
        """
        if not isinstance(qualifiers_dict, dict):
            raise TypeError(f"setQualifiers() expects a dict, got {type(qualifiers_dict).__name__}")

        for key, value in qualifiers_dict.items():
            self.set_qualifier(key, value)

    def set(self, values):
        """
        Set attributes from dict or CData object, unset others.
        Uses smart assignment to avoid overwriting CData objects.

        Args:
            values: None to clear/unset, dict of attributes, or another CData object to copy from
        """
        # Handle None: clear all children and mark as unset
        if values is None:
            for child in self.children():
                if isinstance(child, CData):
                    child.set(None)  # Recursively clear children
            self.unSet()
            return

        # If values is a CData object, extract its attributes as a dict
        if isinstance(values, CData):
            source_obj = values  # Keep reference to original CData object
            if hasattr(values, 'get') and callable(values.get):
                values_dict = values.get()
                # Filter to only include fields that are explicitly set in the source object
                # This prevents validation errors when copying unset fields with constraints
                values = {}
                for k, v in values_dict.items():
                    source_field = getattr(source_obj, k, None)
                    # Only include if the field is set in the source
                    if hasattr(source_field, 'isSet') and callable(source_field.isSet):
                        if source_field.isSet():
                            values[k] = v
                    else:
                        # For non-CData fields, include them
                        values[k] = v
            else:
                # Fallback: create dict from CData attributes
                values = {
                    k: getattr(values, k)
                    for k in dir(values)
                    if not k.startswith('_') and not callable(getattr(values, k))
                    and k not in ['parent', 'name', 'children', 'signals']
                }

        # Ensure values is a dict at this point
        if not isinstance(values, dict):
            raise TypeError(f"set() expects dict or CData, got {type(values).__name__}")

        # Get all fields (metadata-aware)
        metadata = None
        try:
            from .metadata_system import MetadataRegistry
            metadata = MetadataRegistry.get_class_metadata(self.__class__.__name__)
        except Exception:
            pass

        if metadata:
            all_fields = list(metadata.fields.keys())
        else:
            # Fallback: use hierarchical children to get field names
            # This ensures we only process actual CData children, not arbitrary __dict__ entries
            all_fields = []
            for child in self.children():
                if isinstance(child, CData):
                    child_name = child.objectName() if hasattr(child, 'objectName') else None
                    if child_name:
                        all_fields.append(child_name)

        # Smart assignment for fields
        for k in all_fields:
            if k in values:
                current = getattr(self, k, None)
                new_value = values[k]

                # Check if current has a set() method (CLists, CContainers, etc.)
                # If so, use it instead of direct assignment
                if hasattr(current, 'set') and callable(getattr(current, 'set')):
                    # Use the object's set() method for proper deep copying
                    current.set(new_value)
                # Smart assignment: if current is a simple CData value type, update its value
                elif hasattr(current, 'value') and not isinstance(new_value, type(current)):
                    current.value = new_value
                else:
                    setattr(self, k, new_value)
                self._value_states[k] = ValueState.EXPLICITLY_SET
            else:
                # For CContainer types, DON'T unset fields that aren't in the source
                # This allows merging data from a source that has fewer fields
                # (e.g., phaser_pipeline.inputData → phaser_MR_AUTO.inputData)
                # For non-container types, unset fields not in values (original behavior)
                from .ccontainer import CContainer
                if not isinstance(self, CContainer):
                    if hasattr(self, k):
                        self.unSet(k)

        # Mark the object itself as set when values are provided
        # This is important for legacy code that checks `object.isSet()`
        if values and hasattr(self, '_value_states'):
            self._value_states['value'] = ValueState.EXPLICITLY_SET

    def get(self, child_name=None):
        """Get child by name or all attributes as dict. Compatible with old CCP4i2 API.

        If child_name is provided, returns the child object with that name.
        Otherwise, returns a dict of all CData attributes and their values.

        Args:
            child_name: Optional name of child to retrieve

        Returns:
            If child_name provided: The child object with that name
            If child_name is None: Dict of attribute names to values
        """
        # Legacy API: get(childName) retrieves child by name
        if child_name is not None:
            return getattr(self, child_name, None)

        # Modern API: get() returns dict representation
        result = {}

        # Get all fields (metadata-aware)
        metadata = None
        try:
            from .metadata_system import MetadataRegistry
            metadata = MetadataRegistry.get_class_metadata(self.__class__.__name__)
        except Exception:
            pass

        if metadata:
            all_fields = list(metadata.fields.keys())
        else:
            # Fallback: use hierarchical children to get field names
            # This ensures we only process actual CData children, not arbitrary __dict__ entries
            all_fields = []
            for child in self.children():
                if isinstance(child, CData):
                    child_name = child.objectName() if hasattr(child, 'objectName') else None
                    if child_name:
                        all_fields.append(child_name)

        for field_name in all_fields:
            if not hasattr(self, field_name):
                continue

            value = getattr(self, field_name)

            # Handle CData objects recursively
            if isinstance(value, CData):
                # For CList and CContainer, return the object itself
                # (so set() can call their set() methods)
                from .fundamental_types import CList
                from .ccontainer import CContainer
                if isinstance(value, (CList, CContainer)):
                    result[field_name] = value
                # If it has a value attribute, use that
                elif hasattr(value, 'value'):
                    result[field_name] = value.value
                else:
                    # Otherwise recurse
                    result[field_name] = value.get()
            else:
                result[field_name] = value

        return result

    def update(self, values: dict):
        """Update only provided attributes. Uses smart assignment to avoid overwriting CData objects."""
        for k, v in values.items():
            current = getattr(self, k, None)
            # Smart assignment: if current is a CData value type, update its value
            if hasattr(current, 'value') and not isinstance(v, type(current)):
                if v is None:
                    # Handle None by unsetting the CData object instead of assigning None
                    # This prevents TypeError when trying to assign None to typed values like CInt
                    if hasattr(current, 'unSet'):
                        current.unSet()
                else:
                    current.value = v
            else:
                setattr(self, k, v)

            # Mark as explicitly set (or not set if v is None)
            if hasattr(self, '_value_states'):
                if v is None:
                    self._value_states[k] = ValueState.NOT_SET
                else:
                    self._value_states[k] = ValueState.EXPLICITLY_SET

    # objectPath() is now inherited from HierarchicalObject
    # It uses object_path() which properly handles objectName()

    def objectName(self) -> str:
        """Return the name of this object from the HierarchicalObject hierarchy.

        This returns the hierarchical object name (stored in _name), NOT any
        CData 'name' attribute that may be stored in __dict__['name'].

        Returns:
            The hierarchical object's name or empty string
        """
        # Access _name directly to avoid collision with CData 'name' attributes
        return getattr(self, '_name', '')

    @property
    def CONTENTS(self):
        """
        Get list of child CData objects or attributes.

        For CContainer: Returns list of child CData objects (children)
        For other CData types: Returns list of attribute names that are CData objects

        This provides a unified interface for navigating the CData hierarchy,
        similar to legacy CCP4i2's CONTENTS pattern.

        Returns:
            List of CData objects (for CContainer) or list of attribute names (for other CData)
        """
        # Import CContainer here to avoid circular imports
        from .ccontainer import CContainer

        # For CContainer, return actual children objects
        if isinstance(self, CContainer):
            return list(self.get_children())

        # For other CData types, return list of CData attribute names
        cdata_attributes = []

        # Try metadata-aware approach first
        try:
            from .metadata_system import MetadataRegistry
            metadata = MetadataRegistry.get_class_metadata(self.__class__.__name__)
            if metadata:
                for field_name, field_meta in metadata.fields.items():
                    if hasattr(self, field_name):
                        value = getattr(self, field_name)
                        if isinstance(value, CData):
                            cdata_attributes.append(field_name)
                return cdata_attributes
        except Exception:
            pass

        # Fallback: use hierarchical children
        # This ensures we only return actual CData children that are part of
        # the object hierarchy, not arbitrary attributes from __dict__
        for child in self.children():
            if isinstance(child, CData):
                # Get the child's name within this parent
                child_name = child.objectName() if hasattr(child, 'objectName') else None
                if child_name:
                    cdata_attributes.append(child_name)

        return cdata_attributes

    @property
    def CONTENTS_ORDER(self):
        """Get the order of contents for this CData object.

        This provides backward compatibility with legacy CCP4i2 code that
        expects CONTENTS_ORDER to be available on all CData objects.

        For CContainer subclasses, returns _data_order.
        For other CData subclasses, returns contents_order from metadata if available,
        otherwise returns list of child attribute names in alphabetical order.

        Returns:
            List of child names in order
        """
        # Import CContainer here to avoid circular imports
        from .ccontainer import CContainer

        # For CContainer, delegate to _data_order (already handles metadata)
        if isinstance(self, CContainer):
            return self._data_order

        # For other CData types, check metadata
        if hasattr(self, '_metadata') and hasattr(self._metadata, 'contents_order'):
            contents_order = self._metadata.contents_order
            if contents_order:
                return list(contents_order)

        # Fallback: return CONTENTS (attribute names) in alphabetical order
        return sorted(self.CONTENTS)

    def isSet(self, field_name: str = None, allowUndefined: bool = False,
              allowDefault: bool = False, allSet: bool = True) -> bool:
        """Check if a field has been explicitly set.

        Args:
            field_name: Name of the field to check. If None, checks the 'value' attribute.
            allowUndefined: If True, allow None/undefined values to be considered "set"
                          if the qualifier 'allowUndefined' is True
            allowDefault: If False, consider values that equal the default as "not set"
            allSet: For container types - if True, all children must be set;
                   if False, at least one child must be set

        Returns:
            True if field has been set according to the criteria, False otherwise
        """
        if field_name is None:
            # For container objects without a 'value' attribute, check if ANY child is set
            if not hasattr(self, 'value'):
                # Check if any child attribute is set
                for _, state in self._value_states.items():
                    if state == ValueState.EXPLICITLY_SET:
                        return True
                    if state == ValueState.DEFAULT and allowDefault:
                        return True

                # Also check CData children recursively (e.g., CList, nested CContainer)
                # This handles cases like container.UNMERGEDFILES.append(item)
                # where UNMERGEDFILES itself is set but not explicitly assigned
                # Use children() to only check actual hierarchical children,
                # not all attributes from dir(self)
                for child in self.children():
                    if isinstance(child, CData) and hasattr(child, 'isSet'):
                        if child.isSet(field_name=None, allowUndefined=allowUndefined,
                                      allowDefault=allowDefault, allSet=allSet):
                            return True

                return False
            field_name = "value"

        # If the attribute exists and is a CData object, delegate to it
        # This includes CDataFile objects which have special isSet() logic that checks baseName
        if hasattr(self, field_name):
            attr = getattr(self, field_name, None)
            if isinstance(attr, CData) and hasattr(attr, 'isSet'):
                # Delegate to the child object to check if it's set
                return attr.isSet(field_name=None, allowUndefined=allowUndefined,
                                allowDefault=allowDefault, allSet=allSet)

        # Check basic set state from parent's tracking
        state = self._value_states.get(field_name, ValueState.NOT_SET)

        # If checking for explicit set only
        if state == ValueState.NOT_SET:
            # Special handling for allowUndefined
            if allowUndefined and self.get_qualifier('allowUndefined', False):
                return True
            return False

        # If state is DEFAULT, check if we allow defaults
        if state == ValueState.DEFAULT:
            if not allowDefault:
                return False

        # Check the actual value
        if hasattr(self, field_name):
            value = getattr(self, field_name, None)

            # Handle None values
            if value is None or value is NotImplemented:
                if not allowUndefined or not self.get_qualifier('allowUndefined', False):
                    return False
                return True

            # Check if value equals default (if allowDefault is False)
            if not allowDefault:
                default = self.get_qualifier('default')
                if default is not None and value == default:
                    return False

        return True

    def _has_any_set_children(self) -> bool:
        """Check if any child CData attributes have been explicitly set.

        This is used by __setattr__ to determine if a container object
        (like CRefinementPerformance) has any meaningful data to copy,
        even if the parent's isSet() returns False because no 'value' is set.

        Returns:
            True if any child CData attribute is set, False otherwise.
        """
        # Use children() to only check actual hierarchical children,
        # not all attributes from dir(self)
        for child in self.children():
            if isinstance(child, CData) and hasattr(child, 'isSet'):
                if child.isSet():
                    return True
        return False

    def unSet(self, field_name: str = 'value'):
        """Return a field to its not-set state.

        Args:
            field_name: Name of the field to unset (defaults to 'value' for fundamental types)
        """
        # Special case for 'value' property in fundamental types
        if field_name == 'value' and hasattr(self, '_value'):
            # Reset the internal _value to None
            self._value = None
        elif hasattr(self, field_name):
            attr = getattr(self, field_name)
            # If the attribute is a CData (or HierarchicalObject), call destroy before deleting
            if isinstance(attr, HierarchicalObject):
                try:
                    attr.destroy()
                except Exception:
                    pass
            try:
                delattr(self, field_name)
            except AttributeError:
                # Can't delete properties without deleters - that's OK
                pass

        # Mark as not set
        self._value_states[field_name] = ValueState.NOT_SET

    def getValueState(self, field_name: str = "value") -> ValueState:
        """Get the current state of a field.

        Args:
            field_name: Name of the field

        Returns:
            ValueState indicating current state
        """
        return self._value_states.get(field_name, ValueState.NOT_SET)

    def isDefault(self, field_name: str = 'value') -> bool:
        """Check if a field is at its default value.

        Args:
            field_name: Name of the field to check (default: 'value')

        Returns:
            True if the field is at its default value, False otherwise
        """
        return self.getValueState(field_name) == ValueState.DEFAULT

    def setToDefault(self, field_name: str):
        """Set a field to its default value.

        Args:
            field_name: Name of the field to set to default
        """
        if field_name in self._default_values:
            # Set without triggering "explicitly set" state
            self._value_states[field_name] = ValueState.DEFAULT
            super().__setattr__(field_name, self._default_values[field_name])
        else:
            # No default available, unset it
            self.unSet(field_name)

    def setDefault(self, value: Any):
        """Set the default value for this object (old API compatibility).

        Args:
            value: The default value to set
        """
        # For simple value types, set the value and mark as DEFAULT
        if hasattr(self, 'value'):
            super().__setattr__('_value', value)
            self._value_states['value'] = ValueState.DEFAULT
            self._default_values['value'] = value
        else:
            # For complex types, store in default values
            self._default_values['default'] = value

    def getDefaultValue(self, field_name: str) -> Any:
        """Get the default value for a field.

        Args:
            field_name: Name of the field

        Returns:
            Default value or None if no default exists
        """
        return self._default_values.get(field_name)

    def validity(self):
        """Validate this object and return an error report.

        This method checks the object's state and validates it against
        qualifiers (min, max, enumerators, etc.). It also recursively
        validates all children in the hierarchy.

        If this object has allowUndefined=True and is not set, child validation
        errors are downgraded to warnings (they should not block execution).

        Returns:
            CErrorReport containing any validation errors/warnings
        """
        from .error_reporting import CErrorReport, SEVERITY_ERROR

        report = CErrorReport()

        # Check if this object is optional and not set
        # If so, child validation errors should be warnings, not errors
        # Note: allowUndefined=None (not explicitly set) is treated as True (optional)
        # Only allowUndefined=False explicitly marks a field as required
        allow_undefined = self.get_qualifier('allowUndefined')
        is_optional = allow_undefined is None or allow_undefined is True
        is_optional_and_unset = (
            is_optional and
            not self.isSet(allowUndefined=False, allowDefault=False)
        )

        # Base CData validation - can be extended by subclasses
        # Check if object has required qualifiers
        if hasattr(self, '_qualifiers') and self._qualifiers:
            # Basic validation is done here
            # Subclasses will add their own validation
            pass

        # Recursively validate all children in the hierarchy
        # This ensures that CData subclasses with child fields (like CImportUnmerged
        # with cell, wavelength, file, etc.) have their children validated
        if hasattr(self, 'children'):
            for child in self.children():
                if hasattr(child, 'validity'):
                    try:
                        child_report = child.validity()
                        if child_report:
                            # If this object is optional and not set, downgrade
                            # child errors to warnings (they should not block execution)
                            if is_optional_and_unset:
                                child_report.downgrade_to_warnings()
                            report.extend(child_report)
                    except Exception:
                        # Don't let one child's validation failure stop others
                        pass

        return report

    def dataOrder(self) -> list:
        """Return the order of child data items for serialization.

        This method provides the canonical ordering for child items, used by
        JSON serialization (CCP4i2JsonEncoder) and other operations that need
        consistent ordering.

        The ordering priority is:
        1. Explicit CONTENTS_ORDER or CONTENT_ORDER attribute (if defined)
        2. Auto-generated order from @cdata_class attributes metadata (MRO walk)
        3. Fallback to children() order (preserves addition order)

        For fundamental types (CInt, CFloat, CString, CBoolean), this returns
        an empty list since they have no child data items.

        Container types (CContainer, CList) may override this method for
        specialized ordering behavior.

        Returns:
            List of child names in serialization order
        """
        # Get all actual child names from children()
        all_children = set()
        for child in self.children():
            if hasattr(child, 'objectName'):
                name = child.objectName()
                if name:
                    all_children.add(name)

        # Also include children from _data_order that may not be in children() yet
        # This ensures defined but unset children are included in serialization
        if hasattr(self, '_data_order') and self._data_order:
            all_children.update(self._data_order)

        # Get preferred ordering from various sources
        preferred_order = []

        # 1. Check for explicit CONTENTS_ORDER property
        if hasattr(self, 'CONTENTS_ORDER'):
            try:
                val = self.CONTENTS_ORDER
                if isinstance(val, (list, tuple)) and val:
                    preferred_order = list(val)
            except Exception:
                pass

        # 2. Check for explicit CONTENT_ORDER attribute (if no CONTENTS_ORDER)
        if not preferred_order and hasattr(self, 'CONTENT_ORDER'):
            val = getattr(self, 'CONTENT_ORDER', None)
            if isinstance(val, (list, tuple)) and val:
                preferred_order = list(val)

        # 3. Auto-generate from @cdata_class attributes metadata (if no explicit order)
        if not preferred_order:
            for cls in reversed(type(self).__mro__):
                if cls is object:
                    continue
                metadata = getattr(cls, "_metadata", None)
                if metadata and hasattr(metadata, 'attributes'):
                    for attr_name in metadata.attributes.keys():
                        if attr_name not in preferred_order:
                            preferred_order.append(attr_name)

        # Build complete ordering: preferred order first (filtering to actual children),
        # then remaining children not in preferred order
        result = []
        seen = set()

        # Add items from preferred_order that exist as children
        for name in preferred_order:
            if name in all_children and name not in seen:
                result.append(name)
                seen.add(name)

        # Add remaining children not in preferred_order
        for child in self.children():
            if hasattr(child, 'objectName'):
                name = child.objectName()
                if name and name not in seen:
                    result.append(name)
                    seen.add(name)

        return result

    def _check_all_registered_attributes_set(self) -> bool:
        """Check if all registered attributes (from @cdata_class) are set.

        This method inspects the class metadata to find all attributes
        defined in the @cdata_class decorator and checks if each one
        has been explicitly set.

        Returns:
            True if all registered attributes are set, False otherwise.
            Returns True if no metadata is found (backward compatibility).
        """
        # Get class metadata
        metadata = getattr(self.__class__, '_metadata', None)
        if metadata is None:
            # No metadata - assume all set for backward compatibility
            return True

        # Get registered attributes from metadata
        registered_attrs = metadata.attributes if hasattr(metadata, 'attributes') else {}
        if not registered_attrs:
            # No registered attributes - assume all set
            return True

        # Check each registered attribute
        for attr_name in registered_attrs.keys():
            if not hasattr(self, attr_name):
                # Attribute doesn't exist - not set
                return False

            attr_value = getattr(self, attr_name)
            if isinstance(attr_value, CData):
                # For CData children, check if they're set
                if not attr_value.isSet(allowDefault=False):
                    return False
            else:
                # For non-CData attributes, check value state
                state = self._value_states.get(attr_name, ValueState.NOT_SET)
                if state == ValueState.NOT_SET:
                    return False

        return True

    def getEtree(self, name: str = None, excludeUnset: bool = False, allSet: bool = False):
        """Serialize this object to an XML ElementTree element.

        Args:
            name: Optional name for the XML element (uses self.objectName() if not provided)
            excludeUnset: If True, only serialize fields that have been explicitly set
            allSet: If True, only serialize if ALL registered attributes are set.
                   This checks all attributes defined in the @cdata_class decorator.

        Returns:
            xml.etree.ElementTree.Element representing this object, or empty element
            if allSet=True and not all attributes are set
        """
        import xml.etree.ElementTree as ET

        # Use objectName() to get hierarchical object name, not CData 'name' attribute
        element_name = name if name is not None else (self.objectName() if hasattr(self, 'objectName') else 'data')
        elem = ET.Element(element_name)

        # If allSet is True, check if ALL registered attributes are set
        if allSet:
            if not self._check_all_registered_attributes_set():
                return elem  # Return empty element

        # For simple value types, store the value as text
        if hasattr(self, 'value'):
            value = getattr(self, 'value')

            if value is not None:
                # Check if we should exclude this based on set state
                if excludeUnset and hasattr(self, 'isSet'):
                    # Use allowDefault=False to exclude default values
                    if not self.isSet('value', allowDefault=False):
                        return elem  # Return empty element for unset values
                elem.text = str(value)

        # For containers and complex types, serialize hierarchical children
        # Use dataOrder() for canonical ordering - it returns a complete list
        # of all children in the correct serialization order
        children_by_name = {}
        for child in self.children():
            if isinstance(child, CData):
                child_name = child.objectName() if hasattr(child, 'objectName') else None
                if child_name:
                    children_by_name[child_name] = child

        # dataOrder() returns complete ordering of all children
        ordered_names = self.dataOrder() if hasattr(self, 'dataOrder') else []

        for attr_name in ordered_names:
            if attr_name not in children_by_name:
                continue
            child = children_by_name[attr_name]

            # Check if child should be excluded
            if excludeUnset and hasattr(child, 'isSet'):
                # For CData objects, check if they have any set values
                if not child.isSet(allowDefault=False, allSet=False):
                    continue  # Skip unset children

            # Pass allSet and excludeUnset to children
            child_elem = child.getEtree(attr_name, excludeUnset=excludeUnset, allSet=allSet)

            # Only append if the child element has content
            if excludeUnset or allSet:
                has_content = child_elem.text or len(child_elem) > 0
                if has_content:
                    elem.append(child_elem)
            else:
                elem.append(child_elem)

        return elem

    def setEtree(self, element, ignore_missing: bool = False, preserve_state: bool = False):
        """Deserialize from an XML ElementTree element.

        Args:
            element: xml.etree.ElementTree.Element to deserialize from
            ignore_missing: If True, ignore attributes in XML that don't exist on this object
            preserve_state: If True, don't mark values as EXPLICITLY_SET when deserializing.
                           This is useful for copyData() to preserve the original value states.
        """
        # For simple value types, read from text
        if hasattr(self, 'value') and element.text:
            # Temporarily disable validation during deserialization
            # This prevents validation errors when loading values from XML that may
            # not meet current min/max constraints but were valid when saved
            old_skip_validation = getattr(self, '_skip_validation', False)
            self._skip_validation = True

            # Save old state if preserving
            old_state = None
            if preserve_state and hasattr(self, '_value_states'):
                old_state = self._value_states.get('value', ValueState.NOT_SET)

            try:
                # Determine the type and convert
                if hasattr(self, '__class__'):
                    class_name = self.__class__.__name__
                    if class_name == 'CInt':
                        self.value = int(element.text)
                    elif class_name == 'CFloat':
                        self.value = float(element.text)
                    elif class_name == 'CBoolean':
                        self.value = element.text.lower() in ('true', '1', 'yes')
                    elif class_name == 'CString':
                        self.value = element.text
                    else:
                        self.value = element.text
            finally:
                # Restore previous validation state
                self._skip_validation = old_skip_validation

                # Restore previous value state if preserving
                if preserve_state and old_state is not None and hasattr(self, '_value_states'):
                    self._value_states['value'] = old_state

        # For containers, deserialize children
        for child in element:
            child_name = child.tag
            if hasattr(self, child_name):
                child_obj = getattr(self, child_name)
                if isinstance(child_obj, CData):
                    child_obj.setEtree(child, ignore_missing=ignore_missing, preserve_state=preserve_state)
            elif not ignore_missing:
                raise AttributeError(f"Object has no attribute '{child_name}'")

    def getQualifiersEtree(self):
        """Serialize qualifiers to an XML ElementTree element.

        Returns:
            xml.etree.ElementTree.Element containing qualifiers
        """
        import xml.etree.ElementTree as ET

        qualifiers_elem = ET.Element('qualifiers')

        if hasattr(self, '_qualifiers') and self._qualifiers:
            for key, value in self._qualifiers.items():
                qual_elem = ET.Element(key)
                if value is not None:
                    if isinstance(value, bool):
                        qual_elem.text = 'true' if value else 'false'
                    elif isinstance(value, list):
                        qual_elem.text = ','.join(str(v) for v in value)
                    else:
                        qual_elem.text = str(value)
                qualifiers_elem.append(qual_elem)

        return qualifiers_elem

    def setQualifiersEtree(self, qualifiers_element):
        """Deserialize qualifiers from an XML ElementTree element.

        Args:
            qualifiers_element: xml.etree.ElementTree.Element containing qualifiers
        """
        if qualifiers_element is None:
            return

        for child in qualifiers_element:
            key = child.tag
            text = child.text

            # Parse the value based on content
            if text is None:
                value = None
            elif text.lower() in ('true', 'false'):
                value = text.lower() == 'true'
            elif ',' in text:
                # List of values
                value = [item.strip() for item in text.split(',')]
            else:
                # Try to parse as number, otherwise keep as string
                try:
                    if '.' in text:
                        value = float(text)
                    else:
                        value = int(text)
                except (ValueError, AttributeError):
                    value = text

            self.set_qualifier(key, value)

    def _is_value_type(self) -> bool:
        """Check if this is a simple value type (like CString, CInt, etc.)."""
        # Simple heuristic: if class has only basic Python types as attributes
        # and no complex CData children, it's likely a value type
        value_type_patterns = ["String", "Int", "Float", "Bool", "OneWord"]
        return any(
            pattern in self.__class__.__name__ for pattern in value_type_patterns
        )

    def _smart_assign_from_dict(self, source_dict: dict):
        """Assign values from dictionary to object attributes."""
        for key, value in source_dict.items():
            setattr(self, key, value)

    def _smart_assign_from_cdata(self, source: "CData"):
        """Handle smart assignment from another CData object.

        For containers, this merges children from source into self.
        For value types, this copies the value AND the value state.
        """
        if self._is_value_type() and source._is_value_type():
            # Value assignment: copy the underlying value AND preserve value state
            # For simple types, copy their primary value attribute
            primary_attrs = ["value", "text", "string", "content"]
            for attr in primary_attrs:
                if hasattr(source, attr):
                    # Set the value
                    setattr(self, attr, getattr(source, attr))
                    # IMPORTANT: Copy the value state from source AFTER setting value
                    # This overwrites the EXPLICITLY_SET state that setattr() created
                    if hasattr(source, '_value_states') and attr in source._value_states:
                        if hasattr(self, '_value_states'):
                            self._value_states[attr] = source._value_states[attr]
                    return

            # If no primary attribute found, copy hierarchical children from source
            # Use children() to get only actual CData children, not arbitrary __dict__ entries
            for child in source.children():
                if isinstance(child, CData):
                    key = child.objectName() if hasattr(child, 'objectName') else None
                    if key:
                        setattr(self, key, child)
        else:
            # Complex type assignment (like containers and lists)
            # For CList types, copy items from source to self
            from .fundamental_types import CList
            if isinstance(self, CList) and isinstance(source, CList):
                # Clear existing items
                self.clear()
                # Copy items from source to self
                for item in source:
                    self.append(item)
                return

            # For containers, merge attributes from source instead of replacing
            from .ccontainer import CContainer
            if isinstance(self, CContainer) and isinstance(source, CContainer):
                # Container merging: iterate over source's hierarchical children
                # Use children() to get only actual CData children, not all dir() entries
                for child in source.children():
                    # Only copy CData objects
                    if not isinstance(child, CData):
                        continue

                    # Get the child's name within the source
                    key = child.objectName() if hasattr(child, 'objectName') else None
                    if not key:
                        continue

                    # Check if we already have an attribute with this name
                    if hasattr(self, key):
                        # Attribute already exists - recursively merge if both are containers
                        existing_child = getattr(self, key)
                        if isinstance(existing_child, CContainer) and isinstance(child, CContainer):
                            existing_child._smart_assign_from_cdata(child)
                        else:
                            # Not containers - replace with new one
                            setattr(self, key, child)
                    else:
                        # New attribute - add it
                        # Need to reparent the child from source to self
                        child.set_parent(self)
                        setattr(self, key, child)
            else:
                # Non-container complex type: copy hierarchical children, but only if explicitly set
                # This prevents NOT_SET fields from being marked as EXPLICITLY_SET
                # Use children() to get only actual hierarchical children, not all __dict__ entries

                # For CDataFile types, preserve existing path-related fields (relPath, baseName, project)
                # when the destination already has them explicitly set. This prevents subjob paths
                # from overwriting pre-configured main job destination paths during RegisterSubOutputAsMain.
                # Path fields define WHERE a file should be located, while content fields (contentFlag,
                # subType, annotation) define WHAT the file contains.
                #
                # IMPORTANT: dbFileId is a special case - it should NEVER be copied from source to dest.
                # The dbFileId points to a specific database record for the SOURCE file. If we copy it,
                # getFullPath() will use it to look up the path in the database, which returns the
                # SOURCE file's path, not the DESTINATION file's intended path. The destination should
                # get its own dbFileId when it's registered in the database with its own path.
                from .cdata_file import CDataFile
                is_file_type = isinstance(self, CDataFile)
                # Fields to preserve if destination has them set
                path_fields_preserve = {'relPath', 'baseName', 'project'}
                # Fields to NEVER copy (database identity fields that must be unique per file)
                path_fields_never_copy = {'dbFileId'}

                for child in source.children():
                    # Get the child's name within the source
                    key = child.objectName() if hasattr(child, 'objectName') else None
                    if not key:
                        continue

                    # For CDataFile, NEVER copy dbFileId - it's the source's database identity
                    if is_file_type and key in path_fields_never_copy:
                        continue

                    # For CDataFile, preserve existing path fields if already explicitly set
                    # This allows content to be copied without overwriting the destination path
                    if is_file_type and key in path_fields_preserve:
                        existing = getattr(self, key, None)
                        if existing is not None and hasattr(existing, 'isSet') and existing.isSet(allowDefault=False):
                            # Destination already has this path field set - preserve it
                            continue

                    # For CData attributes, only copy if explicitly set
                    if isinstance(child, CData) and hasattr(child, 'isSet'):
                        if child.isSet(allowDefault=False):
                            setattr(self, key, child)

    def _deep_copy_from(self, source: "CData"):
        """Deep copy values and value states from source, without creating references.

        This method performs a true deep copy:
        - For simple value types (CInt, CFloat, CString, CBoolean): copies the value and value state
        - For containers/complex types: recursively deep copies all children

        Unlike _smart_assign_from_cdata, this never creates references to source objects.
        This is intended for use by copyData() to ensure independent copies.

        Args:
            source: The source CData object to copy from
        """
        from .fundamental_types import CList, CString, CInt, CFloat, CBoolean

        # Handle simple value types - copy value and state directly
        if isinstance(self, (CString, CInt, CFloat, CBoolean)) and isinstance(source, (CString, CInt, CFloat, CBoolean)):
            # Copy the value
            if hasattr(source, '_value'):
                self._value = source._value
            elif hasattr(source, 'value'):
                # Use direct attribute access, not property, to avoid triggering setters
                object.__setattr__(self, '_value', getattr(source, '_value', source.value))

            # Copy the value state
            if hasattr(source, '_value_states') and hasattr(self, '_value_states'):
                self._value_states = source._value_states.copy()
            return

        # Handle CList - deep copy each item
        if isinstance(self, CList) and isinstance(source, CList):
            self.clear()
            for source_item in source:
                # Create a new item of the SAME TYPE as the source item
                # This is critical because makeItem() uses the destination's subItem
                # qualifier which may default to CString, but the source may contain
                # CMtzColumn or other specialized types
                source_item_type = type(source_item)
                try:
                    new_item = source_item_type(parent=self, name="[?]")
                except Exception:
                    # Fallback to makeItem() if source type can't be instantiated
                    new_item = self.makeItem()

                if hasattr(new_item, '_deep_copy_from') and isinstance(source_item, CData):
                    new_item._deep_copy_from(source_item)
                elif hasattr(source_item, 'value'):
                    # Simple value type without _deep_copy_from
                    new_item.value = source_item.value
                    if hasattr(source_item, '_value_states') and hasattr(new_item, '_value_states'):
                        new_item._value_states = source_item._value_states.copy()
                self.append(new_item)
            return

        # Handle containers/complex types - recursively copy children
        for source_child in source.children():
            if not isinstance(source_child, CData):
                continue

            child_name = source_child.objectName() if hasattr(source_child, 'objectName') else None
            if not child_name:
                continue

            # Get or create the destination child
            if hasattr(self, child_name):
                dest_child = getattr(self, child_name)
                if isinstance(dest_child, CData) and isinstance(source_child, CData):
                    # Recursively deep copy
                    dest_child._deep_copy_from(source_child)
            # Note: We don't create new children that don't exist in destination
            # copyData's dataList controls which items to copy

    def _setup_hierarchy_for_value(self, key: str, value: Any):
        """Set up hierarchical relationships for attribute values.

        Only sets parent if not already set, respecting explicit parent assignments.
        Only sets name if not already set.

        When an attribute name differs from the child's _name (e.g., self.pointless = plugin
        where plugin._name = 'aimless_pipe_1'), this also adds an alias in the parent's
        _children_by_name and _child_storage so the child can be found by the attribute name.

        Also handles the case where a CData with an existing parent is assigned as a reference:
        e.g., self.AUTOCUTOFF = self.container.controlParameters.AUTOCUTOFF
        The original parent is kept, but we add a reference so it can be looked up by the key.
        """
        if isinstance(value, CData):
            existing_parent = value.parent()
            # Only set parent if not already set (respect explicit parent assignment)
            if existing_parent is None:
                value.set_parent(self)
            # Only set hierarchical name if not already set
            if not value._name:
                value._name = key

            # Add a lookup reference if:
            # 1. The attribute name differs from the child's _name (alias case), OR
            # 2. The object has a different parent (reference/shortcut case)
            # This handles both:
            # - self.pointless = self.makePluginObject('pointless') where plugin._name differs
            # - self.AUTOCUTOFF = self.container.controlParameters.AUTOCUTOFF (same name, different parent)
            if value._name != key or (existing_parent is not None and existing_parent is not self):
                import weakref
                child_ref = weakref.ref(value)
                self._children_by_name[key] = child_ref
                self._child_storage[key] = value
        elif isinstance(value, list):
            # Handle list of CData objects
            for i, item in enumerate(value):
                if isinstance(item, CData):
                    # Only set parent if not already set
                    if item.parent() is None:
                        item.set_parent(self)
                    if not item._name:
                        item._name = f"{key}[{i}]"

    def __setattr__(self, name: str, value: Any):
        """Override setattr to handle smart assignment and hierarchical relationships."""

        # Allow setting internal attributes normally during initialization
        # Note: 'name' is NOT in this list - it's a regular CData attribute that uses smart assignment
        # HierarchicalObject's hierarchical name is stored separately in '_name'
        if (
            name.startswith("_")
            or name in ["parent", "children", "signals"]
            or not hasattr(self, "_hierarchy_initialized")
        ):
            super().__setattr__(name, value)
            return

        # Special handling for metadata attributes: assign directly if dict or list
        if name in ["qualifiers", "qualifiers_order", "qualifiers_definition", "CONTENT_ORDER", "subitem"]:
            if isinstance(value, (dict, list)):
                object.__setattr__(self, name, value)
                return

        # Handle smart assignment patterns
        # Since HierarchicalObject uses _name (not name), there's no collision with CData 'name' attributes
        existing_attr = getattr(self, name, None)

        # DEBUG: Print all contentFlag assignments
        if name == 'contentFlag':
            import traceback
            logger.debug(
                "[SETATTR] %s.%s = %s (type: %s)\n  existing: %s (type: %s)\n  Stack: %s",
                self.__class__.__name__, name, value, type(value).__name__,
                existing_attr, type(existing_attr).__name__,
                ''.join(traceback.format_stack()[-4:-1])
            )

        if isinstance(value, dict):
            # Dictionary assignment: update object attributes from dictionary
            if existing_attr is None:
                # Create new CData object and populate from dict
                new_obj = CData(name=name)
                new_obj._smart_assign_from_dict(value)
                self._setup_hierarchy_for_value(name, new_obj)
                super().__setattr__(name, new_obj)
                return
            elif isinstance(existing_attr, CData):
                # Update existing CData object from dictionary
                existing_attr._smart_assign_from_dict(value)
                # Mark as explicitly set since we're assigning new values
                if hasattr(self, "_value_states"):
                    self._value_states[name] = ValueState.EXPLICITLY_SET
                return  # Don't replace the object, just update it

        elif isinstance(value, CData) and isinstance(existing_attr, CData):
            # CData to CData assignment: use smart assignment logic
            # Only perform assignment if the source value is actually set
            # This prevents copying unset/default values which would mark them as explicitly set
            #
            # For container objects (like CRefinementPerformance), we need to check if ANY
            # child attribute is set, not just the parent's isSet() which checks the 'value'
            # attribute. The _has_any_set_children() helper handles this case.
            if hasattr(value, 'isSet') and not value.isSet():
                # Parent's isSet() returns False - but check if any children are set
                if not value._has_any_set_children():
                    # Source value is not set AND no children are set - skip the assignment
                    # This prevents legacy code like: obj.FRAC = source.FREER_FRACTION
                    # from marking FRAC as set when FREER_FRACTION is just a default
                    return

            # IMPORTANT: Detect plugin reassignment (different instances)
            # When SubstituteLigand does: self.aimlessPlugin = self.makePluginObject('aimless_pipe')
            # twice, we want to REPLACE the reference, not update the old instance in-place.
            # Otherwise the old instance's state (_childJobCounter=5) gets preserved.
            #
            # Check if both are CPluginScript instances (or subclasses) using isinstance().
            # Use deferred import to avoid circular dependency since CCP4PluginScript imports from base_object.
            is_plugin_reassignment = False
            if id(existing_attr) != id(value):
                try:
                    from ccp4i2.core.CCP4PluginScript import CPluginScript
                    is_plugin_reassignment = isinstance(existing_attr, CPluginScript) and isinstance(value, CPluginScript)
                except ImportError:
                    # CPluginScript not yet loaded, fall back to smart assignment
                    pass

            if is_plugin_reassignment:
                # Different plugin instances - bypass smart assignment and replace the reference
                # Fall through to normal assignment logic below
                pass
            else:
                # Same instance or non-plugin objects - use smart assignment
                existing_attr._smart_assign_from_cdata(value)
                # Mark as explicitly set since we're assigning a new value
                if hasattr(self, "_value_states"):
                    self._value_states[name] = ValueState.EXPLICITLY_SET
                return  # Don't replace the object, just update it

        elif (
            hasattr(existing_attr, 'setFullPath')
            and callable(getattr(existing_attr, 'setFullPath'))
            and isinstance(value, str)
        ):
            # Special case: Assigning a string path to an existing CDataFile
            # e.g., task.container.inputData.HKLIN = "/path/to/file.mtz"
            # Use duck typing to avoid circular import
            existing_attr.setFullPath(value)
            # Mark as explicitly set
            if hasattr(self, "_value_states"):
                self._value_states[name] = ValueState.EXPLICITLY_SET
            return  # Don't replace the object, just update its path

        elif (
            isinstance(value, list)
            and existing_attr is not None
            and hasattr(existing_attr, 'set')
            and hasattr(existing_attr, '_items')  # Duck-type check for CList
        ):
            # Special case: Assigning a list to an existing CList
            # e.g., columnGroup.columnList = [dict1, dict2, ...]
            # Use CList.set() to properly handle dict-to-object conversion via subItem
            existing_attr.set(value)
            # Mark as explicitly set
            if hasattr(self, "_value_states"):
                self._value_states[name] = ValueState.EXPLICITLY_SET
            return  # Don't replace the object, just update its items

        elif (
            existing_attr is not None and isinstance(existing_attr, CData)
            and existing_attr._is_value_type()
        ):
            if name in ['contentFlag', 'subType']:  # DEBUG
                logger.debug("Branch: Value type smart assign for %s", name)  # DEBUG
            # Primitive value assignment to existing CData value type
            # e.g., ctrl.NCYCLES = 25 where ctrl.NCYCLES is a CInt
            if isinstance(value, (int, float, str, bool)):
                # Check type compatibility and perform conversions if needed
                type_compatible = False
                converted_value = value
                if hasattr(existing_attr, "value"):
                    # Import types locally to avoid circular import
                    from .fundamental_types import CInt, CFloat, CBoolean
                    try:
                        from .fundamental_types import CString
                    except ImportError:
                        CString = None

                    # Handle CInt: accept int or string convertible to int
                    if isinstance(existing_attr, CInt):
                        if isinstance(value, int):
                            type_compatible = True
                        elif isinstance(value, str):
                            # Legacy code often uses str(int_value) - convert back
                            try:
                                converted_value = int(value)
                                type_compatible = True
                            except (ValueError, TypeError):
                                pass

                    # Handle CFloat: accept float, int, or string convertible to float
                    elif isinstance(existing_attr, CFloat):
                        if isinstance(value, (int, float)):
                            type_compatible = True
                        elif isinstance(value, str):
                            try:
                                converted_value = float(value)
                                type_compatible = True
                            except (ValueError, TypeError):
                                pass

                    # Handle CBoolean: accept bool or string convertible to bool
                    elif isinstance(existing_attr, CBoolean):
                        if isinstance(value, bool):
                            type_compatible = True
                        elif isinstance(value, str):
                            # Handle common boolean string representations
                            if value.lower() in ('true', 'yes', '1'):
                                converted_value = True
                                type_compatible = True
                            elif value.lower() in ('false', 'no', '0'):
                                converted_value = False
                                type_compatible = True

                    # Handle CString: always accept strings
                    elif CString is not None and isinstance(existing_attr, CString) and isinstance(value, str):
                        type_compatible = True

                if type_compatible:
                    # Update the value attribute of the existing CData object
                    if name in ['contentFlag', 'subType']:  # DEBUG
                        logger.debug("[SMART ASSIGN] Updating %s.value = %s, keeping %s", name, converted_value, type(existing_attr).__name__)  # DEBUG
                    existing_attr.value = converted_value
                    # Mark as explicitly set (no longer needed with delegation, but keeping for non-CData attributes)
                    if hasattr(self, "_value_states"):
                        self._value_states[name] = ValueState.EXPLICITLY_SET
                    return  # Don't replace the object, just update its value

        # For new attributes or non-smart assignment, handle hierarchy
        self._setup_hierarchy_for_value(name, value)

        # IMPORTANT: CData children are stored ONLY in hierarchy (via set_parent),
        # NOT in __dict__. This keeps hierarchy as the single source of truth.
        # The __getattribute__ override provides O(1) access via _children_by_name.
        # Non-CData values (primitives, lists of primitives) are stored in __dict__.
        if isinstance(value, CData):
            # CData children: already added to hierarchy via _setup_hierarchy_for_value
            # Don't store in __dict__ - __getattribute__ will find it via _children_by_name
            pass
        else:
            # Non-CData values: store normally in __dict__
            super().__setattr__(name, value)

        # Track that this value has been explicitly set (unless it's internal)
        # IMPORTANT: Don't mark 'value' attribute here for types with @property setters
        # (CInt, CFloat, CBoolean) because those setters handle state tracking themselves.
        # This prevents parameters without defaults from being incorrectly marked as EXPLICITLY_SET
        # during .def.xml loading and container merging operations.
        # However, DO track for CString which doesn't use a @property setter.
        from .fundamental_types import CInt, CFloat, CBoolean
        has_value_property = isinstance(self, (CInt, CFloat, CBoolean))
        skip_value_tracking = (name == "value" and has_value_property)

        # Exclusion list for value state tracking
        # Note: 'name' is NOT excluded - HierarchicalObject uses _name, so no collision
        exclusion_list = ["parent", "children", "signals"]

        if (hasattr(self, "_value_states") and not name.startswith("_")
            and name not in exclusion_list
            and not skip_value_tracking):
            self._value_states[name] = ValueState.EXPLICITLY_SET

    def __getattr__(self, name: str):
        """Auto-create metadata-defined attributes when accessed.

        This is called when normal attribute lookup fails (after __getattribute__).
        If the attribute is defined in metadata but hasn't been instantiated yet,
        create it and add to hierarchy via set_parent().
        """
        # Avoid infinite recursion - only process if object is fully initialized
        if name.startswith("_") or not hasattr(self, "_hierarchy_initialized"):
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

        # Check if this attribute is defined in metadata
        metadata = getattr(self.__class__, "_metadata", None)
        if metadata and name in metadata.attributes:
            # Auto-create the attribute from metadata
            from .class_metadata import MetadataAttributeFactory
            attr_def = metadata.attributes[name]
            attr_obj = MetadataAttributeFactory.create_attribute(name, attr_def, self)

            # Add to hierarchy via set_parent (populates _children_by_name for O(1) access)
            # The __getattribute__ override will find it via _children_by_name
            # DON'T store in __dict__ - hierarchy is the single source of truth for CData
            if isinstance(attr_obj, CData):
                attr_obj._name = name  # Set hierarchical name
                attr_obj.set_parent(self)  # This adds to _children_by_name
            else:
                # Non-CData (shouldn't happen for metadata attributes, but handle gracefully)
                self.__dict__[name] = attr_obj

            # Mark as NOT_SET initially (will be marked EXPLICITLY_SET when assigned)
            if hasattr(self, "_value_states"):
                self._value_states[name] = ValueState.NOT_SET

            return attr_obj

        # Not a metadata attribute - raise AttributeError
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def __delattr__(self, name: str):
        """Override delattr to handle hierarchy-stored children.

        CData children are stored in the hierarchy (_children_by_name, _child_storage)
        rather than __dict__. This method handles deletion from both locations.

        Note: This silently succeeds if the attribute doesn't exist, to support
        defensive deletion patterns in legacy code (e.g., container.deleteObject('MAYBE_EXISTS')).
        """
        found = False

        # First try to remove from hierarchy (for CData children)
        if name in self._children_by_name:
            child_ref = self._children_by_name.pop(name, None)
            self._child_storage.pop(name, None)
            found = True

            # If there's a live child, clean it up
            if child_ref is not None:
                child = child_ref()
                if child is not None and child.parent() is self:
                    # Remove from parent's _children set as well
                    self._remove_child(child)
                    # Note: don't call child.destroy() here - deleteObject already does that

        # Also try to remove from __dict__ (for non-CData attributes)
        if name in self.__dict__:
            del self.__dict__[name]
            found = True

        # Silently succeed if not found - this supports defensive delete patterns
        # Legacy code often does: container.deleteObject('SOMETHING') without checking first
        # Note: We only raise if the attribute DEFINITELY should exist (internal names)
        if not found and name.startswith('_'):
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def find_children_by_type(self, target_type):
        """Find all descendant objects of a specific type.

        Replaces server/ccp4i2/lib/cdata_utils.find_objects_by_type() as a core method.

        Args:
            target_type: Type to search for (e.g., CPerformanceIndicator, CInt)

        Returns:
            List of objects matching the target type

        Example:
            >>> from ccp4i2.core.CCP4XtalData import CPerformanceIndicator
            >>> kpis = plugin.outputData.find_children_by_type(CPerformanceIndicator)
            >>> for kpi in kpis:
            ...     print(f"KPI: {kpi.object_path()}")
        """
        objects = []
        visited = set()  # Prevent cycles

        def traverse(obj):
            obj_id = id(obj)
            if obj_id in visited:
                return
            visited.add(obj_id)

            if isinstance(obj, target_type):
                objects.append(obj)

            # Use children() method from HierarchicalObject for traversal
            if hasattr(obj, 'children'):
                try:
                    for child in obj.children():
                        if child is not None:
                            traverse(child)
                except Exception:
                    pass

        traverse(self)
        return objects

    def find_children_matching(self, predicate):
        """Find all descendant objects matching a predicate function.

        Replaces server/ccp4i2/lib/cdata_utils.find_objects_matching() as a core method.

        Args:
            predicate: Function that takes an object and returns True for matches

        Returns:
            List of objects matching the predicate

        Example:
            >>> # Find all set data files
            >>> from .cdata_file import CDataFile
            >>> set_files = plugin.outputData.find_children_matching(
            ...     lambda obj: isinstance(obj, CDataFile) and obj.isSet()
            ... )
        """
        matches = []
        visited = set()

        def traverse(obj):
            obj_id = id(obj)
            if obj_id in visited:
                return
            visited.add(obj_id)

            try:
                if predicate(obj):
                    matches.append(obj)
            except Exception:
                pass

            # Use children() method for traversal
            if hasattr(obj, 'children'):
                try:
                    for child in obj.children():
                        if child is not None:
                            traverse(child)
                except Exception:
                    pass

        traverse(self)
        return matches

    def __bool__(self):
        """Return True if this parameter has a value set (explicitly or via default).

        The semantic meaning of `if param:` is "does this parameter have a value?"
        rather than checking the truthiness of the underlying value. This matches
        the common usage pattern in plugin code:
            if plugin.controlParameters.SOMEPARAM:
                # Use the parameter's value

        The query above intends to check whether SOMEPARAM has been configured (is set),
        not whether its value is truthy for the underlying Python type.

        For code that needs to check the actual value's truthiness, convert to the
        native type first: `if bool(int(param))` or `if float(param) != 0:`

        For code that needs to check explicit vs default: use `isSet(allowDefault=False)`

        Note: CBoolean overrides this to return both set AND truthy, since boolean
        parameters typically want to check both presence and value.

        Returns:
            True if the value is set (explicitly or via default), False otherwise
        """
        return self.isSet(allowDefault=True)

