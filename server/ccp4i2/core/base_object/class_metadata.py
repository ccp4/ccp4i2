"""Decorator-based metadata system for CData classes.

This module provides decorators and utilities to embed metadata directly
in class definitions, making the system more maintainable than external JSON files.
"""

from typing import Dict, Any, Optional, Type, List
from dataclasses import dataclass, field
from enum import Enum


class AttributeType(Enum):
    """Types of attributes that can be created.

    Only fundamental types (CInt, CFloat, CBoolean, CString) have their own enum values.
    All other types (CFilePath, CUUID, CList, etc.) should use CUSTOM with custom_class parameter.
    """

    INT = "CInt"
    FLOAT = "CFloat"
    BOOLEAN = "CBoolean"
    STRING = "CString"
    CUSTOM = "Custom"


@dataclass
class AttributeDefinition:
    """Definition of a class attribute.

    All attribute constraints (min/max/default/enumerators/etc.) should be
    defined in the class-level qualifiers, not here. This class only defines
    what TYPE of attribute to create.
    """

    attr_type: AttributeType
    custom_class: Optional[str] = None  # For AttributeType.CUSTOM types


@dataclass
class ClassMetadata:
    """Complete metadata for a CData class."""

    attributes: Dict[str, AttributeDefinition] = field(default_factory=dict)
    qualifiers: Dict[str, Any] = field(default_factory=dict)
    error_codes: Dict[int, str] = field(default_factory=dict)
    docstring: Optional[str] = None
    file_extensions: Optional[List[str]] = None
    mime_type: Optional[str] = None
    gui_label: Optional[str] = None
    contents_order: Optional[List[str]] = None
    qualifiers_order: Optional[List[str]] = None
    qualifiers_definition: Optional[Dict[str, Any]] = None
    content_qualifiers: Optional[Dict[str, Dict[str, Any]]] = None  # Per-field qualifiers


# Global registry of class metadata
_CLASS_METADATA_REGISTRY: Dict[str, ClassMetadata] = {}


def attribute(attr_type: AttributeType, custom_class: Optional[str] = None) -> AttributeDefinition:
    """Helper function to create attribute definitions.

    All attribute constraints (min/max/default/enumerators/etc.) should be
    defined in class-level qualifiers, not here.

    Args:
        attr_type: The type of attribute (INT, FLOAT, STRING, CUSTOM, etc.)
        custom_class: For AttributeType.CUSTOM, the class name to instantiate

    Returns:
        AttributeDefinition instance

    Example:
        project = attribute(AttributeType.CUSTOM, custom_class="CProjectId")
        label = attribute(AttributeType.CUSTOM, custom_class="COneWord")
        items = attribute(AttributeType.CUSTOM, custom_class="CList")
    """
    return AttributeDefinition(attr_type=attr_type, custom_class=custom_class)


def cdata_class(
    attributes: Optional[Dict[str, AttributeDefinition]] = None,
    qualifiers: Optional[Dict[str, Any]] = None,
    error_codes: Optional[Dict[int, str]] = None,
    file_extensions: Optional[List[str]] = None,
    mime_type: Optional[str] = None,
    gui_label: Optional[str] = None,
    contents_order: Optional[List[str]] = None,
    qualifiers_order: Optional[List[str]] = None,
    qualifiers_definition: Optional[Dict[str, Any]] = None,
    content_qualifiers: Optional[Dict[str, Dict[str, Any]]] = None,
):
    """Class decorator to add metadata to CData classes.

    Args:
        attributes: Dictionary of attribute name -> AttributeDefinition
        qualifiers: Dictionary of class qualifiers
        error_codes: Dictionary of error code -> message
        file_extensions: List of supported file extensions
        mime_type: MIME type for file classes
        gui_label: Label for GUI display
        contents_order: List specifying display order of attributes in UI
        qualifiers_order: List specifying display order of qualifiers
        qualifiers_definition: Dictionary of qualifier type definitions
        content_qualifiers: Per-field qualifiers for child attributes (from CONTENTS)

    Example:
        @cdata_class(
            attributes={
                'project': attribute(AttributeType.CUSTOM, custom_class="CProjectId"),
                'baseName': attribute(AttributeType.CUSTOM, custom_class="CFilePath"),
                'size': attribute(AttributeType.INT, default=0, min_value=0)
            },
            file_extensions=['dat', 'txt'],
            mime_type='text/plain'
        )
        class CDataFile(CData):
            '''A data file with embedded metadata.'''
            pass
    """

    def decorator(cls: Type) -> Type:
        # Create metadata object
        metadata = ClassMetadata(
            attributes=attributes or {},
            qualifiers=qualifiers or {},
            error_codes=error_codes or {},
            docstring=cls.__doc__,
            file_extensions=file_extensions,
            mime_type=mime_type,
            gui_label=gui_label,
            contents_order=contents_order,
            qualifiers_order=qualifiers_order,
            qualifiers_definition=qualifiers_definition,
            content_qualifiers=content_qualifiers,
        )

        # Store in global registry
        _CLASS_METADATA_REGISTRY[cls.__name__] = metadata

        # Store as class attribute for easy access
        cls._metadata = metadata

        # ALSO set direct class attributes for backward compatibility with CData.__init__
        # Use _class_* naming to avoid shadowing instance methods like qualifiers()
        if qualifiers:
            cls._class_qualifiers = qualifiers
        if qualifiers_order:
            cls.qualifiers_order = qualifiers_order
        if qualifiers_definition:
            cls.qualifiers_definition = qualifiers_definition
        if contents_order:
            cls.CONTENT_ORDER = contents_order
        if error_codes:
            cls.ERROR_CODES = error_codes

        return cls

    return decorator


def get_class_metadata(class_name: str) -> Optional[ClassMetadata]:
    """Get metadata for a class by name.

    Args:
        class_name: Name of the class

    Returns:
        ClassMetadata instance or None if not found
    """
    return _CLASS_METADATA_REGISTRY.get(class_name)


def get_class_metadata_by_type(cls: Type) -> Optional[ClassMetadata]:
    """Get metadata for a class by type.

    Args:
        cls: The class type

    Returns:
        ClassMetadata instance or None if not found
    """
    return getattr(cls, "_metadata", None)


class MetadataAttributeFactory:
    """Factory for creating attribute objects from metadata definitions."""

    @classmethod
    def create_attribute(
        cls, name: str, attr_def: AttributeDefinition, parent_obj
    ) -> Any:
        """Create an attribute object from definition, sourcing qualifiers from class-level metadata."""
        from .base_classes import ValueState

        # Get class-level qualifiers from parent_obj's class metadata
        qualifiers = {}
        meta = getattr(parent_obj.__class__, '_metadata', None)
        if meta and hasattr(meta, 'qualifiers') and meta.qualifiers:
            qualifiers = meta.qualifiers

        # Helper to get qualifier value, fallback to attribute definition
        def q(key, default=None):
            return qualifiers.get(key, getattr(attr_def, key, default))

        # Patch: pass qualifiers to attribute creation
        if attr_def.attr_type == AttributeType.INT:
            return cls._create_int_attribute(name, attr_def, parent_obj, qualifiers)
        elif attr_def.attr_type == AttributeType.FLOAT:
            return cls._create_float_attribute(name, attr_def, parent_obj, qualifiers)
        elif attr_def.attr_type == AttributeType.BOOLEAN:
            return cls._create_boolean_attribute(name, attr_def, parent_obj, qualifiers)
        elif attr_def.attr_type == AttributeType.STRING:
            return cls._create_string_attribute(name, attr_def, parent_obj, qualifiers)
        elif attr_def.attr_type == AttributeType.CUSTOM:
            return cls._create_custom_attribute(name, attr_def, parent_obj, qualifiers)
        else:
            raise ValueError(f"Unknown attribute type: {attr_def.attr_type}")

    @classmethod
    def _create_int_attribute(
        cls, name: str, attr_def: AttributeDefinition, parent_obj, qualifiers
    ):
        """Create an integer attribute using CInt."""
        from .fundamental_types import CInt, ValueState

        # Create CInt object without a value (so it stays NOT_SET)
        attr = CInt(parent=parent_obj, name=name)

        # Set default value from qualifiers if explicitly provided
        # IMPORTANT: Set as DEFAULT state, not EXPLICITLY_SET
        default_value = qualifiers.get('default')
        if default_value is not None:
            attr._value = int(default_value)
            attr._value_states['value'] = ValueState.DEFAULT

        # CInt already has all necessary methods (_is_value_type, __int__, __str__, etc.)
        # The min/max validation will be handled by qualifiers at the class level

        return attr

    @classmethod
    def _create_float_attribute(
        cls, name: str, attr_def: AttributeDefinition, parent_obj, qualifiers
    ):
        """Create a float attribute using CFloat."""
        from .fundamental_types import CFloat, ValueState

        # Create CFloat object without a value (so it stays NOT_SET)
        attr = CFloat(parent=parent_obj, name=name)

        # Set default value from qualifiers if explicitly provided
        # IMPORTANT: Set as DEFAULT state, not EXPLICITLY_SET
        default_value = qualifiers.get('default')
        if default_value is not None:
            attr._value = float(default_value)
            attr._value_states['value'] = ValueState.DEFAULT

        # CFloat already has all necessary methods (_is_value_type, __float__, __str__, etc.)
        # The min/max validation will be handled by qualifiers at the class level

        return attr

    @classmethod
    def _create_boolean_attribute(
        cls, name: str, attr_def: AttributeDefinition, parent_obj, qualifiers
    ):
        """Create a boolean attribute using CBoolean."""
        from .fundamental_types import CBoolean, ValueState

        # Create CBoolean object without a value (so it stays NOT_SET)
        attr = CBoolean(parent=parent_obj, name=name)

        # Set default value from qualifiers if explicitly provided
        # IMPORTANT: Set as DEFAULT state, not EXPLICITLY_SET
        default_value = qualifiers.get('default')
        if default_value is not None:
            attr._value = bool(default_value)
            attr._value_states['value'] = ValueState.DEFAULT

        # CBoolean already has all necessary methods (_is_value_type, __bool__, __str__, etc.)

        return attr

    @classmethod
    def _create_string_attribute(
        cls, name: str, attr_def: AttributeDefinition, parent_obj, qualifiers
    ):
        """Create a string-type attribute, sourcing default from qualifiers."""
        from .fundamental_types import CString

        attr = CString(parent=parent_obj, name=name)

        # Set default value from qualifiers if provided
        default_value = qualifiers.get('default')
        if default_value is not None:
            attr.value = str(default_value)

        return attr

    @classmethod
    def _create_custom_attribute(
        cls, name: str, attr_def: AttributeDefinition, parent_obj, qualifiers
    ):
        """Create a custom attribute type using the custom_class specification."""
        from .base_classes import CData, ValueState

        # Get the custom class name
        custom_class_name = attr_def.custom_class
        if not custom_class_name:
            # Fallback to CString if no custom class specified
            return cls._create_string_attribute(name, attr_def, parent_obj, qualifiers)

        # Build a class registry (similar to DEF XML parser)
        custom_class = cls._get_class_from_registry(custom_class_name)

        if custom_class is None:
            # Class not found, fallback to CString
            print(f"Warning: Custom class '{custom_class_name}' not found for attribute '{name}', using CString")
            return cls._create_string_attribute(name, attr_def, parent_obj, qualifiers)

        # Create instance of the custom class
        try:
            obj = custom_class(parent=parent_obj, name=name)

            # Apply qualifiers
            if qualifiers:
                # Ensure _qualifiers attribute exists
                if not hasattr(obj, '_qualifiers') or obj._qualifiers is None:
                    obj._qualifiers = {}
                # Update qualifiers
                if isinstance(obj._qualifiers, dict):
                    obj._qualifiers.update(qualifiers)
                else:
                    obj._qualifiers = qualifiers

            # Set default value from qualifiers if provided
            default_value = qualifiers.get('default')
            if default_value is not None and hasattr(obj, 'value'):
                obj.value = default_value
                if hasattr(parent_obj, "_value_states"):
                    parent_obj._value_states[name] = ValueState.DEFAULT

            return obj
        except Exception as e:
            print(f"Warning: Failed to create custom class '{custom_class_name}': {e}")
            return cls._create_string_attribute(name, attr_def, parent_obj, qualifiers)

    @classmethod
    def _get_class_from_registry(cls, class_name: str):
        """Get a class from the registry, building it if needed.

        If class_name ends with 'Stub', tries to find the implementation class first
        (without the Stub suffix), falling back to the stub if not found.
        """
        # Import here to avoid circular dependencies
        from .fundamental_types import CInt, CFloat, CBoolean, CString, CList
        from .base_classes import CContainer

        # Build a simple registry - only fundamental types and container
        # Other classes are imported dynamically below
        registry = {
            "CInt": CInt,
            "CFloat": CFloat,
            "CBoolean": CBoolean,
            "CString": CString,
            "CContainer": CContainer,
            "CList": CList,
        }

        # Try to get from basic registry first
        if class_name in registry:
            return registry[class_name]

        # If class_name ends with "Stub", try to get implementation class first
        impl_class_name = class_name
        if class_name.endswith('Stub'):
            impl_class_name = class_name[:-4]  # Remove "Stub" suffix

        # Try to import from core implementation and stub modules
        try:
            import importlib
            # Try importing from various possible locations
            # Order matters - try implementation first, then stubs
            possible_modules = [
                # Try implementation classes first (without Stub suffix)
                ('ccp4i2.core.CCP4Data', impl_class_name),
                ('ccp4i2.core.CCP4ModelData', impl_class_name),
                ('ccp4i2.core.CCP4File', impl_class_name),
                ('ccp4i2.core.CCP4XtalData', impl_class_name),
                ('ccp4i2.core.CCP4Annotation', impl_class_name),
                ('ccp4i2.core.CCP4RefmacData', impl_class_name),
                ('ccp4i2.core.CCP4MathsData', impl_class_name),
                (f'ccp4i2.core.{impl_class_name}', impl_class_name),

                # Then try stub classes
                ('ccp4i2.core.cdata_stubs.CCP4Data', class_name),
                ('ccp4i2.core.cdata_stubs.CCP4ModelData', class_name),
                ('ccp4i2.core.cdata_stubs.CCP4File', class_name),
                ('ccp4i2.core.cdata_stubs.CCP4XtalData', class_name),
                (f'ccp4i2.core.cdata_stubs.{class_name}', class_name),
            ]

            for module_path, lookup_name in possible_modules:
                try:
                    module = importlib.import_module(module_path)
                    if hasattr(module, lookup_name):
                        return getattr(module, lookup_name)
                except (ImportError, AttributeError):
                    continue
        except Exception:
            pass

        return None


def apply_metadata_to_instance(instance):
    """Apply metadata-defined attributes to a class instance.

    Args:
        instance: The instance to apply metadata to
    """
    # Collect attributes and content_qualifiers from all ancestor classes with metadata
    # Walk MRO in REVERSE order so that child classes override parent classes
    merged_attributes = {}
    merged_content_qualifiers = {}
    for cls in reversed(instance.__class__.__mro__):
        if cls is object:
            continue
        metadata = getattr(cls, "_metadata", None)
        if metadata:
            # Parent attributes are added first (because we're in reverse),
            # then child overrides them
            merged_attributes.update(metadata.attributes)
            # Same for content_qualifiers (per-field qualifiers from CONTENTS)
            if metadata.content_qualifiers:
                merged_content_qualifiers.update(metadata.content_qualifiers)

    # Create attributes from merged metadata
    for attr_name, attr_def in merged_attributes.items():
        # Check if attribute exists in hierarchy (via _children_by_name cache)
        # NOT in __dict__ - CData children are stored in hierarchy, not __dict__
        # This is important because generated classes have type annotations like
        # `label: Optional[COneWord] = None` which creates a class attribute,
        # but we want to replace it with an actual COneWord instance
        children_by_name = getattr(instance, '_children_by_name', {})
        if attr_name not in children_by_name:
            attr_obj = MetadataAttributeFactory.create_attribute(
                attr_name, attr_def, instance
            )
            # Add to hierarchy via set_parent (populates _children_by_name for O(1) access)
            # DON'T store in __dict__ - hierarchy is the single source of truth
            if hasattr(attr_obj, 'set_parent'):
                attr_obj._name = attr_name  # Set hierarchical name
                attr_obj.set_parent(instance)  # This adds to _children_by_name
            else:
                # Fallback for non-CData objects (shouldn't happen for metadata attributes)
                instance.__dict__[attr_name] = attr_obj

            # Apply per-field qualifiers from content_qualifiers
            if attr_name in merged_content_qualifiers and attr_obj is not None:
                field_qualifiers = merged_content_qualifiers[attr_name]
                if hasattr(attr_obj, 'set_qualifier'):
                    for qual_name, qual_value in field_qualifiers.items():
                        attr_obj.set_qualifier(qual_name, qual_value)

                # If there's a default value, apply it to the attribute's value
                # This ensures the UI shows the correct default (e.g., polymerType = "PROTEIN")
                # and the value is properly initialized for validation
                if 'default' in field_qualifiers and hasattr(attr_obj, 'value'):
                    default_value = field_qualifiers['default']
                    if default_value is not None:
                        # Set value but mark as DEFAULT state, not EXPLICITLY_SET
                        # This allows isSet(allowDefault=False) to correctly identify
                        # fields that haven't been user-modified
                        attr_obj._value = default_value if hasattr(attr_obj, '_value') else None
                        if hasattr(attr_obj, '_value_states'):
                            from .base_classes import ValueState
                            attr_obj._value_states['value'] = ValueState.DEFAULT

            if hasattr(instance, "_value_states"):
                from .base_classes import ValueState

                instance._value_states[attr_name] = ValueState.NOT_SET
