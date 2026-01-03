"""
Modern metadata system for CData classes.
Provides a Pythonic way to capture and use the rich metadata from CCP4i2 qualifiers.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union, Type, get_type_hints
from enum import Enum
import inspect


@dataclass
class FieldMetadata:
    """Metadata for a single field/attribute."""

    # Basic properties
    name: str
    data_type: Type = Any
    default: Any = None

    # Validation properties
    minlength: Optional[int] = None
    maxlength: Optional[int] = None
    minimum: Optional[Union[int, float]] = None
    maximum: Optional[Union[int, float]] = None

    # Enumeration properties
    enumerators: Optional[List[str]] = None
    menu_text: Optional[List[str]] = None
    only_enumerators: bool = False

    # UI properties
    tooltip: Optional[str] = None
    gui_label: Optional[str] = None
    widget_type: Optional[str] = None

    # Business logic properties
    required: bool = False
    readonly: bool = False
    hidden: bool = False

    def validate(self, value: Any) -> List[str]:
        """Validate a value against this field's metadata."""
        errors = []

        if self.required and (value is None or value == ""):
            errors.append(f"{self.name} is required")

        if value is not None:
            # Length validation for strings
            if isinstance(value, str):
                if self.minlength and len(value) < self.minlength:
                    errors.append(
                        f"{self.name} must be at least {self.minlength} characters"
                    )
                if self.maxlength and len(value) > self.maxlength:
                    errors.append(
                        f"{self.name} must be at most {self.maxlength} characters"
                    )

            # Numeric validation
            if isinstance(value, (int, float)):
                if self.minimum is not None and value < self.minimum:
                    errors.append(f"{self.name} must be at least {self.minimum}")
                if self.maximum is not None and value > self.maximum:
                    errors.append(f"{self.name} must be at most {self.maximum}")

            # Enumeration validation
            if self.only_enumerators and self.enumerators:
                if str(value) not in self.enumerators:
                    errors.append(
                        f"{self.name} must be one of: {', '.join(self.enumerators)}"
                    )

        return errors


@dataclass
class ClassMetadata:
    """Metadata for an entire class."""

    name: str
    docstring: Optional[str] = None
    fields: Dict[str, FieldMetadata] = field(default_factory=dict)
    error_codes: Dict[int, str] = field(default_factory=dict)

    # Inheritance information
    base_classes: List[str] = field(default_factory=list)
    immediate_parent: Optional[str] = None

    # File/module information
    source_file: Optional[str] = None
    source_module: Optional[str] = None

    def get_field_metadata(self, field_name: str) -> Optional[FieldMetadata]:
        """Get metadata for a specific field."""
        return self.fields.get(field_name)

    def get_required_fields(self) -> List[str]:
        """Get list of required field names."""
        return [name for name, meta in self.fields.items() if meta.required]

    def get_enumerated_fields(self) -> Dict[str, List[str]]:
        """Get fields that have enumeration constraints."""
        return {
            name: meta.enumerators
            for name, meta in self.fields.items()
            if meta.enumerators
        }

    def validate_instance(self, instance: Any) -> List[str]:
        """Validate an entire class instance."""
        errors = []

        for field_name, field_meta in self.fields.items():
            value = getattr(instance, field_name, None)
            field_errors = field_meta.validate(value)
            errors.extend(field_errors)

        return errors


class MetadataRegistry:
    """Global registry for class metadata."""

    _metadata: Dict[str, ClassMetadata] = {}

    @classmethod
    def register(cls, class_name: str, metadata: ClassMetadata):
        """Register metadata for a class."""
        cls._metadata[class_name] = metadata

    @classmethod
    def get(cls, class_name: str) -> Optional[ClassMetadata]:
        """Get metadata for a class."""
        return cls._metadata.get(class_name)

    @classmethod
    def get_by_class(cls, class_obj: Type) -> Optional[ClassMetadata]:
        """Get metadata by class object."""
        return cls._metadata.get(class_obj.__name__)

    @classmethod
    def list_classes(cls) -> List[str]:
        """List all registered classes."""
        return list(cls._metadata.keys())


def metadata_field(
    default: Any = None,
    tooltip: str = None,
    required: bool = False,
    minlength: int = None,
    maxlength: int = None,
    minimum: Union[int, float] = None,
    maximum: Union[int, float] = None,
    enumerators: List[str] = None,
    menu_text: List[str] = None,
    only_enumerators: bool = False,
    readonly: bool = False,
    hidden: bool = False,
    widget_type: str = None,
    gui_label: str = None,
):
    """
    Decorator-like function to attach metadata to fields.
    Usage: field_name: str = metadata_field(default="test", tooltip="Field description")
    """
    return {
        "default": default,
        "tooltip": tooltip,
        "required": required,
        "minlength": minlength,
        "maxlength": maxlength,
        "minimum": minimum,
        "maximum": maximum,
        "enumerators": enumerators,
        "menu_text": menu_text,
        "only_enumerators": only_enumerators,
        "readonly": readonly,
        "hidden": hidden,
        "widget_type": widget_type,
        "gui_label": gui_label,
    }


def with_metadata(cls):
    """
    Class decorator to automatically extract and register metadata.
    """
    # Extract field metadata from type hints and defaults
    fields = {}
    type_hints = get_type_hints(cls)

    for attr_name, attr_type in type_hints.items():
        if not attr_name.startswith("_"):
            # Get default value and metadata
            default_value = getattr(cls, attr_name, None)

            field_meta = FieldMetadata(
                name=attr_name, data_type=attr_type, default=default_value
            )

            # If default_value is a metadata dict, extract it
            if isinstance(default_value, dict) and "tooltip" in default_value:
                for key, value in default_value.items():
                    if hasattr(field_meta, key):
                        setattr(field_meta, key, value)
                # Set the actual default
                field_meta.default = default_value.get("default")

            fields[attr_name] = field_meta

    # Create and register class metadata
    class_meta = ClassMetadata(name=cls.__name__, docstring=cls.__doc__, fields=fields)

    MetadataRegistry.register(cls.__name__, class_meta)

    # Add convenience methods to the class
    def get_metadata(self):
        return MetadataRegistry.get_by_class(self.__class__)

    def validate(self):
        metadata = self.get_metadata()
        if metadata:
            return metadata.validate_instance(self)
        return []

    def get_field_info(self, field_name: str):
        metadata = self.get_metadata()
        if metadata:
            return metadata.get_field_metadata(field_name)
        return None

    # Add methods to class
    cls.get_metadata = get_metadata
    cls.validate = validate
    cls.get_field_info = get_field_info

    return cls


# Error code system
class ErrorCodes(Enum):
    """Standard error codes for validation."""

    REQUIRED_FIELD_MISSING = 100
    VALUE_TOO_SHORT = 101
    VALUE_TOO_LONG = 102
    VALUE_TOO_SMALL = 103
    VALUE_TOO_LARGE = 104
    INVALID_ENUMERATION = 105
    INVALID_TYPE = 106
    READONLY_FIELD = 107


def get_error_message(error_code: Union[int, ErrorCodes]) -> str:
    """Get human-readable error message for error code."""
    error_messages = {
        ErrorCodes.REQUIRED_FIELD_MISSING: "Required field is missing or empty",
        ErrorCodes.VALUE_TOO_SHORT: "Value is shorter than minimum length",
        ErrorCodes.VALUE_TOO_LONG: "Value is longer than maximum length",
        ErrorCodes.VALUE_TOO_SMALL: "Value is smaller than minimum",
        ErrorCodes.VALUE_TOO_LARGE: "Value is larger than maximum",
        ErrorCodes.INVALID_ENUMERATION: "Value is not in allowed list",
        ErrorCodes.INVALID_TYPE: "Value is not of expected type",
        ErrorCodes.READONLY_FIELD: "Field is read-only and cannot be modified",
    }

    if isinstance(error_code, int):
        # Try to find matching enum
        for enum_val in ErrorCodes:
            if enum_val.value == error_code:
                return error_messages.get(enum_val, f"Unknown error code: {error_code}")

    return error_messages.get(error_code, f"Unknown error: {error_code}")
