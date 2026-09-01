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
    """Complete metadata for a CData class.

    Note: Qualifier values are stored on the class as _qualifiers_template,
    NOT in this metadata object. This metadata holds structural information
    (attributes, ordering, error codes) that doesn't change per-instance.
    The ``qualifiers`` property provides read access to the owning class's
    ``_qualifiers_template`` for backwards compatibility.
    """

    attributes: Dict[str, AttributeDefinition] = field(default_factory=dict)
    error_codes: Dict[int, str] = field(default_factory=dict)
    docstring: Optional[str] = None
    file_extensions: Optional[List[str]] = None
    mime_type: Optional[str] = None
    gui_label: Optional[str] = None
    contents_order: Optional[List[str]] = None
    qualifiers_order: Optional[List[str]] = None
    content_qualifiers: Optional[Dict[str, Dict[str, Any]]] = None  # Per-field qualifiers
    _owner_class: Optional[Type] = field(default=None, repr=False)

    @property
    def qualifiers(self) -> Optional[Dict[str, Any]]:
        """Read-through to the owning class's _qualifiers_template."""
        if self._owner_class is not None:
            return getattr(self._owner_class, '_qualifiers_template', None)
        return None


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
    content_qualifiers: Optional[Dict[str, Dict[str, Any]]] = None,
):
    """Class decorator to add metadata to CData classes.

    Qualifier system:
        Qualifiers (constraints, defaults, display hints) are stored as a single
        class-level template dict (cls._qualifiers_template). At instance creation,
        CData.__init__ copies this to self._qualifiers, which is the sole runtime
        source of truth. No other copies are made.

    Args:
        attributes: Dictionary of attribute name -> AttributeDefinition
        qualifiers: Dictionary of class qualifiers (stored as _qualifiers_template)
        error_codes: Dictionary of error code -> message
        file_extensions: List of supported file extensions
        mime_type: MIME type for file classes
        gui_label: Label for GUI display
        contents_order: List specifying display order of attributes in UI
        qualifiers_order: List specifying display order of qualifiers
        content_qualifiers: Per-field qualifiers for child attributes (from CONTENTS)

    Example:
        @cdata_class(
            attributes={
                'project': attribute(AttributeType.CUSTOM, custom_class="CProjectId"),
                'baseName': attribute(AttributeType.CUSTOM, custom_class="CFilePath"),
            },
            qualifiers={'min': None, 'max': None, 'default': 50},
            file_extensions=['dat', 'txt'],
            mime_type='text/plain'
        )
        class CDataFile(CData):
            '''A data file with embedded metadata.'''
            pass
    """

    def decorator(cls: Type) -> Type:
        # Create metadata object (structural info only, no qualifier values)
        metadata = ClassMetadata(
            attributes=attributes or {},
            error_codes=error_codes or {},
            docstring=cls.__doc__,
            file_extensions=file_extensions,
            mime_type=mime_type,
            gui_label=gui_label,
            contents_order=contents_order,
            qualifiers_order=qualifiers_order,
            content_qualifiers=content_qualifiers,
        )

        # Store in global registry
        _CLASS_METADATA_REGISTRY[cls.__name__] = metadata

        # Store structural metadata as class attribute (with back-reference)
        metadata._owner_class = cls
        cls._metadata = metadata

        # Single canonical source for qualifier default values.
        # CData.__init__ copies this to self._qualifiers per-instance.
        if qualifiers:
            cls._qualifiers_template = qualifiers

        # Other class-level attributes
        if contents_order:
            cls.CONTENT_ORDER = contents_order
        # ERROR_CODES is a plain class attribute, so a subclass declaring its
        # own used to *shadow* its ancestors' entirely --- declare one code and
        # you lost all the inherited ones. That is why every generated stub
        # carries a complete copy: 2,889 entries of which only 346 are
        # distinct, the same code repeated up to 44 times.
        #
        # Merging along the MRO makes inheritance work, so a class can declare
        # only what it adds. It is purely additive: measured over the stubs it
        # changes nothing for 173 classes, adds ancestor codes to 66 that were
        # missing them --- `self.ERROR_CODES[105]` would have raised there ---
        # and redefines nothing anywhere, so no lookup can change its answer.
        merged = {}
        for base in reversed(cls.__mro__[1:]):
            merged.update(base.__dict__.get("ERROR_CODES", {}) or {})
        if error_codes:
            merged.update(error_codes)
        # A class may also declare ERROR_CODES in its own body, by hand. The
        # decorator runs after the body, so assigning here would discard it ---
        # and silently, since the codes look present until one is looked up.
        # CRangeSelection does exactly this: the decorator carries "201" as a
        # string and the body carries 201 as an int, and validity() reads the
        # int. Its own body wins.
        merged.update(cls.__dict__.get("ERROR_CODES", {}) or {})
        if merged:
            cls.ERROR_CODES = merged

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
        """Create an attribute object from definition, sourcing qualifiers from class-level template."""
        from .base_classes import ValueState

        # Get class-level qualifiers from the single canonical source
        qualifiers = getattr(parent_obj.__class__, '_qualifiers_template', {})

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
                    from .cdata import _coerced_qualifier
                    obj._qualifiers.update(
                        {k: _coerced_qualifier(k, v) for k, v in qualifiers.items()})
                else:
                    from .cdata import _coerced_qualifier
                    obj._qualifiers = {k: _coerced_qualifier(k, v)
                                       for k, v in qualifiers.items()}

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

        Class names arrive as strings --- from `content()`, from a quoted
        annotation, from a def.xml `<className>` --- and are resolved here, which
        is what lets a field name a class defined later or in a module that
        would import circularly.

        The 'Stub' suffix is still stripped on the way in. No class carries it
        any more and nothing in the tree writes one, but a saved job or an
        out-of-tree wrapper may name one, and mapping it to the merged class is
        exactly right.
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
                ('ccp4i2.core.CCP4Data', class_name),
                ('ccp4i2.core.CCP4ModelData', class_name),
                ('ccp4i2.core.CCP4File', class_name),
                ('ccp4i2.core.CCP4XtalData', class_name),
                (f'ccp4i2.core.{class_name}', class_name),
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


#: The annotation on a stub is already a complete declaration of its fields:
#: `structure: Optional[CPdbDataFile] = None`. Measured across the 122
#: stub classes, the annotation set and the decorator's `attributes=` agree
#: exactly --- 621 fields, and for the 278 naming a custom class the annotation
#: names the same one, while the other 343 map one-to-one
#: (STRING->CString, INT->CInt, FLOAT->CFloat, BOOLEAN->CBoolean).
#:
#: So the decorator's attributes= is a second copy of something Python already
#: records, and a dataclass would use the annotation. This reads it, so the
#: duplicate can go and the annotation becomes the single declaration.
_ANNOTATION_KINDS = {
    "CString": AttributeType.STRING,
    "CInt": AttributeType.INT,
    "CFloat": AttributeType.FLOAT,
    "CBoolean": AttributeType.BOOLEAN,
}


def _class_named_in(annotation):
    """The CData class an annotation names, e.g. Optional[CPdbDataFile].

    `typing.Optional[X]` is not a class and its __name__ is 'Optional', so it
    has to be unwrapped rather than stringified. Forward references arrive as
    strings and keep their text.
    """
    import typing
    if annotation is None:
        return None
    for arg in typing.get_args(annotation) or ():
        if arg is type(None):
            continue
        named = _class_named_in(arg)
        if named:
            return named
    if isinstance(annotation, str):
        return annotation.split("[")[-1].rstrip("]").split(".")[-1].strip("'\" ") or None
    if isinstance(annotation, typing.ForwardRef):
        return annotation.__forward_arg__.split(".")[-1]
    name = getattr(annotation, "__name__", None)
    return name if name and name not in ("Optional", "Union", "NoneType") else None


class Content:
    """One declared field of a CData class: its type and its own qualifiers.

    Deliberately *not* a `dataclasses.field`. The declaration reads similarly
    and that is the point --- but a dataclass makes promises CData cannot keep,
    and importing the syntax would import the expectations with it:

      - `CCell() == CCell()` would compare field by field; CData compares by
        identity, and `isSet()` compares a value against its default.
      - `hash()` would be disabled, and CData objects are hashable.
      - `asdict()` would recurse into a tree of live parented objects.
      - `a: CCellLength` would say "assign a CCellLength here", when in fact
        `cell.a = 55.0` coerces and the field stays a CCellLength. That is the
        plugin API, at 231 call sites.

    So `is_dataclass()` on a CData class stays False, and a reader who knows
    dataclasses is not invited to assume the rest of the protocol.

    What this *does* carry is the pair that belongs together and used to be
    written apart --- the field's class, and the qualifiers for that field,
    which lived in a class-level `content_qualifiers` dict keyed by name.

    Qualifiers given here are the *declaration*. They seed the instance's own
    `_qualifiers`, which stays mutable and private: 56 sites change a qualifier
    at runtime, and `test_qualifiers_are_per_instance` pins the isolation.
    """

    __slots__ = ("cls", "qualifiers", "order", "name")
    _counter = 0

    def __init__(self, cls, **qualifiers):
        self.cls = cls
        self.qualifiers = qualifiers
        Content._counter += 1
        self.order = Content._counter          # declaration order, for py<3.7 safety

    def __set_name__(self, owner, name):
        """Record the declaration on the class as the class body runs.

        Not for elegance: a plain class attribute can be overwritten by a later
        statement in the same body, and several classes do exactly that ---
        CResolutionRange declares `low` and `high` and then defines properties
        of the same name, which replace them in the class dict. Annotations
        survive that because they live in `__annotations__`, a separate
        mapping; assignments do not.

        `__set_name__` fires at the moment of assignment, so the declaration is
        captured before anything can shadow it.
        """
        self.name = name
        # Must be the class's *own* mapping, not one inherited from a base ---
        # otherwise a subclass's fields would be recorded onto its parent.
        if "__cdata_contents__" not in vars(owner):
            setattr(owner, "__cdata_contents__", {})
        vars(owner)["__cdata_contents__"][name] = self

    def __repr__(self):
        name = self.cls if isinstance(self.cls, str) else getattr(self.cls, "__name__", self.cls)
        return f"content({name}{', ' if self.qualifiers else ''}{', '.join(f'{k}={v!r}' for k, v in self.qualifiers.items())})"


def content(cls, **qualifiers):
    """Declare a field of a CData class.

        class CCell(CData):
            a = content(CCellLength, toolTip='Cell length a in A', guiLabel='a')

    `cls` may be the class or its name; a name is resolved when the field is
    built, which is how a field can refer to a class defined later or in a
    module that would import circularly.
    """
    return Content(cls, **qualifiers)


def contents_from_declarations(cls):
    """The `content()` declarations made by `cls` itself, in declaration order."""
    declared = vars(cls).get("__cdata_contents__") or {}
    return dict(sorted(declared.items(), key=lambda kv: kv[1].order))


def attributes_from_annotations(cls) -> "Dict[str, AttributeDefinition]":
    """Field declarations read from the class's own annotations."""
    declared = {}
    for name, annotation in getattr(cls, "__annotations__", {}).items():
        if name.startswith("_"):
            continue
        named = _class_named_in(annotation)
        if not named:
            continue
        kind = _ANNOTATION_KINDS.get(named)
        if kind is not None:
            declared[name] = AttributeDefinition(attr_type=kind)
        else:
            declared[name] = AttributeDefinition(
                attr_type=AttributeType.CUSTOM,
                custom_class=named,
            )
    return declared


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
            # then child overrides them.
            #
            # The class's own annotations are preferred where it has them,
            # because that is the declaration a dataclass would use and the
            # decorator's attributes= is a second copy of it. Measured across
            # the 122 generated stubs the two agree exactly; the fallback is
            # for hand-written classes such as CDataFile, which carry the
            # decorator but no annotations.
            # A class may declare its fields with content(), which carries the
            # field's class and its own qualifiers together. Where it does,
            # that is the declaration; annotations and the decorator's
            # attributes= are the two older spellings of the same thing.
            # Both spellings are read, so a class may use either or both.
            # content() wins for the names it declares, because it is the one
            # that carries the field's qualifiers alongside its class.
            from_annotations = attributes_from_annotations(cls)
            merged_attributes.update(from_annotations or metadata.attributes)
            for _name, _decl in contents_from_declarations(cls).items():
                _cls = _decl.cls if isinstance(_decl.cls, str) else _decl.cls.__name__
                # The fundamental types get their own AttributeType, exactly as
                # attributes_from_annotations gives them. They are not merely a
                # different spelling of CUSTOM: the two construction paths treat
                # qualifiers differently, and routing a CString through CUSTOM
                # puts the *class's* qualifiers onto the child --- CAnnotation's
                # label and toolTip landed on its `text` field that way.
                _kind = _ANNOTATION_KINDS.get(_cls)
                merged_attributes[_name] = (
                    attribute(_kind) if _kind is not None
                    else attribute(AttributeType.CUSTOM, custom_class=_cls))
                if _decl.qualifiers:
                    # Replace, as a class-level content_qualifiers entry does:
                    # the declaration here is this class's whole statement about
                    # the field, not an overlay on an ancestor's.
                    merged_content_qualifiers[_name] = dict(_decl.qualifiers)
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
        # Children live in __dict__ now, keyed by name.
        from .hierarchy_system import HierarchicalObject
        existing_child = instance.__dict__.get(attr_name)
        if not isinstance(existing_child, HierarchicalObject):
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
