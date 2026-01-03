import logging
from ccp4i2.core import CCP4Data
from ccp4i2.core.base_object.cdata import CData
from ccp4i2.core.base_object.fundamental_types import CList, CFloat, CInt, CString, CBoolean
from ccp4i2.core.base_object.base_classes import CContainer

logger = logging.getLogger(f"ccp4i2:{__name__}")


def value_dict_for_object(ccp4i2_object):
    """
    Convert a CData object tree to nested Python primitives (dicts, lists, str, int, float, bool).

    Rules:
    - Fundamental types (CInt, CFloat, CString, CBoolean) → Python primitive (int, float, str, bool)
    - CList → Python list of converted items
    - CContainer/CData with children → dict of {child_name: converted_child}
    - Regular CData objects → dict of {attribute_name: converted_attribute}
    - Plain Python types (dict, list, str, int, float, bool, None) → pass through
    """
    logger.debug("In value_dict_for_object: %s", type(ccp4i2_object).__name__)

    # Handle None
    if ccp4i2_object is None:
        return None

    # Handle plain Python primitives - pass through
    if isinstance(ccp4i2_object, (str, int, float, bool)):
        return ccp4i2_object

    # Handle plain Python dict
    if isinstance(ccp4i2_object, dict):
        return handle_dict(ccp4i2_object)

    # Handle plain Python list
    if isinstance(ccp4i2_object, list):
        return handle_list(ccp4i2_object)

    # Handle fundamental CData types - extract Python value
    if isinstance(ccp4i2_object, (CInt, CFloat, CString, CBoolean)):
        return handle_fundamental_type(ccp4i2_object)

    # Handle CList - convert to Python list
    if isinstance(ccp4i2_object, (CList, CCP4Data.CList)):
        return handle_clist(ccp4i2_object)

    # Handle CContainer - convert children to dict
    if isinstance(ccp4i2_object, CContainer):
        return handle_container(ccp4i2_object)

    # Handle general CData objects - convert attributes to dict
    if isinstance(ccp4i2_object, CData):
        # Check if object has a custom to_dict() method (e.g., CDataFileContent subclasses)
        if hasattr(ccp4i2_object, 'to_dict') and callable(ccp4i2_object.to_dict):
            custom_dict = ccp4i2_object.to_dict()
            if custom_dict is not None:
                return custom_dict
        # Fall back to generic CData handling
        return handle_cdata(ccp4i2_object)

    # Handle other objects with __dict__ (e.g., CPdbDataComposition)
    if hasattr(ccp4i2_object, '__dict__'):
        return handle_plain_object(ccp4i2_object)

    # Fallback for unknown types
    logger.warning(
        "Unknown type %s, returning string representation",
        type(ccp4i2_object)
    )
    return str(ccp4i2_object)


def handle_fundamental_type(ccp4i2_object):
    """Extract Python primitive from fundamental CData types (CInt, CFloat, CString, CBoolean)."""
    logger.debug("Handling fundamental type: %s", type(ccp4i2_object).__name__)

    # Use the .value property if available (modern CData)
    if hasattr(ccp4i2_object, 'value'):
        return ccp4i2_object.value

    # Fallback to _value for older code
    if hasattr(ccp4i2_object, '_value'):
        return ccp4i2_object._value

    # Last resort: convert to string
    return str(ccp4i2_object)


def handle_clist(ccp4i2_object):
    """Convert CList to Python list, recursively converting items."""
    logger.debug("Handling CList: %s", type(ccp4i2_object).__name__)
    result = []

    # CList should be iterable
    try:
        for item in ccp4i2_object:
            result.append(value_dict_for_object(item))
    except Exception as e:
        logger.warning("Failed to iterate CList: %s", e)
        return []

    return result


def handle_container(ccp4i2_object):
    """Convert CContainer to dict of {child_name: child_value_dict}."""
    logger.debug("Handling CContainer: %s", type(ccp4i2_object).__name__)
    result = {}

    # Get all children from the container
    # Try both children() (HierarchicalObject) and get_children() (legacy) methods
    try:
        if hasattr(ccp4i2_object, 'children') and callable(ccp4i2_object.children):
            children = ccp4i2_object.children()
        elif hasattr(ccp4i2_object, 'get_children'):
            children = ccp4i2_object.get_children()
        else:
            children = []

        for child in children:
            # Get child name - try objectName() first, then object_name(), then name attribute
            if hasattr(child, 'objectName') and callable(child.objectName):
                name = child.objectName()
            elif hasattr(child, 'object_name') and callable(child.object_name):
                name = child.object_name()
            elif hasattr(child, 'name'):
                name = child.name
            else:
                name = str(child)
            result[name] = value_dict_for_object(child)
    except Exception as e:
        logger.warning("Failed to get children from container: %s", e)

    return result if result else {}


def handle_cdata(ccp4i2_object):
    """Convert regular CData object to dict of its attributes."""
    logger.debug("Handling CData: %s", type(ccp4i2_object).__name__)
    result = {}

    # Get registered attributes from metadata
    if hasattr(ccp4i2_object, 'get_merged_metadata'):
        try:
            metadata = ccp4i2_object.get_merged_metadata('attributes')
            if metadata:
                for attr_name in metadata.keys():
                    if hasattr(ccp4i2_object, attr_name):
                        attr_value = getattr(ccp4i2_object, attr_name)
                        # Only include if the attribute is set
                        if hasattr(ccp4i2_object, 'isSet') and ccp4i2_object.isSet(attr_name):
                            result[attr_name] = value_dict_for_object(attr_value)
        except Exception as e:
            logger.debug("Failed to get metadata attributes: %s", e)

    # Second strategy: try to get children (in case it has a hierarchical structure)
    # Try both children() (HierarchicalObject) and get_children() (legacy) methods
    if not result:
        try:
            children = None
            if hasattr(ccp4i2_object, 'children') and callable(ccp4i2_object.children):
                children = ccp4i2_object.children()
            elif hasattr(ccp4i2_object, 'get_children'):
                children = ccp4i2_object.get_children()

            if children:
                for child in children:
                    # Get child name - try objectName() first, then object_name(), then name attribute
                    if hasattr(child, 'objectName') and callable(child.objectName):
                        name = child.objectName()
                    elif hasattr(child, 'object_name') and callable(child.object_name):
                        name = child.object_name()
                    elif hasattr(child, 'name'):
                        name = child.name
                    else:
                        name = str(child)
                    result[name] = value_dict_for_object(child)
        except Exception as e:
            logger.debug("Failed to get children: %s", e)

    # Third strategy: inspect __dict__ for public attributes
    # This handles objects where data is stored directly
    if not result:
        try:
            skip_attrs = {
                'destroyed', 'parent_changed', 'child_added', 'child_removed',
                'qualifiers_order', 'qualifiers_definition'
            }
            for attr_name, attr_value in ccp4i2_object.__dict__.items():
                # Skip private attributes, signals, and metadata
                if attr_name.startswith('_'):
                    continue
                if attr_name in skip_attrs:
                    continue
                # Include CData, basic types, lists, dicts, or custom objects
                result[attr_name] = value_dict_for_object(attr_value)
        except Exception as e:
            logger.debug("Failed to inspect __dict__: %s", e)

    # Fourth strategy: Check class CONTENTS for known attributes
    # This handles CDataFileContent classes like CPdbData
    if not result:
        try:
            contents = getattr(type(ccp4i2_object), 'CONTENTS', None)
            if contents:
                for attr_name in contents:
                    if hasattr(ccp4i2_object, attr_name):
                        attr_value = getattr(ccp4i2_object, attr_name)
                        if attr_value is not None:
                            result[attr_name] = value_dict_for_object(attr_value)
        except Exception as e:
            logger.debug("Failed to get CONTENTS attributes: %s", e)

    # Fifth strategy: Check known property names for file content classes
    if not result:
        known_content_props = ['composition', 'sequences', 'cell', 'spaceGroup']
        for prop_name in known_content_props:
            if hasattr(ccp4i2_object, prop_name):
                try:
                    prop_value = getattr(ccp4i2_object, prop_name)
                    if prop_value is not None:
                        result[prop_name] = value_dict_for_object(prop_value)
                except Exception:
                    pass

    return result if result else {}


def handle_plain_object(obj):
    """Convert a plain Python object with __dict__ to a dictionary."""
    logger.debug("Handling plain object: %s", type(obj).__name__)
    result = {}
    try:
        for attr_name, attr_value in obj.__dict__.items():
            if attr_name.startswith('_'):
                continue
            result[attr_name] = value_dict_for_object(attr_value)
    except Exception as e:
        logger.warning("Failed to convert plain object: %s", e)
    return result if result else {}


def handle_dict(ccp4i2_object):
    """Convert plain Python dict, recursively converting values."""
    logger.debug("Handling dict")
    result = {}
    for key in ccp4i2_object:
        result[key] = value_dict_for_object(ccp4i2_object[key])
    return result if result else {}


def handle_list(ccp4i2_object):
    """Convert plain Python list, recursively converting items."""
    logger.debug("Handling list")
    result = []
    for value in ccp4i2_object:
        result.append(value_dict_for_object(value))
    return result
