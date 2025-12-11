"""JSON encoder for CData container serialization."""
import json
import logging
from ccp4i2.core import CCP4Container
from ccp4i2.core import CCP4File
from ccp4i2.core.base_object import CData
from ccp4i2.core.base_object.fundamental_types import (
    CInt, CFloat, CString, CBoolean, CList
)


logger = logging.getLogger(f"ccp4i2:{__name__}")


def base_class(o):
    """Determine the base class name for a CData object."""
    if isinstance(o, CCP4File.CDataFile):
        result = "CDataFile"
    elif isinstance(o, CList):
        result = "CList"
    elif isinstance(o, CCP4Container.CContainer):
        result = "CContainer"
    elif isinstance(o, CString):
        result = "CString"
    elif isinstance(o, CInt):
        result = "CInt"
    elif isinstance(o, CFloat):
        result = "CFloat"
    elif isinstance(o, CBoolean):
        result = "CBoolean"
    else:
        result = "CData"
    return result


class CCP4i2JsonEncoder(json.JSONEncoder):
    """JSON encoder that serializes CData objects with full metadata."""

    def default(self, o):
        if isinstance(o, CData):
            qualifiers = {}
            # Safely get qualifiers from all possible sources
            # 1. Legacy QUALIFIERS class attribute (old style)
            if hasattr(type(o), 'QUALIFIERS'):
                qualifiers.update(type(o).QUALIFIERS)
            # 2. New metadata system: _class_qualifiers (@cdata_class)
            if hasattr(type(o), '_class_qualifiers'):
                class_quals = type(o)._class_qualifiers
                if isinstance(class_quals, dict):
                    qualifiers.update(class_quals)
            # 3. Instance-level _qualifiers (copied from class or overridden)
            if hasattr(o, '_qualifiers') and o._qualifiers:
                qualifiers.update(o._qualifiers)
            # Filter out NotImplemented values
            qualifiers = {
                k: v for k, v in qualifiers.items()
                if v is not NotImplemented
            }

            # Get CONTENTS_ORDER using dataOrder() - the single source of truth
            # for child ordering. This handles:
            # - Explicit CONTENTS_ORDER/CONTENT_ORDER attributes
            # - Auto-generated order from @cdata_class metadata
            # - Fallback to children() order
            contents_order = o.dataOrder() if hasattr(o, 'dataOrder') else []

            obj_path = ""
            if hasattr(o, 'objectPath'):
                obj_path = o.objectPath()

            # Determine _value based on object type
            # - Primitives (CInt, CFloat, CString, CBoolean): use _value attribute
            # - CList: use ordered array of children
            # - All other CData (CContainer, CDataFile, etc.): serialize children as dict
            #
            # Note: All CData objects can have children (via HierarchicalObject).
            # CDataFile and other @cdata_class decorated classes have pre-defined
            # children from their metadata. CContainer adds children programmatically.
            if isinstance(o, (CInt, CFloat, CString, CBoolean)):
                # Fundamental types - use _value
                value = getattr(o, '_value', None)
            elif isinstance(o, CList):
                # CList stores items in ordered array - they'll be serialized recursively
                value = list(o) if hasattr(o, '__iter__') else []
            else:
                # All other CData (CContainer, CDataFile, and any @cdata_class objects)
                # Build dict of children by name from children()
                children_by_name = {}
                for child in o.children():
                    if isinstance(child, CData):
                        child_name = child.objectName()
                        if child_name:
                            children_by_name[child_name] = child

                # Also check for any children defined in _data_order or CONTENT_ORDER
                # that might not be in children() yet (e.g., unset CDataFile objects)
                # This ensures the UI can display widgets for all defined parameters
                defined_children = set()
                if hasattr(o, '_data_order') and o._data_order:
                    defined_children.update(o._data_order)
                if hasattr(o, 'CONTENT_ORDER') and o.CONTENT_ORDER:
                    defined_children.update(o.CONTENT_ORDER)

                # For any defined child not yet in children_by_name, try to access it
                # via getattr which may trigger lazy initialization
                for name in defined_children:
                    if name not in children_by_name:
                        try:
                            child = getattr(o, name, None)
                            if child is not None and isinstance(child, CData):
                                children_by_name[name] = child
                        except Exception:
                            pass  # Skip any errors during lazy access

                # Use dataOrder() for complete ordering - it already handles
                # preferred order + remaining children
                value_dict = {}
                for name in contents_order:
                    if name in children_by_name:
                        value_dict[name] = children_by_name[name]

                value = value_dict

            result = {
                "_class": type(o).__name__,
                "_value": value,
                "_qualifiers": qualifiers,
                "_CONTENTS_ORDER": contents_order,
                "_objectPath": obj_path,
            }
            result["_baseClass"] = base_class(o)
            if isinstance(o, CList) and hasattr(o, 'makeItem'):
                result["_subItem"] = o.makeItem()
            return result
        if o is NotImplemented:
            return None
        if o.__class__.__name__ == "ObjectType":
            # Hack for container objects rooted as QObjects
            return None
        # Handle class/type objects (e.g., ABCMeta in columnGroupClassList)
        if isinstance(o, type):
            # Return class name as string for class references
            if o.__module__:
                return f"{o.__module__}.{o.__name__}"
            return o.__name__
        try:
            result = json.dumps(o)
        except TypeError:
            # Fall back to string representation without logging (reduces noise)
            result = str(o)
        return result
