import logging
import warnings
import re
from xml.etree import ElementTree as ET
from ccp4i2.core import CCP4Container
from ccp4i2.core import CCP4ModelData
from ccp4i2.core import CCP4Data
from ccp4i2.core.base_object.cdata import CData
from ccp4i2.core import CCP4File
from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.base_object.fundamental_types import CList

logger = logging.getLogger(f"ccp4x:{__name__}")


def find_objects(within, func, multiple=False, growing_list=None, growing_name=None):
    """
    Recursively searches for objects within a container or list that match a given condition.

    Args:
        within (CCP4Container.CContainer or list): The container or list to search within.
        func (callable): A function that takes an object and returns True if the object matches the condition.
        multiple (bool, optional): If True, find all matching objects. If False, stop after finding the first match. Defaults to True.
        growing_list (list, optional): A list to accumulate the matching objects. If None, a new list is created. Defaults to None.

    Returns:
        list: A list of objects that match the given condition.
    """
    if growing_list is None:
        growing_list = []
    original_length = len(growing_list)

    if growing_name is None:
        growing_name = within.objectPath()

    search_domain = (
        within._value
        if isinstance(within, (CCP4Data.CList, CList))
        else within.CONTENTS
    )
    for ichild, child_ref in enumerate(search_domain):
        child = (
            child_ref
            if isinstance(within, (CCP4Data.CList, CList))
            else getattr(within, child_ref, None)
        )
        current_name = (
            f"{growing_name}[{ichild}]"
            if isinstance(within._value, (CCP4Data.CList, CList, list))
            else f"{growing_name}.{child.objectName()}"
        )

        if func(child):
            growing_list.append(child)
            if not multiple:
                logger.debug("Match for %s", child.objectName())
                return growing_list
        elif isinstance(child, (CCP4Data.CList, CList, list)) or hasattr(
            child, "CONTENTS"
        ):
            find_objects(child, func, multiple, growing_list, current_name)
            if not multiple and len(growing_list) > original_length:
                return growing_list

    return growing_list


def find_object_by_path(base_element: CData, object_path: str):
    """
    Find a descendent CData item from a root element given a dot-separated path.

    DEPRECATED: Use base_element.find_by_path(object_path) instead.
    This function is now a thin wrapper around CContainer.find_by_path().

    Args:
        base_element: The container to search in (usually plugin.container)
        object_path: Dot-separated path (e.g., "prosmart_refmac.controlParameters.NCYCLES")
                    The first element (task name) is skipped since base_element is already
                    the container, not the plugin.

    Returns:
        The found CData object

    Raises:
        AttributeError: If the path is not found

    Example:
        >>> container = plugin.container
        >>> # Old way (deprecated):
        >>> obj = find_object_by_path(container, "prosmart_refmac.controlParameters.NCYCLES")
        >>> # New way (preferred):
        >>> obj = container.find_by_path("prosmart_refmac.controlParameters.NCYCLES")
    """
    warnings.warn(
        "find_object_by_path() is deprecated. Use container.find_by_path() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    # Delegate to the new core method
    if hasattr(base_element, 'find_by_path'):
        return base_element.find_by_path(object_path, skip_first=True)

    # Fallback for non-container objects (shouldn't happen, but be safe)
    path_elements = object_path.split(".")
    if len(path_elements) > 1:
        path_to_search = path_elements[1:]
    else:
        path_to_search = [object_path]

    current = base_element
    for segment in path_to_search:
        next_obj = getattr(current, segment, None)
        if next_obj is None and hasattr(current, 'find'):
            next_obj = current.find(segment)
        if next_obj is None:
            raise AttributeError(f"Element '{segment}' not found in '{current}'")
        current = next_obj

    return current
