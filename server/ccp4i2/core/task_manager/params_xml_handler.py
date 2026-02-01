"""
Params XML Handler for CCP4i2 Task Parameters

This module handles .params.xml files that capture user-specified parameter values
for CCP4i2 tasks. It provides bidirectional functionality:

1. Export: Create .params.xml from a task hierarchy (only explicitly set parameters)
2. Import: Parse .params.xml and overlay values onto a task hierarchy

The params XML format uses simple tag names that correspond to the parameter names
in the task hierarchy, with values stored as text content or structured data for files.
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime
import socket
import os
import logging

from ..base_object.base_classes import CData, CContainer
from ..base_object.fundamental_types import *

logger = logging.getLogger(__name__)


class ParamsXmlHandler:
    """Handler for CCP4i2 .params.xml files."""

    def __init__(self):
        self.namespace = "http://www.ccp4.ac.uk/ccp4ns"
        self.namespace_prefix = "ccp4"

    def export_params_xml(
        self, task: CData, output_path: str, user_id: str = None, exclude_unset: bool = True
    ) -> bool:
        """
        Export explicitly set parameters from a task hierarchy to a .params.xml file.

        Args:
            task: The root task object (from .def.xml parsing)
            output_path: Path where to save the .params.xml file
            user_id: User ID for the header (defaults to current user)
            exclude_unset: If True, only export explicitly set parameters (default: True)

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Create root element with namespace
            root = ET.Element(f"{{{self.namespace}}}ccp4i2")
            root.set(f"xmlns:{self.namespace_prefix}", self.namespace)

            # Create header
            header = ET.SubElement(root, "ccp4i2_header")

            # Add header information
            ET.SubElement(header, "function").text = "PARAMS"
            ET.SubElement(header, "userId").text = user_id or os.getenv(
                "USER", "unknown"
            )
            ET.SubElement(header, "hostName").text = socket.gethostname()
            ET.SubElement(header, "creationTime").text = datetime.now().strftime(
                "%H:%M %d/%b/%y"
            )
            ET.SubElement(header, "pluginName").text = getattr(
                task, "name", "unknown_task"
            )

            # Create body - the container's CHILDREN go directly into ccp4i2_body
            # (not wrapped in a <container> element, per original CCP4i2 format)
            container = task.getEtree(excludeUnset=True)
            body = ET.Element('ccp4i2_body')
            # Append each child of the container directly to body (unwrap the container)
            for child in container:
                body.append(child)
            root.append(body)

            # Export all containers and their explicitly set parameters
            # No parent_context at top level (task is the root container)
            #self._export_container(task, body, parent_context=None, exclude_unset=exclude_unset)

            # Write to file with proper formatting
            self._write_formatted_xml(root, output_path)

            print(f"✅ Exported params to: {output_path}")
            return True

        except Exception as e:
            print(f"❌ Error exporting params XML: {e}")
            import traceback

            traceback.print_exc()
            return False

    def import_params_xml(self, task: CData, params_xml_path: str) -> bool:
        """
        Import parameter values from a .params.xml file and overlay them onto a task hierarchy.

        Args:
            task: The root task object (from .def.xml parsing)
            params_xml_path: Path to the .params.xml file

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            pass  # DEBUG: print(f"[DEBUG import_params_xml] UPDATED CODE - importing from: {params_xml_path}")
            if not Path(params_xml_path).exists():
                print(f"❌ Params file not found: {params_xml_path}")
                return False

            # Parse the XML
            tree = ET.parse(params_xml_path)
            root = tree.getroot()

            # Debug: show what we have
            logger.debug(f"[DEBUG] Root tag: {root.tag}")
            logger.debug(f"[DEBUG] Root children: {[child.tag for child in root]}")

            # Find the body element (handle both with and without namespace)
            # Check if root itself is the container (legacy format)
            if 'container' in root.tag.lower():
                logger.debug(f"[DEBUG] Legacy format detected - root is container")
                body = root
            else:
                # Try with namespace first
                body = root.find(".//{http://www.ccp4.ac.uk/ccp4ns}ccp4i2_body")
                if body is None:
                    # Try without namespace
                    body = root.find(".//ccp4i2_body")
                if body is None:
                    # Try legacy <body> tag
                    body = root.find(".//body")
                if body is None:
                    # Try direct child access
                    for child in root:
                        if 'body' in child.tag.lower():
                            body = child
                            logger.debug(f"[DEBUG] Found body via direct iteration: {child.tag}")
                            break

            if body is None:
                print("❌ No ccp4i2_body, body, or container found in params XML")
                logger.debug(f"[DEBUG] Searched in: {params_xml_path}")
                return False

            # Legacy format: <body><container><inputData>...
            # Modern format: <ccp4i2_body><inputData>...
            # If body has a <container> child, skip the wrapper
            container_elem = body.find("container")
            if container_elem is not None:
                logger.debug(f"[DEBUG] Skipping legacy <container> wrapper")
                body = container_elem

            # Import all parameter values
            imported_count = self._import_container_values(body, task)

            logger.debug("Imported %d parameters from: %s", imported_count, params_xml_path)
            return True

        except Exception as e:
            logger.error("Error importing params XML: %s", e, exc_info=True)
            return False

    def _export_container(self, container: CData, parent_elem: ET.Element, parent_context: str = None, exclude_unset: bool = True):
        """Recursively export container contents.

        Uses the new object-oriented getEtree() approach where each CData type
        (CContainer, CList, fundamental types) handles its own serialization.

        Args:
            container: CData container to export
            parent_elem: XML element to append to
            parent_context: Name of parent container (e.g., 'inputData', 'outputData') for context-aware export
            exclude_unset: If True, only export explicitly set parameters (default: True)
        """
        from ccp4i2.core.base_object.fundamental_types import CList

        # DEBUG: Print what container we're exporting
        container_name = getattr(container, 'name', 'unnamed')
        container_type = type(container).__name__
        pass # print(f"[DEBUG _export_container] Exporting container: {container_name} (type: {container_type})")

        # Get all attributes that are CData objects
        skip_attrs = [
            "CONTENTS",  # Class-level metadata, not an instance attribute
            "child_added",
            "child_removed",
            "content",  # For CDataFile - skip file content, only export metadata
            "fileContent",  # Alternative name for file content
            "destroyed",
            "finished",  # Skip signal objects
            "object_info",
            "parent_changed",
            "parent",  # Skip parent to avoid circular references
            "state",
        ]

        cdata_attrs_found = 0
        for attr_name in sorted(dir(container)):
            # Skip private, signal, and special attributes
            if attr_name.startswith("_") or attr_name in skip_attrs:
                continue

            # Try to get the attribute - may not exist or may be callable
            try:
                attr = getattr(container, attr_name)

                # Skip callables
                if callable(attr):
                    continue

                # Only export CData objects (not Signal, ObjectInfo, or other internal objects)
                if not isinstance(attr, CData):
                    continue

                cdata_attrs_found += 1

                if hasattr(attr, "name"):  # It's a CData object

                    if isinstance(attr, CList):
                        # Use CList's object-oriented getEtree() method
                        # CList.getEtree() handles serialization with proper class name tags
                        if hasattr(attr, 'getEtree'):
                            list_etree = attr.getEtree(excludeUnset=exclude_unset)
                            # Only append if the list has content
                            if list_etree.text or len(list_etree) > 0:
                                parent_elem.append(list_etree)

                    elif isinstance(attr, CContainer):
                        # Use CContainer's object-oriented getEtree() method
                        # CContainer.getEtree() handles recursion and serialization
                        container_elem = ET.SubElement(parent_elem, attr_name)

                        # Get the container's serialized content
                        # Note: we don't pass a name here because we already created the parent element
                        container_etree = attr.getEtree(name=None, excludeUnset=exclude_unset)

                        # Copy the children from the etree to our container element
                        for child in container_etree:
                            container_elem.append(child)

                        # Remove empty containers EXCEPT for standard sub-containers
                        # Legacy code expects inputData, outputData, and controlParameters
                        # to always be present, even if empty
                        standard_containers = {'inputData', 'outputData', 'controlParameters'}
                        if len(container_elem) == 0 and attr_name not in standard_containers:
                            parent_elem.remove(container_elem)

                    else:
                        # For simple CData types (CInt, CFloat, CString, CDataFile, etc.)
                        # Check if it should be exported based on context and value state
                        should_export = False
                        attr_type = type(attr).__name__

                        if parent_context == 'inputData':
                            # In inputData, export CDataFile objects if they have a baseName
                            try:
                                from ccp4i2.core.base_object.base_classes import CDataFile
                                if isinstance(attr, CDataFile):
                                    has_basename = hasattr(attr, 'baseName')
                                    basename_isset = attr.baseName.isSet() if has_basename else False
                                    basename_val = attr.baseName.value if (has_basename and basename_isset) else None

                                    if has_basename and basename_isset:
                                        if basename_val and str(basename_val).strip():
                                            should_export = True
                            except ImportError:
                                pass

                        # Default: check if explicitly set
                        if not should_export:
                            is_explicitly_set = self._is_explicitly_set(attr)
                            if is_explicitly_set:
                                should_export = True

                        if should_export:
                            param_elem = ET.SubElement(parent_elem, attr_name)
                            self._export_parameter_value(attr, param_elem)

            except Exception as e:
                print(f"Warning: Error processing {attr_name}: {e}")
                continue

        pass #print(f"[DEBUG _export_container] Finished exporting {container_name}: found {cdata_attrs_found} CData attributes")

    def _is_explicitly_set(self, param: CData) -> bool:
        """Check if a parameter has been explicitly set by the user."""
        try:
            # Special handling for CDataFile: check if baseName attribute is set
            # For CDataFiles, isSet is true if baseName is set AND not empty
            # Import here to avoid circular dependency
            try:
                from ccp4i2.core.base_object.base_classes import CDataFile
                if isinstance(param, CDataFile):
                    if hasattr(param, 'baseName'):
                        basename_obj = param.baseName
                        if hasattr(basename_obj, 'isSet') and basename_obj.isSet():
                            # Also check if the value is not empty
                            if hasattr(basename_obj, 'value'):
                                value = basename_obj.value
                                if isinstance(value, str) and len(value.strip()) > 0:
                                    return True
                    return False
            except ImportError:
                pass

            # Check if it has set state tracking
            if hasattr(param, "getValueState"):
                from ccp4i2.core.base_object.base_classes import ValueState

                state = param.getValueState("value")
                return state == ValueState.EXPLICITLY_SET

            # Fallback: check if it has a value attribute that's not default
            if hasattr(param, "value"):
                # For now, consider any non-empty value as explicitly set
                # This could be refined with better default tracking
                value = param.value
                if isinstance(value, bool):
                    return True  # Booleans are always considered set if present
                elif isinstance(value, (int, float)):
                    return True  # Numbers are always considered set if present
                elif isinstance(value, str):
                    return len(value.strip()) > 0
                else:
                    return value is not None

            return False

        except Exception:
            return False

    def _export_parameter_value(self, param: CData, elem: ET.Element):
        """Export the value of a parameter to an XML element."""
        try:
            param_name = getattr(param, 'name', 'unnamed')

            if hasattr(param, "value"):
                value = param.value

                # Handle different data types
                if isinstance(value, bool):
                    elem.text = str(value)
                elif isinstance(value, (int, float)):
                    elem.text = str(value)
                elif isinstance(value, str):
                    elem.text = value
                else:
                    elem.text = str(value)

            # Handle file objects (they have structured data)
            elif self._is_file_object(param):
                pass  # DEBUG: print(f"[DEBUG _export_parameter_value] Exporting file object '{param_name}', type: {type(param).__name__}")
                self._export_file_data(param, elem)

            # Handle performance indicators and other complex objects
            elif hasattr(param, "__dict__"):
                self._export_complex_object(param, elem)

        except Exception as e:
            param_name = param.objectName() if hasattr(param, 'objectName') else str(type(param).__name__)
            print(f"Warning: Error exporting value for {param_name}: {e}")
            elem.text = ""

    def _is_file_object(self, param: CData) -> bool:
        """Check if this is a file-type object."""
        type_name = type(param).__name__
        return "File" in type_name or "DataFile" in type_name

    def _export_file_data(self, file_obj: CData, elem: ET.Element):
        """Export file object data in the structured format."""
        file_name = getattr(file_obj, 'name', 'unnamed')
        pass  # DEBUG: print(f"[DEBUG _export_file_data] Exporting file '{file_name}'")

        # For file objects, we need to export structured data
        # Only export attributes that have been explicitly set (non-empty values)
        file_attrs = [
            "project",
            "baseName",
            "relPath",
            "annotation",
            "dbFileId",
            "subType",
            "contentFlag",
        ]

        exported_count = 0
        for attr in file_attrs:
            if hasattr(file_obj, attr):
                attr_obj = getattr(file_obj, attr)
                if attr_obj is not None:
                    # Get the string value
                    value_str = str(attr_obj)
                    pass  # DEBUG: print(f"[DEBUG _export_file_data]   {attr}: '{value_str}' (len={len(value_str.strip())})")

                    # Only export if:
                    # 1. It's a non-empty string, OR
                    # 2. It's a numeric value (subType, contentFlag)
                    # This prevents exporting empty <project />, <baseName />, etc.
                    if value_str and len(value_str.strip()) > 0:
                        sub_elem = ET.SubElement(elem, attr)
                        sub_elem.text = value_str
                        exported_count += 1
                    elif attr in ["subType", "contentFlag"]:
                        # Always export numeric flags even if 0
                        sub_elem = ET.SubElement(elem, attr)
                        sub_elem.text = value_str
                        exported_count += 1

        pass  # DEBUG: print(f"[DEBUG _export_file_data] Exported {exported_count} attributes for file '{file_name}'")

    def _export_complex_object(self, obj: CData, elem: ET.Element):
        """Export complex objects with multiple attributes."""
        # For objects like PERFORMANCEINDICATOR
        skip_attrs = [
            "CONTENTS",  # Class-level metadata
            "child_added",
            "child_removed",
            "destroyed",
            "finished",  # Skip signal objects
            "object_info",
            "parent_changed",
            "parent",  # Skip parent to avoid circular references
            "state",
            "name",
        ]
        for attr_name in dir(obj):
            if (
                not attr_name.startswith("_")
                and not callable(getattr(obj, attr_name))
                and attr_name not in skip_attrs
            ):

                try:
                    value = getattr(obj, attr_name)
                    if value is not None and not isinstance(value, CData):
                        sub_elem = ET.SubElement(elem, attr_name)
                        sub_elem.text = str(value)
                except Exception:
                    continue

    def _import_container_values(
        self, xml_container: ET.Element, cdata_container: CData
    ) -> int:
        """Recursively import values from XML into CData container."""
        from ccp4i2.core.base_object.fundamental_types import CList
        imported_count = 0

        # Special handling for CList: elements with class names indicate list items
        if isinstance(cdata_container, CList):
            # For CList, child elements are tagged with the item's class name (CCP4i2 convention)
            # e.g., <CImportUnmerged> for items in CImportUnmergedList
            # We create new items and populate them from the XML
            for child_elem in xml_container:
                elem_name = child_elem.tag

                # Create a new list item using makeItem()
                # This creates an item of the correct type based on the CList's subItem qualifier
                new_item = cdata_container.makeItem()

                # Import the item's contents - use appropriate method based on item type
                # For simple value types (CString, CInt, CFloat, CBoolean), the XML element
                # has text content directly (e.g., <CString>CA</CString>)
                # For complex types (CContainer), we need to recurse into child elements
                from ccp4i2.core.base_object.fundamental_types import CString, CInt, CFloat, CBoolean
                if isinstance(new_item, (CString, CInt, CFloat, CBoolean)):
                    # Simple value type - use _import_parameter_value which handles text content
                    self._import_parameter_value(child_elem, new_item)
                else:
                    # Complex type - recurse to handle nested structure
                    self._import_container_values(child_elem, new_item)

                # Append the item to the list
                cdata_container.append(new_item)
            return imported_count

        # Regular container handling
        for child_elem in xml_container:
            elem_name = child_elem.tag

            # Check if this corresponds to an attribute in the CData container
            if hasattr(cdata_container, elem_name):
                attr = getattr(cdata_container, elem_name)

                if isinstance(attr, CList):
                    # CList: Recurse to handle list items
                    imported_count += self._import_container_values(child_elem, attr)

                elif isinstance(attr, CContainer):
                    # Recurse into nested containers
                    imported_count += self._import_container_values(child_elem, attr)

                else:
                    # Import parameter value
                    if self._import_parameter_value(child_elem, attr):
                        imported_count += 1

        return imported_count

    def _import_parameter_value(self, xml_elem: ET.Element, param: CData) -> bool:
        """Import a value from XML element into a CData parameter."""
        try:
            # Handle simple values
            if xml_elem.text is not None and len(list(xml_elem)) == 0:
                # Simple text value
                text_value = xml_elem.text.strip()

                if hasattr(param, "value"):
                    # Convert to appropriate type
                    if isinstance(param, CBoolean):
                        param.value = text_value.lower() in ("true", "1", "yes")
                    elif isinstance(param, CInt):
                        param.value = int(text_value)
                    elif isinstance(param, CFloat):
                        param.value = float(text_value)
                    elif isinstance(param, CString):
                        param.value = text_value
                    else:
                        # Generic assignment
                        param.value = text_value

                    return True

            # Handle structured data (like file objects)
            elif len(list(xml_elem)) > 0:
                return self._import_structured_data(xml_elem, param)

            return False

        except Exception as e:
            param_name = param.objectName() if hasattr(param, 'objectName') else str(type(param).__name__)
            print(f"Warning: Error importing value for {param_name}: {e}")
            return False

    def _import_structured_data(self, xml_elem: ET.Element, param: CData) -> bool:
        """Import structured data from XML into a parameter."""
        from ccp4i2.core.CCP4File import CDataFile
        from ccp4i2.core.base_object.fundamental_types import CList
        try:
            imported_any = False

            # Special handling for CDataFile: use setFullPath() if baseName is present
            if isinstance(param, CDataFile):
                basename_elem = xml_elem.find('baseName')
                if basename_elem is not None and basename_elem.text:
                    basename_value = basename_elem.text.strip()
                    # Use setFullPath to properly initialize the file
                    param.setFullPath(basename_value)
                    imported_any = True
                # Continue to import other attributes (relPath, project, etc.)

            for child in xml_elem:
                attr_name = child.tag

                # Skip baseName if we already handled it above for CDataFile
                if isinstance(param, CDataFile) and attr_name == 'baseName':
                    continue

                # Check if this attribute exists on the param
                if not hasattr(param, attr_name):
                    continue

                # Get the attribute
                attr = getattr(param, attr_name)

                # Special handling for CList - need to use makeItem() and append()
                if isinstance(attr, CList):
                    self._import_container_values(child, attr)
                    imported_any = True
                    continue

                # If the child has nested elements, it's a CData object - recurse
                if len(list(child)) > 0:
                    if isinstance(attr, CData):
                        # Recursively import nested structure
                        imported_any |= self._import_structured_data(child, attr)
                    continue

                # Otherwise, it's a simple text value
                if child.text is not None:
                    value = child.text.strip()

                    # Skip importing contentFlag=0 from params.xml when DEF has set a non-zero value
                    # This prevents old params.xml files (with contentFlag=0) from overwriting
                    # the correct DEF default (e.g., contentFlag=2 for mmCIF files like CIFFILE)
                    # Legacy params.xml files always had contentFlag=0 because the DEF default
                    # wasn't being applied properly before the fix.
                    if attr_name == 'contentFlag' and value == '0':
                        existing_attr = getattr(param, attr_name, None)
                        if existing_attr is not None:
                            current_value = getattr(existing_attr, 'value', existing_attr)
                            if current_value and current_value != 0:
                                # Don't overwrite non-zero DEF default with zero from old params.xml
                                continue

                    # Set attribute if it exists
                    if hasattr(param, attr_name):
                        # Convert to appropriate type if needed
                        try:
                            existing_attr = getattr(param, attr_name)

                            # Use .set() method if available to ensure ValueState.EXPLICITLY_SET
                            # This is critical for parameter overlay - when we load from XML,
                            # we need to mark values as explicitly set so they get saved again
                            if hasattr(existing_attr, 'set'):
                                # It's a CData object with .set() method
                                if isinstance(existing_attr, bool):
                                    existing_attr.set(value.lower() in ("true", "1", "yes"))
                                elif hasattr(existing_attr, 'value'):
                                    # Check the type of the value to convert appropriately
                                    if isinstance(existing_attr.value, bool):
                                        existing_attr.set(value.lower() in ("true", "1", "yes"))
                                    elif isinstance(existing_attr.value, int):
                                        existing_attr.set(int(value))
                                    elif isinstance(existing_attr.value, float):
                                        existing_attr.set(float(value))
                                    else:
                                        existing_attr.set(value)
                                else:
                                    existing_attr.set(value)
                            else:
                                # Fallback to setattr for non-CData attributes
                                if isinstance(existing_attr, bool):
                                    setattr(param, attr_name, value.lower() in ("true", "1", "yes"))
                                elif isinstance(existing_attr, int):
                                    setattr(param, attr_name, int(value))
                                elif isinstance(existing_attr, float):
                                    setattr(param, attr_name, float(value))
                                else:
                                    setattr(param, attr_name, value)

                            imported_any = True

                        except (ValueError, TypeError):
                            # Just set as string if conversion fails
                            attr = getattr(param, attr_name)
                            if hasattr(attr, 'set'):
                                attr.set(value)
                            else:
                                setattr(param, attr_name, value)
                            imported_any = True

            return imported_any

        except Exception as e:
            print(f"Warning: Error importing structured data: {e}")
            return False

    def _write_formatted_xml(self, root: ET.Element, output_path: str):
        """Write XML with proper formatting and indentation."""
        # Create the tree
        tree = ET.ElementTree(root)

        # Add XML declaration and format
        self._indent_xml(root)

        # Write to file
        tree.write(output_path, encoding="utf-8", xml_declaration=True)

    def _indent_xml(self, elem: ET.Element, level: int = 0):
        """Add proper indentation to XML elements."""
        indent = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = indent + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = indent
            for child in elem:
                self._indent_xml(child, level + 1)
            if not child.tail or not child.tail.strip():
                child.tail = indent
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = indent


def export_task_params(task: CData, output_path: str, user_id: str = None) -> bool:
    """
    Convenience function to export task parameters to .params.xml file.

    Args:
        task: Task hierarchy (from .def.xml parsing)
        output_path: Where to save the .params.xml file
        user_id: User ID for the header

    Returns:
        bool: True if successful
    """
    handler = ParamsXmlHandler()
    return handler.export_params_xml(task, output_path, user_id)


def import_task_params(task: CData, params_xml_path: str) -> bool:
    """
    Convenience function to import parameters from .params.xml file into task.

    Args:
        task: Task hierarchy (from .def.xml parsing)
        params_xml_path: Path to .params.xml file

    Returns:
        bool: True if successful
    """
    handler = ParamsXmlHandler()
    return handler.import_params_xml(task, params_xml_path)


# Example usage
if __name__ == "__main__":
    # This would typically be used with a task loaded from .def.xml
    print("Params XML Handler - use export_task_params() and import_task_params()")
