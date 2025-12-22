"""
DEF XML Parser for CCP4i2 Task Definitions

This module parses .def.xml files to create dynamic CData object hierarchies
that represent task control parameters, input data, and output data structures.

The parser creates a complete in-memory representation using our modern CData
classes with hierarchical relationships, smart assignment, and set state tracking.
"""

import xml.etree.ElementTree as ET
from typing import Dict, Any, Optional, Union, Type
from pathlib import Path

from ..base_object.base_classes import CData, CContainer, ValueState
from ..base_object.fundamental_types import (
    CInt, CFloat, CBoolean, CString, CList
)
from ..base_object.metadata_system import (
    FieldMetadata, MetadataRegistry
)

# Import load_nested_xml for handling .def.xml inheritance
from ccp4i2.lib.utils.parameters.load_xml import load_nested_xml


class DefXmlParser:
    """Parser for CCP4i2 .def.xml task definition files."""

    def __init__(self):
        self.class_registry = self._build_class_registry()
        self.metadata_registry = MetadataRegistry()

    def _build_class_registry(self) -> Dict[str, Type[CData]]:
        """Build registry of available CData classes."""
        registry = {}

        # Add fundamental types
        registry.update(
            {
                "CInt": CInt,
                "CFloat": CFloat,
                "CBoolean": CBoolean,
                "CString": CString,
                "CContainer": CContainer,
                "CList": CList,
            }
        )

        # Add all implementation classes from core/
        import importlib
        import pkgutil

        try:
            # Import from the core package (implementation classes)
            from ccp4i2 import core

            # Get all implementation module files
            implementation_modules = [
                'CCP4Annotation',
                'CCP4ComFilePatchManager',
                'CCP4CootData',
                'CCP4CustomTaskManager',
                'CCP4Data',
                'CCP4File',
                'CCP4ImportedJobManager',
                'CCP4MathsData',
                'CCP4ModelData',
                'CCP4PerformanceData',
                'CCP4Preferences',
                'CCP4RefmacData',
                'CCP4XtalData',
            ]

            for module_name in implementation_modules:
                try:
                    module = importlib.import_module(f'ccp4i2.core.{module_name}')
                    for attr_name in dir(module):
                        if attr_name.startswith('_'):
                            continue
                        attr = getattr(module, attr_name)
                        if (
                            isinstance(attr, type)
                            and issubclass(attr, CData)
                            and attr is not CData
                            and attr is not CContainer
                        ):
                            registry[attr.__name__] = attr
                except ImportError as e:
                    print(f"Note: Could not import {module_name}: {e}")
                    continue

        except ImportError as e:
            print(f"Note: Implementation classes not yet available: {e}")

        return registry

    def parse_def_xml(self, xml_path: Union[str, Path]) -> CData:
        """
        Parse a .def.xml file and create the corresponding CData hierarchy.

        This method handles .def.xml inheritance by expanding <file> tags that reference
        parent .def.xml files (e.g., prosmart_refmac inheriting from refmac).

        Args:
            xml_path: Path to the .def.xml file

        Returns:
            Root CData object representing the task definition
        """
        xml_path = Path(xml_path)
        if not xml_path.exists():
            raise FileNotFoundError(f"DEF XML file not found: {xml_path}")

        # Parse XML
        tree = ET.parse(xml_path)
        root = tree.getroot()

        # Expand file references (inheritance) using load_nested_xml
        # This handles cases like prosmart_refmac inheriting from refmac
        root = load_nested_xml(root)

        # Extract task name from ccp4i2_body id or filename
        task_name = self._extract_task_name(root, xml_path)

        # Create root container
        root_container = CContainer()
        root_container._name = task_name

        # Parse the main ccp4i2_body structure
        body = root.find(".//ccp4i2_body[@id]")
        if body is not None:
            self._parse_body(body, root_container)

        # Parse any additional containers at root level (skip those already processed)
        all_containers = root.findall(".//container[@id]")
        processed = set()

        if body is not None:
            for nested_body in body.findall(".//ccp4i2_body"):
                for c in nested_body.findall("./container[@id]"):
                    processed.add(c.get("id"))

        for container in all_containers:
            container_id = container.get("id")
            if container_id not in processed:
                self._parse_container(container, root_container)

        # Re-enable validation after parsing is complete
        # This allows runtime validation when users set values
        self._reenable_validation(root_container)

        return root_container

    def _reenable_validation(self, obj: CData) -> None:
        """Recursively re-enable validation on all CData objects after parsing."""
        # Re-enable validation on this object
        if hasattr(obj, '_skip_validation'):
            obj._skip_validation = False

        # Recursively process all children in the hierarchy
        if hasattr(obj, 'children'):
            for child in obj.children():
                if isinstance(child, CData):
                    self._reenable_validation(child)

    def _extract_task_name(self, root: ET.Element, xml_path: Path) -> str:
        """Extract task name from XML structure or filename."""
        # Try to get from ccp4i2_body id
        body = root.find(".//ccp4i2_body[@id]")
        if body is not None:
            task_id = body.get("id")
            if task_id:
                return task_id

        # Try to get from pluginName in header
        plugin_name = root.find(".//pluginName")
        if plugin_name is not None and plugin_name.text:
            return plugin_name.text

        # Fall back to filename
        return xml_path.stem.replace(".def", "")

    def _parse_body(self, body: ET.Element, parent: CData) -> None:
        """Parse a ccp4i2_body element."""
        # Handle nested ccp4i2_body (some have nested structure)
        nested_body = body.find("./ccp4i2_body")
        if nested_body is not None:
            self._parse_body(nested_body, parent)
            return

        # Parse all containers in this body
        for container in body.findall("./container[@id]"):
            self._parse_container(container, parent)

    def _parse_container(self, container: ET.Element, parent: CData) -> None:
        """Parse a container element and add it to parent."""
        container_id = container.get("id")
        if not container_id:
            return

        # Create container object
        container_obj = CContainer()
        container_obj._name = container_id

        # Parse all content elements
        for content in container.findall("./content[@id]"):
            self._parse_content(content, container_obj)

        # Parse nested containers
        for nested_container in container.findall("./container[@id]"):
            self._parse_container(nested_container, container_obj)

        # Parse nested ccp4i2_body
        for nested_body in container.findall("./ccp4i2_body"):
            self._parse_body(nested_body, container_obj)

        # Add to parent - setattr will trigger hierarchy setup via __setattr__
        setattr(parent, container_id, container_obj)

    def _parse_content(self, content: ET.Element, parent: CData) -> None:
        """Parse a content element and add it to parent."""
        content_id = content.get("id")
        class_name = content.find("./className")

        if not content_id or class_name is None:
            return

        class_name_str = class_name.text
        if not class_name_str:
            return

        # Get qualifiers
        qualifiers = self._parse_qualifiers(content.find("./qualifiers"))

        # Create the appropriate object
        obj = self._create_object(class_name_str, qualifiers, content)

        if obj is not None:
            obj._name = content_id

            # Handle subItem for lists
            sub_item = content.find("./subItem")
            if sub_item is not None and isinstance(obj, CList):
                sub_class_name_elem = sub_item.find("./className")
                sub_qualifiers = self._parse_qualifiers(sub_item.find("./qualifiers"))
                if sub_class_name_elem is not None and sub_class_name_elem.text:
                    # Get the actual class object from registry
                    sub_class_name = sub_class_name_elem.text
                    sub_class = self.class_registry.get(sub_class_name)

                    # Set the subItem qualifier that makeItem() expects
                    if sub_class:
                        obj.set_qualifier('subItem', {
                            'class': sub_class,
                            'qualifiers': sub_qualifiers
                        })

            # Add to parent - setattr will trigger hierarchy setup via __setattr__
            if content_id == 'SMILESIN':
                import logging
                logger = logging.getLogger(__name__)
                logger.info(f"[SMILESIN DEBUG] About to setattr SMILESIN to {parent.objectName() if hasattr(parent, 'objectName') else parent}")
                logger.info(f"[SMILESIN DEBUG] obj={obj}, type={type(obj)}")
            setattr(parent, content_id, obj)
            if content_id == 'SMILESIN':
                import logging
                logger = logging.getLogger(__name__)
                logger.info(f"[SMILESIN DEBUG] After setattr, checking if SMILESIN is in parent's children...")
                if hasattr(parent, 'children'):
                    children_names = [c.objectName() for c in parent.children()]
                    logger.info(f"[SMILESIN DEBUG] children: {children_names}")
                    logger.info(f"[SMILESIN DEBUG] SMILESIN in children: {'SMILESIN' in children_names}")

    def _parse_qualifiers(self, qualifiers: Optional[ET.Element]) -> Dict[str, Any]:
        """Parse qualifiers element into a dictionary."""
        if qualifiers is None:
            return {}

        result = {}
        for child in qualifiers:
            tag = child.tag
            text = child.text

            # Handle different value types
            if text is None:
                result[tag] = None
            elif text.lower() in ("true", "false"):
                result[tag] = text.lower() == "true"
            elif tag in ("enumerators", "menuText", "fileExtensions"):
                # These are always arrays, even with a single value
                if "," in text:
                    result[tag] = [item.strip() for item in text.split(",")]
                else:
                    result[tag] = [text.strip()]
            elif tag in ("min", "max", "default") and self._is_number(text):
                result[tag] = self._parse_number(text)
            else:
                result[tag] = text

        # Handle special nested elements like <default>
        default_elem = qualifiers.find("./default")
        if default_elem is not None:
            default_dict = {}
            for child in default_elem:
                child_text = child.text
                if child_text is not None:
                    if self._is_number(child_text):
                        default_dict[child.tag] = self._parse_number(child_text)
                    else:
                        default_dict[child.tag] = child_text
            if default_dict:
                result["default"] = default_dict

        return result

    def _is_number(self, text: str) -> bool:
        """Check if text represents a number."""
        try:
            float(text)
            return True
        except (ValueError, TypeError):
            return False

    def _parse_number(self, text: str) -> Union[int, float]:
        """Parse a number string to int or float."""
        try:
            if "." in text:
                return float(text)
            else:
                return int(text)
        except (ValueError, TypeError):
            return text

    def _create_object(
        self, class_name: str, qualifiers: Dict[str, Any], content: ET.Element
    ) -> Optional[CData]:
        """Create an object of the specified class with given qualifiers."""
        # Handle special cases
        if class_name == "CList":
            return self._create_list_object(qualifiers, content)

        # Always construct CString for className 'CString'
        if class_name == "CString":
            obj = CString()
            # Skip validation during .def.xml parsing
            if hasattr(obj, '_skip_validation'):
                obj._skip_validation = True
            self._apply_qualifiers(obj, qualifiers, class_name)
            return obj

        # Get class from registry
        cls = self.class_registry.get(class_name)
        if cls is None:
            print(f"Warning: Unknown class '{class_name}', using CString as fallback")
            cls = CString

        # Create object
        try:
            obj = cls()

            # Skip validation during .def.xml parsing to avoid constraint violations
            # on default values (e.g., default=0.0 with min=1.0)
            if hasattr(obj, '_skip_validation'):
                obj._skip_validation = True

            # Apply qualifiers as metadata and initial values
            self._apply_qualifiers(obj, qualifiers, class_name)

            return obj

        except Exception as e:
            print(f"Error creating object of class '{class_name}': {e}")
            # Fallback to CString
            obj = CString()
            # Skip validation during .def.xml parsing
            if hasattr(obj, '_skip_validation'):
                obj._skip_validation = True
            self._apply_qualifiers(obj, qualifiers, "CString")
            return obj

    def _create_list_object(
        self, qualifiers: Dict[str, Any], content: ET.Element  # noqa: ARG002
    ) -> CList:
        """Create a CList object with proper item type."""
        # CList is already imported, just create an instance
        # qualifiers and content may be used in future for subItem type hints
        return CList()

    def _apply_qualifiers(
        self, obj: CData, qualifiers: Dict[str, Any], class_name: str
    ) -> None:
        """Apply qualifiers to an object."""
        # Create metadata for this field
        metadata = FieldMetadata(
            name=getattr(obj, "name", "unknown"), data_type=class_name
        )

        # Apply constraints and defaults
        for key, value in qualifiers.items():
            if key == "default":
                # Skip if value is the string "None" (from XML parsing)
                # This happens when .def.xml has <default>None</default>
                if value == "None":
                    continue

                if isinstance(value, dict):
                    # Complex default - store in metadata AND apply to object attributes
                    # This handles cases like CPdbDataFile with:
                    #   <default><subType>1</subType><contentFlag>2</contentFlag></default>
                    # These attributes need to be set on the object for fileExtensions() to work
                    metadata.default = value
                    for attr_name, attr_value in value.items():
                        if hasattr(obj, attr_name):
                            try:
                                attr = getattr(obj, attr_name)
                                if hasattr(attr, 'set'):
                                    attr.set(attr_value)
                                elif hasattr(attr, 'value'):
                                    attr.value = attr_value
                                else:
                                    setattr(obj, attr_name, attr_value)
                            except Exception as e:
                                print(f"Error setting default {attr_name}={attr_value}: {e}")
                else:
                    # Simple default - set the value and metadata
                    metadata.default = value
                    try:
                        if hasattr(obj, "_set_default_value"):
                            obj._set_default_value(value)
                        else:
                            # For fundamental types, set the actual value
                            if hasattr(obj, "value"):
                                obj.value = value
                                if hasattr(obj, "_value_states"):
                                    obj._value_states["value"] = ValueState.DEFAULT
                    except Exception as e:
                        # Silently skip errors for None-like values
                        if value not in [None, "", "None"]:
                            print(
                                f"Error setting default value {value} for {class_name}: {e}"
                            )

            elif key == "min":
                metadata.minimum = value
                # Also set on object for runtime validation
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("min", value)

            elif key == "max":
                metadata.maximum = value
                # Also set on object for runtime validation
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("max", value)

            elif key == "minLength":
                metadata.minlength = value
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("minLength", value)

            elif key == "maxLength":
                metadata.maxlength = value
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("maxLength", value)

            elif key == "enumerators":
                # Coerce string enumerators to numeric type for CInt/CFloat
                # This ensures validation comparisons work correctly
                if isinstance(value, list) and class_name in ("CInt", "CFloat"):
                    coerced_value = []
                    for item in value:
                        if isinstance(item, str):
                            try:
                                if class_name == "CInt":
                                    coerced_value.append(int(item))
                                else:  # CFloat
                                    coerced_value.append(float(item))
                            except (ValueError, TypeError):
                                coerced_value.append(item)  # Keep original if not convertible
                        else:
                            coerced_value.append(item)
                    value = coerced_value
                metadata.enumerators = value
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("enumerators", value)

            elif key == "menuText":
                metadata.menu_text = value

            elif key == "onlyEnumerators":
                metadata.only_enumerators = value
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("onlyEnumerators", value)

            elif key == "toolTip":
                metadata.help_text = value

            elif key == "mustExist":
                metadata.required = value
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("mustExist", value)

            elif key == "allowUndefined":
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("allowUndefined", value)

            elif key == "fromPreviousJob":
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("fromPreviousJob", value)

            elif key == "requiredSubType":
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("requiredSubType", value)

            elif key == "saveToDb":
                if hasattr(obj, "set_qualifier"):
                    obj.set_qualifier("saveToDb", value)

            # Store all qualifiers for reference
            if not hasattr(metadata, "extra_data"):
                metadata.extra_data = {}
            metadata.extra_data[key] = value

        # Store metadata on the object (simplified for now)
        if hasattr(obj, "name") and obj.name:
            obj._metadata = metadata


def parse_def_xml_file(xml_path: Union[str, Path]) -> CData:
    """
    Convenience function to parse a .def.xml file.

    Args:
        xml_path: Path to the .def.xml file

    Returns:
        Root CData object representing the task definition
    """
    parser = DefXmlParser()
    return parser.parse_def_xml(xml_path)


# Example usage and testing
if __name__ == "__main__":
    # Test with a sample XML string (partial for testing)
    sample_xml = """<?xml version="1.0" encoding="UTF-8"?>
    <ns0:ccp4i2 xmlns:ns0="http://www.ccp4.ac.uk/ccp4ns">
      <ccp4i2_body id="servalcat_pipe">
        <container id="inputData">
          <content id="XYZIN">
            <className>CPdbDataFile</className>
            <qualifiers>
              <ifAtomSelection>True</ifAtomSelection>
              <mustExist>True</mustExist>
              <allowUndefined>False</allowUndefined>
              <toolTip>File containing model coordinates (PDB/mmCIF).</toolTip>
            </qualifiers>
          </content>
        </container>
        <container id="controlParameters">
          <content id="NCYCLES">
            <className>CInt</className>
            <qualifiers>
              <default>10</default>
              <min>0</min>
              <toolTip>Number of refinement cycles to perform.</toolTip>
            </qualifiers>
          </content>
          <content id="ADD_WATERS">
            <className>CBoolean</className>
            <qualifiers>
              <default>False</default>
              <toolTip>Add waters and perform further refinement.</toolTip>
            </qualifiers>
          </content>
        </container>
      </ccp4i2_body>
    </ns0:ccp4i2>"""

    # Save sample and test
    with open("/tmp/test_def.xml", "w") as f:
        f.write(sample_xml)

    try:
        result = parse_def_xml_file("/tmp/test_def.xml")
        print("✅ DEF XML Parser created successfully!")
        print(f"Root object: {result}")
        print(f"Has inputData: {hasattr(result, 'inputData')}")
        print(f"Has controlParameters: {hasattr(result, 'controlParameters')}")

        if hasattr(result, "controlParameters"):
            ctrl = result.controlParameters
            if hasattr(ctrl, "NCYCLES"):
                print(f"NCYCLES value: {ctrl.NCYCLES.value}")
                print(f"NCYCLES is set: {ctrl.NCYCLES.isSet('value')}")
            if hasattr(ctrl, "ADD_WATERS"):
                print(f"ADD_WATERS value: {ctrl.ADD_WATERS.value}")
                print(f"ADD_WATERS is set: {ctrl.ADD_WATERS.isSet('value')}")

    except Exception as e:
        print(f"❌ Error testing DEF XML parser: {e}")
        import traceback

        traceback.print_exc()
