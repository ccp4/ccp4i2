"""
XML Processing Module for Nested File References

This module provides functionality for processing XML documents that contain special
<file> nodes with embedded XML file references. It recursively copies XML structures
while expanding file references and merging their content into the destination document.

The module is designed to handle XML templates and configurations that reference
external XML files, commonly used in scientific workflow systems like CCP4i2.

Typical Usage Example:
    import xml.etree.ElementTree as ET
    from load_nested_xml import load_nested_xml

    # Load an XML file with embedded file references
    tree = ET.parse("input.xml")
    root = tree.getroot()

    # Process the XML, expanding file references and removing file nodes
    processed_root = load_nested_xml(root)

Key Features:
    - Removes <file> nodes while processing their referenced content
    - Merges content from referenced XML files into ccp4i2_body elements
    - Applies text overrides from simple nodes when values differ
    - Handles CCP4I2_TOP project references with relative path resolution
"""

import xml.etree.ElementTree as ET
import pathlib
import logging
from typing import Optional, Dict, List
from core import CCP4File


logger = logging.getLogger(f"ccp4x:{__name__}")


def load_nested_xml(src: ET.Element, dest: Optional[ET.Element] = None) -> ET.Element:
    """
    Copy an etree element to another etree element with special handling for 'file' nodes.

    This is the main entry point for the XML processing system. It recursively copies
    XML structures while applying special processing for <file> nodes that contain
    references to external XML files.

    File nodes are processed for their embedded content but are not included in the
    final XML output. Instead, their referenced content is merged into appropriate
    locations in the destination document.

    Args:
        src (ET.Element): Source etree element to copy from. This can be any XML element,
            but typically represents the root of an XML document.
        dest (ET.Element, optional): Destination etree element to copy to. If None,
            creates a new empty element with the same tag as src. Defaults to None.

    Returns:
        ET.Element: The destination etree element with copied content. All <file> nodes
            will have been removed, and their referenced content will have been merged
            into the appropriate locations.

    Example:
        >>> import xml.etree.ElementTree as ET
        >>> root = ET.fromstring('<root><child>text</child><file>...</file></root>')
        >>> result = load_nested_xml(root)
        >>> # result will contain <child>text</child> but not the <file> element

    Note:
        This function modifies the destination element in-place and also returns it.
        The source element is not modified.
    """
    if dest is None:
        dest = ET.Element(src.tag)

    # Copy attributes
    dest.attrib.update(src.attrib)

    # Copy text content
    if src.text:
        dest.text = src.text

    # Copy tail content (text after the element)
    if src.tail:
        dest.tail = src.tail

    # Handle special case for 'file' nodes - process them but don't add to destination
    if src.tag == "file":
        _handle_file_node(src, dest)
        # Return dest without adding the file node itself
        return dest

    # Recursively copy all child elements (excluding file nodes)
    for child in src:
        if child.tag != "file":
            child_copy = ET.Element(child.tag)
            child_copy = load_nested_xml(child, child_copy)

            # Check if this is a container with an id attribute that might need merging
            child_id = child_copy.get('id')
            if child_copy.tag == 'container' and child_id:
                # Look for existing container with same id in destination
                existing = None
                for dest_child in dest:
                    if dest_child.tag == 'container' and dest_child.get('id') == child_id:
                        existing = dest_child
                        break

                if existing is not None:
                    # Merge children from child_copy into existing container
                    child_count = len(list(child_copy))
                    logger.info(f"[load_nested_xml] Merging container[@id='{child_id}'] with {child_count} children into existing container")
                    for subchild in child_copy:
                        # Check if this is a content element that should be merged with existing
                        subchild_id = subchild.get('id')
                        if subchild.tag == 'content' and subchild_id:
                            existing_content = None
                            for ec in existing.findall(f"./content[@id='{subchild_id}']"):
                                existing_content = ec
                                break
                            if existing_content is not None:
                                # Merge qualifiers from subchild into existing content
                                _merge_content_qualifiers(existing_content, subchild)
                                logger.info(f"[load_nested_xml] Merged qualifiers for content[@id='{subchild_id}']")
                                continue  # Don't append, we merged
                        existing.append(subchild)
                else:
                    # No existing container with this id, append as new
                    logger.debug(f"[load_nested_xml] Appending new container[@id='{child_id}'] to destination")
                    dest.append(child_copy)
            else:
                # Not a container or no id, just append
                dest.append(child_copy)
        else:
            # Process file node but don't add it to destination
            _handle_file_node(child, dest)

    # After all processing is complete, apply text overrides from simple nodes
    _apply_text_overrides(src, dest)

    return dest


def _handle_file_node(file_node: ET.Element, dest_root: ET.Element) -> None:
    """
    Handle special processing for 'file' nodes with CI2XmlDataFile children.

    This function processes <file> nodes that contain CI2XmlDataFile elements, which
    specify references to external XML files. It extracts the file path information,
    loads the referenced XML file, and merges its content into the destination.

    The file node itself is NOT added to the destination - only its referenced content
    is processed and merged.

    Args:
        file_node (ET.Element): The <file> element to process. Expected to contain
            a CI2XmlDataFile child element with project, baseName, and relPath children.
        dest_root (ET.Element): The root destination element where ccp4i2_body children
            from the referenced file will be merged.

    Expected File Node Structure:
        <file>
            <CI2XmlDataFile>
                <project>CCP4I2_TOP</project>
                <relPath>pipelines/some_pipeline/script</relPath>
                <baseName>parameters.xml</baseName>
            </CI2XmlDataFile>
        </file>

    Note:
        Currently only supports project="CCP4I2_TOP". Other project values are ignored.
        The function logs status messages for debugging purposes.
    """
    # Find the CI2XmlDataFile child node
    ci2_xml_data_file = file_node.find("CI2XmlDataFile")
    if ci2_xml_data_file is None:
        logger.debug("File node does not contain CI2XmlDataFile child")
        return

    # Extract project, baseName, and relPath
    project_node = ci2_xml_data_file.find("project")
    base_name_node = ci2_xml_data_file.find("baseName")
    rel_path_node = ci2_xml_data_file.find("relPath")

    if project_node is None or base_name_node is None or rel_path_node is None:
        logger.debug(
            "CI2XmlDataFile missing required elements (project, baseName, or relPath)"
        )
        return

    # Get text values and strip whitespace
    project = project_node.text.strip() if project_node.text else ""
    base_name = base_name_node.text.strip() if base_name_node.text else ""
    rel_path = rel_path_node.text.strip() if rel_path_node.text else ""

    # Check if project is "CCP4I2_TOP"
    if project == "CCP4I2_TOP":
        # Construct the file path
        ccp4_file_parent = pathlib.Path(CCP4File.__file__).parent.parent
        file_path = ccp4_file_parent / rel_path / base_name

        logger.debug(
            f"Processing CCP4I2_TOP file path: {file_path} (file node will be removed)"
        )

        # Parse the XML file and merge ccp4i2_body children
        _parse_and_merge_xml_file(file_path, dest_root)
    else:
        logger.debug(f"Skipping file node with unsupported project: {project}")


def _parse_and_merge_xml_file(file_path: pathlib.Path, dest_root: ET.Element) -> None:
    """
    Parse an XML file and merge its ccp4i2_body children into the destination root.

    This function loads an external XML file, finds all ccp4i2_body elements within it,
    and recursively merges their children into the destination document's ccp4i2_body
    element. If no ccp4i2_body exists in the destination, one is created.

    Args:
        file_path (pathlib.Path): Path to the XML file to parse and merge.
        dest_root (ET.Element): The root destination element where ccp4i2_body
            children will be merged.

    Raises:
        ET.ParseError: If the XML file cannot be parsed due to syntax errors.
        FileNotFoundError: If the specified file does not exist.

    Note:
        This function logs detailed status messages for debugging.
        If the file doesn't exist or cannot be parsed, warnings/errors are logged but
        no exceptions are raised - the function fails gracefully.

    Processing Flow:
        1. Check if file exists
        2. Parse the XML file
        3. Find all ccp4i2_body elements in the parsed file
        4. Find or create ccp4i2_body in destination
        5. Recursively copy children using load_nested_xml (which handles nested files)
    """
    try:
        logger.debug(f"Attempting to parse XML file: {file_path}")

        # Check if file exists
        if not file_path.exists():
            logger.exception(f"File not found: {file_path}")
            return

        # Parse the XML file
        tree = ET.parse(file_path)
        parsed_root = tree.getroot()

        logger.debug(f"Successfully parsed XML. Root element: <{parsed_root.tag}>")

        # Find all ccp4i2_body nodes in the parsed XML
        ccp4i2_body_nodes = parsed_root.findall(".//ccp4i2_body")

        if not ccp4i2_body_nodes:
            logger.debug(f"No ccp4i2_body nodes found in {file_path}")
            return

        logger.debug(f"Found {len(ccp4i2_body_nodes)} ccp4i2_body node(s)")

        # Determine where to merge the included file's children
        # - If dest_root is ccp4i2_body: use it directly
        # - If dest_root is a container: merge directly into container (no nested ccp4i2_body)
        # - Otherwise (e.g., top-level ccp4i2): find or create ccp4i2_body
        if dest_root.tag == "ccp4i2_body":
            dest_ccp4i2_body = dest_root
            logger.debug("Destination root is already ccp4i2_body, using it directly")
        elif dest_root.tag == "container":
            # When including a file inside a container, merge directly into the container
            # This is critical for wrapper inclusion patterns like metalCoordWrapper
            dest_ccp4i2_body = dest_root
            logger.debug(f"Destination root is container[@id='{dest_root.get('id', '')}'], merging directly into it")
        else:
            dest_ccp4i2_body = _find_or_create_ccp4i2_body(dest_root)

        # Merge children from all ccp4i2_body nodes
        total_merged = 0
        total_appended = 0
        for i, body_node in enumerate(ccp4i2_body_nodes):
            children_count = len(list(body_node))
            logger.debug(
                f"Processing ccp4i2_body node {i+1} with {children_count} children"
            )

            # Copy all children of the ccp4i2_body node using our recursive function
            # This will also remove any nested file nodes while processing their content
            for child in body_node:
                # For file nodes, process them directly to merge content into dest_ccp4i2_body
                if child.tag == "file":
                    _handle_file_node(child, dest_ccp4i2_body)
                    continue  # File nodes are not added to destination, only their content

                # For non-file nodes, use load_nested_xml to recursively process
                child_copy = load_nested_xml(child)

                # Check if this is a container with an id attribute
                child_id = child_copy.get('id')
                if child_copy.tag == 'container' and child_id:
                    # Look for existing container with same id in destination
                    existing = None
                    for dest_child in dest_ccp4i2_body:
                        if dest_child.tag == 'container' and dest_child.get('id') == child_id:
                            existing = dest_child
                            break

                    if existing is not None:
                        # Merge children from child_copy into existing container
                        child_count = len(list(child_copy))
                        logger.info(f"Merging container[@id='{child_id}'] children ({child_count} children) into existing container")
                        for subchild in child_copy:
                            subchild_id = subchild.get('id', 'NO_ID')
                            logger.debug(f"  Merging child: <{subchild.tag}> id='{subchild_id}' into container[@id='{child_id}']")
                            existing.append(subchild)
                            total_merged += 1
                    else:
                        # No existing container with this id, append as new
                        dest_ccp4i2_body.append(child_copy)
                        total_appended += 1
                else:
                    # Not a container or no id, just append
                    dest_ccp4i2_body.append(child_copy)
                    total_appended += 1

        logger.debug(
            f"Successfully merged {total_merged} children into existing containers, appended {total_appended} new elements"
        )

    except ET.ParseError as e:
        logger.exception(f"Error parsing XML file {file_path}: {e}")
    except FileNotFoundError as e:
        logger.exception(f"File not found when parsing {file_path}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error processing {file_path}: {e}")


def _merge_content_qualifiers(existing_content: ET.Element, override_content: ET.Element) -> None:
    """
    Merge qualifier values from an override content element into an existing content element.

    This function handles the case where a parent def.xml overrides qualifiers of a content
    element from an included def.xml (via <file> node). The override content typically only
    specifies the qualifiers to override, not the full content definition.

    For example, servalcat_pipe.def.xml includes metalCoord.def.xml but overrides:
        <content id="XYZIN">
            <qualifiers>
                <allowUndefined>True</allowUndefined>
            </qualifiers>
        </content>

    This overrides the XYZIN's allowUndefined=False from metalCoord.def.xml.

    Args:
        existing_content: The content element from the included file (has className, full qualifiers)
        override_content: The override content element (may only have qualifiers to override)

    Algorithm:
        1. Find <qualifiers> element in both existing and override content
        2. If override has no qualifiers, nothing to merge
        3. If existing has no qualifiers but override does, copy entire qualifiers element
        4. Otherwise, merge individual qualifier children from override into existing
           (override values take precedence)
    """
    # Find qualifiers elements
    existing_qualifiers = existing_content.find('qualifiers')
    override_qualifiers = override_content.find('qualifiers')

    # Nothing to merge if override has no qualifiers
    if override_qualifiers is None:
        return

    # If existing has no qualifiers, copy the entire qualifiers element from override
    if existing_qualifiers is None:
        existing_content.append(override_qualifiers)
        logger.debug(f"Added new qualifiers element to content[@id='{existing_content.get('id')}']")
        return

    # Merge individual qualifier values from override into existing
    for override_qualifier in override_qualifiers:
        qualifier_tag = override_qualifier.tag

        # Find matching qualifier in existing
        existing_qualifier = existing_qualifiers.find(qualifier_tag)

        if existing_qualifier is not None:
            # Override the existing value
            old_value = existing_qualifier.text
            existing_qualifier.text = override_qualifier.text
            logger.debug(
                f"Override qualifier '{qualifier_tag}' in content[@id='{existing_content.get('id')}']: "
                f"'{old_value}' -> '{override_qualifier.text}'"
            )
        else:
            # Add the new qualifier
            existing_qualifiers.append(override_qualifier)
            logger.debug(
                f"Added qualifier '{qualifier_tag}' to content[@id='{existing_content.get('id')}']: "
                f"'{override_qualifier.text}'"
            )


def _apply_text_overrides(src: ET.Element, dest: ET.Element) -> None:
    """
    Apply text content overrides from simple nodes in src to matching nodes in dest.

    Matching is based on both tag name and 'id' attribute values throughout the path.
    This ensures that only truly equivalent elements are overridden, not just elements
    with similar xpath locations.

    Args:
        src (ET.Element): Source element to extract simple node text from.
        dest (ET.Element): Destination element to apply text overrides to.

    Algorithm:
        1. Find all simple nodes with text content in source (excluding file nodes)
        2. Build tag+id path-to-element mapping for destination
        3. For each simple node, find matching tag+id path in destination
        4. Compare text values and apply override only if they differ

    Example:
        Source: <container id="input"><content id="file1">parent_value</content></container>
        Dest:   <container id="input"><content id="file1">embedded_value</content></container>
        Result: The content will be overridden because both tag and id match.

        But if dest had <content id="file2">, no override would occur.
    """
    # Get all simple nodes (nodes without children) from source that have text content
    # Exclude simple nodes that are within file nodes
    simple_nodes_with_text = _get_simple_nodes_with_id_path(
        src, exclude_file_nodes=True
    )

    if not simple_nodes_with_text:
        return

    logger.debug(
        f"Found {len(simple_nodes_with_text)} simple nodes with text content to check for overrides (excluding file nodes)"
    )

    # Build tag+id path to element mapping for destination
    dest_path_map = _build_id_path_map(dest)

    # Apply overrides only when text differs
    overrides_applied = 0
    matches_checked = 0

    for id_path, src_text_content in simple_nodes_with_text.items():
        if id_path in dest_path_map:
            dest_elements = dest_path_map[id_path]
            for dest_element in dest_elements:
                matches_checked += 1
                dest_text = dest_element.text.strip() if dest_element.text else ""

                # Only apply override if text values differ
                if dest_text != src_text_content:
                    dest_element.text = src_text_content
                    logger.debug(
                        f"Override applied: {id_path} '{dest_text}' -> '{src_text_content}'"
                    )
                    overrides_applied += 1

    if matches_checked > 0:
        logger.debug(
            f"Checked {matches_checked} matching path(s), applied {overrides_applied} text overrides"
        )
    else:
        logger.debug("No matching tag+id paths found for text overrides")


def _get_simple_nodes_with_id_path(
    element: ET.Element, current_path: str = "", exclude_file_nodes: bool = False
) -> Dict[str, str]:
    """
    Get all simple nodes (leaf nodes with text content) using tag+id path identification.

    A "simple node" is defined as an element that has no child elements but contains
    text content. The path is built using both tag names and id attributes.

    Args:
        element (ET.Element): Element to traverse for simple nodes.
        current_path (str, optional): Current tag+id path being built during traversal.
            Used internally for recursion. Defaults to "".
        exclude_file_nodes (bool, optional): If True, exclude nodes that are within
            <file> elements from the results. Defaults to False.

    Returns:
        Dict[str, str]: Dictionary mapping tag+id path strings to text content.
        Each path uniquely identifies a simple node location based on both tag names
        and id attributes, and the value is the stripped text content.

    Path Format:
        "tag1[@id=value1]/tag2[@id=value2]/tag3[@id=value3]"
        If an element has no id attribute, it's represented as "tag[@id=]"

    Example:
        >>> xml = '<root><container id="main"><content id="file1">text1</content></container></root>'
        >>> element = ET.fromstring(xml)
        >>> result = _get_simple_nodes_with_id_path(element)
        >>> print(result)
        {'root[@id=]/container[@id=main]/content[@id=file1]': 'text1'}
    """
    simple_nodes = {}

    # Build current tag+id path
    element_id = element.get("id", "")
    element_path = f"{element.tag}[@id={element_id}]"

    if current_path:
        id_path = f"{current_path}/{element_path}"
    else:
        id_path = element_path

    # Skip file nodes and their descendants if exclude_file_nodes is True
    if exclude_file_nodes and element.tag == "file":
        return simple_nodes

    # Check if this is a simple node (no children) with text content
    children = list(element)
    if not children and element.text and element.text.strip():
        simple_nodes[id_path] = element.text.strip()

    # Recursively process children
    for child in children:
        child_simple_nodes = _get_simple_nodes_with_id_path(
            child, id_path, exclude_file_nodes
        )
        simple_nodes.update(child_simple_nodes)

    return simple_nodes


def _build_id_path_map(
    element: ET.Element, current_path: str = ""
) -> Dict[str, List[ET.Element]]:
    """
    Build a mapping from tag+id path strings to lists of elements at those paths.

    This function creates a comprehensive index of all elements in an XML tree,
    organized by their tag+id path location. This enables efficient lookup of elements
    by their tag and id combination for operations like text overrides.

    Args:
        element (ET.Element): Element to traverse and map.
        current_path (str, optional): Current tag+id path being built during traversal.
            Used internally for recursion. Defaults to "".

    Returns:
        Dict[str, List[ET.Element]]: Dictionary mapping tag+id path strings to lists of
            elements found at those paths. Multiple elements can theoretically share the
            same path if there are duplicates, though this would be unusual.

    Path Format:
        "tag1[@id=value1]/tag2[@id=value2]/tag3[@id=value3]"
        If an element has no id attribute, it's represented as "tag[@id=]"

    Example:
        >>> xml = '<root><container id="main"><content id="file1">text</content></container></root>'
        >>> element = ET.fromstring(xml)
        >>> path_map = _build_id_path_map(element)
        >>> path = 'root[@id=]/container[@id=main]/content[@id=file1]'
        >>> len(path_map[path])
        1
    """
    path_map = {}

    # Build current tag+id path
    element_id = element.get("id", "")
    element_path = f"{element.tag}[@id={element_id}]"

    if current_path:
        id_path = f"{current_path}/{element_path}"
    else:
        id_path = element_path

    # Add current element to map
    if id_path not in path_map:
        path_map[id_path] = []
    path_map[id_path].append(element)

    # Recursively process children
    for child in element:
        child_path_map = _build_id_path_map(child, id_path)
        # Merge child maps
        for child_id_path, child_elements in child_path_map.items():
            if child_id_path not in path_map:
                path_map[child_id_path] = []
            path_map[child_id_path].extend(child_elements)

    return path_map


def _find_or_create_ccp4i2_body(root: ET.Element) -> ET.Element:
    """
    Find an existing ccp4i2_body element in the root, or create one if it doesn't exist.

    The ccp4i2_body element is a special container used in CCP4i2 XML documents to
    hold the main content. This function ensures that such an element exists in the
    destination document for content merging operations.

    Args:
        root (ET.Element): The root element to search in for ccp4i2_body.

    Returns:
        ET.Element: The ccp4i2_body element, either found or newly created.

    Search Strategy:
        Uses ".//ccp4i2_body" xpath to find ccp4i2_body elements anywhere in the
        document tree, not just as direct children of root.

    Creation Strategy:
        If no ccp4i2_body is found, creates one as a direct child of the root element.

    Note:
        This function logs status messages indicating whether an existing element
        was found or a new one was created.
    """
    # Try to find existing ccp4i2_body
    ccp4i2_body = root.find(".//ccp4i2_body")

    if ccp4i2_body is not None:
        logger.debug("Found existing ccp4i2_body in destination")
        return ccp4i2_body

    # Create new ccp4i2_body if not found
    logger.debug("Creating new ccp4i2_body in destination")
    ccp4i2_body = ET.SubElement(root, "ccp4i2_body")
    return ccp4i2_body


# Example usage and testing
if __name__ == "__main__":
    """
    Example usage demonstrating the load_nested_xml functionality.

    This example loads a real CCP4i2 pipeline definition file and processes it
    to show how file nodes are expanded and content is merged. It provides
    detailed output for debugging and verification purposes.
    """
    # Load XML from the prosmart_refmac.def.xml file
    xml_file_path = (
        pathlib.Path(CCP4File.__file__).parent.parent
        / "pipelines"
        / "prosmart_refmac"
        / "script"
        / "prosmart_refmac.def.xml"
    )

    try:
        print(f"Loading XML from: {xml_file_path}")

        # Check if file exists
        if not xml_file_path.exists():
            print(f"Error: XML file not found at {xml_file_path}")
            exit(1)

        # Parse the XML file
        tree = ET.parse(xml_file_path)
        root = tree.getroot()

        print(f"Successfully loaded XML. Root element: <{root.tag}>")
        print(f"Root attributes: {root.attrib}")

        # Count file nodes in the original XML
        original_file_nodes = root.findall(".//file")
        print(f"Original XML contains {len(original_file_nodes)} 'file' nodes")

        # Test the load_nested_xml function
        print("\n" + "=" * 60)
        print("Testing load_nested_xml with prosmart_refmac.def.xml:")
        print("=" * 60)

        copied_root = load_nested_xml(root)

        print(f"\nCopied successfully. Copied root element: <{copied_root.tag}>")
        print(f"Number of direct children in original: {len(list(root))}")
        print(f"Number of direct children in copy: {len(list(copied_root))}")

        # Verify file nodes have been removed
        result_file_nodes = copied_root.findall(".//file")
        print(f"\nResult contains {len(result_file_nodes)} 'file' nodes (should be 0)")

        # Check for ccp4i2_body in the result
        result_body_nodes = copied_root.findall(".//ccp4i2_body")
        if result_body_nodes:
            print(f"Result contains {len(result_body_nodes)} ccp4i2_body node(s)")
            for i, body_node in enumerate(result_body_nodes):
                children_count = len(list(body_node))
                print(f"  ccp4i2_body {i+1} has {children_count} children")
        else:
            print(f"No ccp4i2_body nodes in result")

        # Show simple nodes that would be candidates for text override
        simple_nodes = _get_simple_nodes_with_text(root, exclude_file_nodes=True)
        if simple_nodes:
            print(
                f"\nFound {len(simple_nodes)} simple nodes with text in original (excluding file nodes):"
            )
            for xpath, text in list(simple_nodes.items())[:10]:  # Show first 10
                print(f"  {xpath}: '{text}'")
            if len(simple_nodes) > 10:
                print(f"  ... and {len(simple_nodes) - 10} more")

        # Optionally print a sample of the XML structure (first few lines)
        xml_string = ET.tostring(copied_root, encoding="unicode")
        lines = xml_string.split("\n")
        print(f"\nFirst 50 lines of expanded XML structure (file nodes removed):")
        for i, line in enumerate(lines[:50]):
            print(f"{i+1:2d}: {line}")
        if len(lines) > 50:
            print(f"... ({len(lines) - 50} more lines)")

    except ET.ParseError as e:
        print(f"Error parsing XML file: {e}")
    except FileNotFoundError as e:
        print(f"File not found: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")
