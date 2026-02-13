"""
CContainer - Container class for heterogeneous collections of CData objects.

Provides container management with:
- addContent() / addObject() for child object creation
- deleteObject() for removal
- dataOrder() for ordered access
- children() for accessing all child items (from HierarchicalObject)
- XML serialization (loadContentsFromXml, saveContentsToXml)
- DEF file support (loadDefFile, saveDefFile)
- Params file support (loadParamsFile, saveParamsFile)

Supports both list-style and dict-style access:
- container[0] - Access by index (via children())
- container.itemName - Access by name (named attributes)
"""

from typing import List, Optional

from .cdata import CData
from .hierarchy_system import HierarchicalObject


class CContainer(CData):
    """Base class for container CData classes."""

    def __init__(self, parent=None, name=None, **kwargs):
        """Initialize a CContainer.

        Args:
            parent: Optional parent object
            name: Optional name for this container
            **kwargs: Additional arguments passed to parent CData class

        Note: Items should be added to containers using addContent(), addObject(),
        or direct attribute assignment. Access children via the children() method
        (inherited from HierarchicalObject) or by index (container[0]).
        """
        # Must call super().__init__() FIRST to set up base CData methods
        super().__init__(parent=parent, name=name, **kwargs)

        # Initialize _data_order from metadata contents_order if available
        if hasattr(self, '_metadata') and hasattr(self._metadata, 'contents_order'):
            contents_order = self._metadata.contents_order
            if contents_order:
                self._data_order = list(contents_order)  # Copy to avoid aliasing
            else:
                self._data_order = []
        else:
            self._data_order = []  # Track order of content items

    def addContent(self, content_class, name: str, **kwargs):
        """Add a new content item to the container (old API compatibility).

        Args:
            content_class: The class type to instantiate
            name: Name for the new content item
            **kwargs: Additional arguments to pass to the constructor

        Returns:
            The newly created content object
        """
        # Create instance of the content class WITHOUT parent (to avoid setattr issues)
        if isinstance(content_class, type):
            new_obj = content_class(name=name, **kwargs)
        else:
            # If it's a string, try to resolve it
            from .fundamental_types import CInt, CFloat, CBoolean, CString, CList
            class_map = {
                'CInt': CInt,
                'CFloat': CFloat,
                'CBoolean': CBoolean,
                'CString': CString,
                'CList': CList,
                'CContainer': CContainer,
            }
            cls = class_map.get(content_class)
            if cls is None:
                raise ValueError(f"Unknown content class: {content_class}")
            new_obj = cls(name=name, **kwargs)

        # Add to container
        setattr(self, name, new_obj)

        # Explicitly set parent relationship (in case setattr doesn't handle it)
        if hasattr(new_obj, 'set_parent'):
            new_obj.set_parent(self)

        self._data_order.append(name)
        return new_obj

    def addObject(self, obj: CData, name: str = None, reparent: bool = True):
        """Add an existing object to the container (old API compatibility).

        Args:
            obj: The CData object to add
            name: Optional name for the object (uses obj.objectName() if not provided)
            reparent: If True (default), set this container as the object's parent.
                     If False, preserve the object's existing parent relationship.
                     Legacy CCP4i2 compatibility parameter.

        Returns:
            The added object
        """
        if not isinstance(obj, CData):
            raise TypeError("Object must be a CData instance")

        obj_name = name if name is not None else obj.objectName()
        if not obj_name:
            raise ValueError("Object must have a name")

        # Set parent relationship (unless reparent=False)
        if reparent:
            obj.set_parent(self)
            obj._name = obj_name

        # Add to container
        setattr(self, obj_name, obj)
        if obj_name not in self._data_order:
            self._data_order.append(obj_name)

        return obj

    def deleteObject(self, name: str):
        """Delete an object from the container (old API compatibility).

        Args:
            name: Name of the object to delete

        Note: Silently succeeds if object doesn't exist (legacy code expects this behavior)
        """
        if not hasattr(self, name):
            return  # Silently succeed (legacy code expects defensive deletes)

        obj = getattr(self, name)

        # Cleanup hierarchy if it's a CData object
        if isinstance(obj, HierarchicalObject):
            try:
                obj.destroy()
            except Exception:
                pass

        # Remove from container
        delattr(self, name)

        # Remove from data order
        if name in self._data_order:
            self._data_order.remove(name)

    def dataOrder(self) -> list:
        """Return complete ordering of all children for serialization.

        Returns a list containing ALL child names, ordered as:
        1. Items from CONTENT_ORDER (if defined), filtered to actual children
        2. Items from _data_order not already included
        3. Any remaining children not in either list

        This ensures dataOrder() always returns a complete list of all children,
        with preferred ordering applied to the subset specified.

        Returns:
            List of all child names in serialization order
        """
        # Get all actual child names
        all_children = set()
        for child in self.children():
            if hasattr(child, 'objectName'):
                name = child.objectName()
                if name:
                    all_children.add(name)

        result = []
        seen = set()

        # 1. First, add items from CONTENT_ORDER (if defined)
        if hasattr(self, 'CONTENT_ORDER') and self.CONTENT_ORDER:
            for name in self.CONTENT_ORDER:
                if name in all_children and name not in seen:
                    result.append(name)
                    seen.add(name)

        # 2. Then add items from _data_order not already included
        if self._data_order:
            for name in self._data_order:
                if name in all_children and name not in seen:
                    result.append(name)
                    seen.add(name)

        # 3. Finally, add any remaining children not in either list
        for child in self.children():
            if hasattr(child, 'objectName'):
                name = child.objectName()
                if name and name not in seen:
                    result.append(name)
                    seen.add(name)

        return result

    @property
    def _dataOrder(self) -> list:
        """Legacy alias for _data_order (camelCase naming from old CCP4i2).

        The old CCP4i2 code directly accessed container._dataOrder. We use
        _data_order (snake_case), so provide this as a property alias.

        Returns:
            List of names in the order they were added
        """
        return self._data_order

    @property
    def CONTENTS_ORDER(self) -> list:
        """Legacy alias for _data_order (old CCP4i2 class attribute).

        In old CCP4i2, CONTENTS_ORDER was a class attribute defined in the
        metadata. In our new system, this is an instance attribute _data_order
        that tracks the order children were added.

        Used by legacy code like ccp4i2crank.py line 294.

        Returns:
            List of names in the order they were added
        """
        return self._data_order


    def copyData(self, otherContainer, dataList: Optional[List[str]] = None):
        """Copy data from another container into this container.

        This performs a deep copy by serializing the source data to XML (via getEtree)
        and deserializing it back into this container (via setEtree). This ensures
        that all nested structures, qualifiers, and metadata are properly copied.

        Args:
            otherContainer: Source CContainer to copy from
            dataList: Optional list of item names to copy. If None, copies all items
                     from the source container.

        Example:
            # Copy all data from source to destination
            dest.copyData(source)

            # Copy specific items only (as shown in pipelines.rst documentation)
            self.pdbset.container.inputData.copyData(self.container.inputData, ['XYZIN'])
        """
        if not isinstance(otherContainer, CContainer):
            raise TypeError(f"otherContainer must be a CContainer, got {type(otherContainer).__name__}")

        # Determine which items to copy
        if dataList is None:
            # Copy all items from source
            # Use dataOrder() if it has items, otherwise find all CData children
            items_from_order = otherContainer.dataOrder()
            if items_from_order:
                items_to_copy = items_from_order
            else:
                # Fallback: Find all CData attributes by iterating children
                # This handles cases where items weren't added to _data_order
                from ccp4i2.core.base_object.base_classes import CData
                items_to_copy = [child.objectName() for child in otherContainer.children()
                                if isinstance(child, CData)]
        else:
            items_to_copy = dataList

        # Copy each item using direct deep copy (no XML serialization)
        # This preserves value states from the source correctly
        for item_name in items_to_copy:
            if not hasattr(otherContainer, item_name):
                # Skip items that don't exist in source
                continue

            source_item = getattr(otherContainer, item_name)

            # Check if we already have this item in the destination
            if hasattr(self, item_name):
                # Item exists - deep copy values from source
                dest_item = getattr(self, item_name)
                if hasattr(dest_item, '_deep_copy_from'):
                    dest_item._deep_copy_from(source_item)
                else:
                    # Fallback: use smart assignment if _deep_copy_from not available
                    if hasattr(dest_item, '_smart_assign_from_cdata'):
                        dest_item._smart_assign_from_cdata(source_item)
            else:
                # Item doesn't exist - need to create it first
                # Get the class of the source item
                source_class = type(source_item)

                # Create a new instance
                new_item = source_class(name=item_name)

                # Add it to this container
                setattr(self, item_name, new_item)
                if hasattr(new_item, 'set_parent'):
                    new_item.set_parent(self)
                if item_name not in self._data_order:
                    self._data_order.append(item_name)

                # Now populate it with the source data using deep copy
                if hasattr(new_item, '_deep_copy_from'):
                    new_item._deep_copy_from(source_item)

    def clear(self):
        """Remove all content items from the container (old API compatibility)."""
        # Get list of all content items
        items_to_remove = list(self._data_order)

        # Delete each one
        for name in items_to_remove:
            try:
                self.deleteObject(name)
            except Exception:
                pass

        # Clear the order list
        self._data_order.clear()

    def loadContentsFromXml(self, xml_file: str):
        """Load container contents from an XML file (old API compatibility).

        Args:
            xml_file: Path to the XML file to load
        """
        import xml.etree.ElementTree as ET
        from pathlib import Path

        xml_path = Path(xml_file)
        if not xml_path.exists():
            raise FileNotFoundError(f"XML file not found: {xml_file}")

        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Use setEtree to deserialize
        self.setEtree(root, ignore_missing=True)

    def saveContentsToXml(self, xml_file: str):
        """Save container contents to an XML file (old API compatibility).

        Args:
            xml_file: Path to the XML file to save
        """
        import xml.etree.ElementTree as ET
        from pathlib import Path

        # Get the XML element tree
        root = self.getEtree()

        # Create the tree and write to file
        tree = ET.ElementTree(root)
        ET.indent(tree, space="  ")  # Pretty print with 2-space indentation
        tree.write(xml_file, encoding='utf-8', xml_declaration=True)

    def loadDataFromXml(self, xml_file: str):
        """Alias for loadContentsFromXml (old API compatibility)."""
        self.loadContentsFromXml(xml_file)

    def saveDataToXml(self, xml_file: str):
        """Alias for saveContentsToXml (old API compatibility)."""
        self.saveContentsToXml(xml_file)

    # Priority 3: DEF and PARAMS file-specific methods
    def loadDefFile(self, filename: str):
        """Load container structure from a .def.xml file (old API compatibility).

        DEF files define the structure and qualifiers of a container,
        but not the actual data values.

        Args:
            filename: Path to the .def.xml file
        """
        # For now, use loadContentsFromXml
        # In future, this could use the DEF XML parser specifically
        self.loadContentsFromXml(filename)

    def saveDefFile(self, filename: str):
        """Save container structure to a .def.xml file (old API compatibility).

        DEF files define the structure and qualifiers of a container,
        but not the actual data values.

        Args:
            filename: Path to the .def.xml file
        """
        # Save structure with qualifiers but without data values
        self.saveContentsToXml(filename)

    def loadParamsFile(self, filename: str):
        """Load container data values from a .params.xml file (old API compatibility).

        PARAMS files contain the actual data values for a container
        whose structure is already defined.

        Args:
            filename: Path to the .params.xml file
        """
        # Load data values into existing structure
        self.loadDataFromXml(filename)

    def saveParamsFile(self, filename: str):
        """Save container data values to a .params.xml file (old API compatibility).

        PARAMS files contain the actual data values for a container.

        Args:
            filename: Path to the .params.xml file
        """
        # Save only data values, not structure
        self.saveDataToXml(filename)

    def __bool__(self):
        """Return True for boolean conversion.

        Containers always evaluate to True when they exist, regardless of whether
        they contain items. This matches normal Python object behavior and supports
        legacy code patterns like:
            if container:  # Check if container exists

        To check if a container is empty, use:
            if len(container) == 0:
        """
        return True

    def __len__(self):
        """Return number of child items in container.

        Uses HierarchicalObject.children() to get the count.
        """
        return len(self.children())

    def __getitem__(self, index):
        """Get child item by index.

        Uses HierarchicalObject.children() to access children by index.
        Note: The order of children may not be deterministic if not added
        via named attributes tracked in _data_order.
        """
        return self.children()[index]

    def __setattr__(self, name: str, value):
        """Override setattr to maintain _data_order list.

        When a CData object is set as an attribute, add its name to _data_order
        to maintain the order of children for iteration.
        """
        # Call parent setattr first
        super().__setattr__(name, value)

        # If this is a CData child being added (not an internal attribute),
        # add it to _data_order if not already there
        if (not name.startswith('_')
            and hasattr(self, '_data_order')
            and isinstance(value, CData)
            and name not in self._data_order):
            self._data_order.append(name)

    def __getattr__(self, name: str):
        """Allow attribute-style access to children by name.

        This enables accessing child objects like container.inputData
        instead of having to search through children manually.

        Supports CCP4i2 truncation convention where long attribute names
        (e.g., RESOLUTION_HIGH) fall back to truncated versions (e.g., RESO_HIGH)
        if the full name is not found.

        Args:
            name: Name of the child to access

        Returns:
            The child object with matching name

        Raises:
            AttributeError: If no child with that name exists
        """
        # Avoid infinite recursion for internal attributes
        if name.startswith('_'):
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

        # Helper function to search for a child by name
        def find_child(search_name):
            # Check children from HierarchicalObject hierarchy
            try:
                children_method = object.__getattribute__(self, 'children')
                children_list = children_method()
                for child in children_list:
                    # Skip destroyed children
                    if hasattr(child, 'state'):
                        from .hierarchy_system import ObjectState
                        if child.state == ObjectState.DESTROYED:
                            continue
                    # Check if name matches
                    if hasattr(child, 'objectName') and child.objectName() == search_name:
                        return child
            except (AttributeError, TypeError):
                pass

            return None

        # First try the exact name
        result = find_child(name)
        if result is not None:
            return result

        # CCP4i2 truncation fallback: try truncating long keywords to 4 characters
        # Examples: RESOLUTION_HIGH -> RESO_HIGH, RESOLUTION_LOW -> RESO_LOW
        # Truncate ALL parts of compound names longer than 4 characters
        # Only apply to uppercase names with underscores (typical CCP4 parameters)
        if '_' in name and name.isupper():
            parts = name.split('_')
            # Truncate each part that's longer than 4 characters
            truncated_parts = [part[:4] if len(part) > 4 else part for part in parts]
            truncated_name = '_'.join(truncated_parts)
            # Only try if truncation actually changed something
            if truncated_name != name:
                result = find_child(truncated_name)
                if result is not None:
                    return result

        # CCP4i2 truncation fallback (REVERSE): try expanding short keywords
        # Examples: INPU_FIXED -> INPUT_FIXED, RESO_HIGH -> RESOLUTION_HIGH
        # For each 4-character part, try expanding to known common keywords
        if '_' in name and name.isupper():
            parts = name.split('_')
            # Common CCP4i2 keyword expansions (4-char -> full word)
            # These are from analyzing legacy .def.xml files
            expansions = {
                'RESO': 'RESOLUTION',
                'SEQU': 'SEQUENCE',
                'COMP': 'COMPOSITION',
                'SPAC': 'SPACEGROUP',
                'INPU': 'INPUT',
                'OUTP': 'OUTPUT',
                'ENSE': 'ENSEMBLE',
                'ATOM': 'ATOM',  # No change but included for consistency
                'CRYS': 'CRYSTAL',
                'SYMM': 'SYMMETRY',
                'REFM': 'REFMAC',
                # Add more as needed
            }

            # Try expanding each 4-character part
            expanded_parts = []
            any_expanded = False
            for part in parts:
                if len(part) == 4 and part in expansions:
                    expanded_parts.append(expansions[part])
                    any_expanded = True
                else:
                    expanded_parts.append(part)

            # Only try if at least one part was expanded
            if any_expanded:
                expanded_name = '_'.join(expanded_parts)
                result = find_child(expanded_name)
                if result is not None:
                    return result

        # Not found - raise AttributeError with debug info
        # List available children for debugging
        try:
            children_method = object.__getattribute__(self, 'children')
            available_names = [child.objectName() for child in children_method() if hasattr(child, 'objectName')]
            debug_info = f" Available children: {available_names}" if available_names else ""
        except:
            debug_info = ""
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'.{debug_info}")

    def find_by_path(self, path: str, skip_first: bool = True):
        """Find a descendant object by dot-separated path.

        Replaces server/ccp4i2/lib/utils/containers/find_objects.find_object_by_path()
        as a core CContainer method.

        Args:
            path: Dot-separated path (e.g., "controlParameters.NCYCLES")
                  Supports array indexing: "ASU_CONTENT[0].source"
            skip_first: If True (default), skip first path element for legacy compatibility.
                       Legacy paths often include task name as first element
                       (e.g., "prosmart_refmac.controlParameters.NCYCLES"), but since
                       we're called on the container, we skip it.

        Returns:
            The found CData object

        Raises:
            AttributeError: If the path is not found

        Example:
            >>> container = plugin.container
            >>> ncycles = container.find_by_path("prosmart_refmac.controlParameters.NCYCLES")
            >>> # Or with modern path (skip_first=False):
            >>> ncycles = container.find_by_path("controlParameters.NCYCLES", skip_first=False)
            >>> # Array indexing:
            >>> seq = container.find_by_path("task.inputData.ASU_CONTENT[0].source")
        """
        import re
        path_elements = path.split(".")

        # Skip first element if requested (legacy task/plugin name)
        if skip_first and len(path_elements) > 1:
            path_to_search = path_elements[1:]
        else:
            path_to_search = path_elements

        # Navigate the path using built-in hierarchy traversal
        current = self
        for segment in path_to_search:
            # Check for array indexing: "name[index]" or just "[index]"
            array_match = re.match(r'^(\w*)(\[(\d+)\])$', segment)
            if array_match:
                name_part = array_match.group(1)  # May be empty
                array_index = int(array_match.group(3))

                # First navigate to the named child (if name provided)
                if name_part:
                    next_obj = getattr(current, name_part, None)
                    if next_obj is None and hasattr(current, 'find'):
                        found = current.find(name_part)
                        if found != -1:
                            next_obj = found
                    if next_obj is None:
                        raise AttributeError(
                            f"Element '{name_part}' not found in path '{path}' "
                            f"(searching from {current.object_path()})"
                        )
                    current = next_obj

                # Then apply array indexing
                try:
                    # For CList, use __getitem__ which returns children by index
                    if hasattr(current, '__getitem__'):
                        current = current[array_index]
                    elif hasattr(current, 'children'):
                        children = current.children()
                        if 0 <= array_index < len(children):
                            current = children[array_index]
                        else:
                            raise IndexError(f"Index {array_index} out of range")
                    else:
                        raise AttributeError(
                            f"Object at '{segment}' does not support indexing"
                        )
                except (IndexError, KeyError) as e:
                    raise AttributeError(
                        f"Index [{array_index}] not valid in path '{path}' "
                        f"(searching from {current.object_path()}): {e}"
                    )
            else:
                # Regular name lookup (no array indexing)
                # Try getattr() first - triggers __getattr__ which searches children
                next_obj = getattr(current, segment, None)

                # Fall back to .find() if available - does depth-first recursive search
                # Note: find() returns -1 for not found (like Python's str.find())
                if next_obj is None and hasattr(current, 'find'):
                    found = current.find(segment)
                    # Only use the result if it's not -1 (not found indicator)
                    if found != -1:
                        next_obj = found

                if next_obj is None:
                    raise AttributeError(
                        f"Element '{segment}' not found in path '{path}' "
                        f"(searching from {current.object_path()})"
                    )

                current = next_obj

        return current

    def find_all_files(self):
        """Find all CDataFile objects in container hierarchy.

        Replaces server/ccp4i2/lib/cdata_utils.find_all_files() as a core method.

        Returns:
            List of all CDataFile objects found in hierarchy (deduplicated by id)

        Example:
            >>> output_files = plugin.outputData.find_all_files()
            >>> for file_obj in output_files:
            ...     print(f"Found {file_obj.objectName()}: {file_obj.object_path()}")
        """
        from .cdata_file import CDataFile

        files = []
        visited = set()  # Track visited object IDs to prevent duplicates/cycles

        def traverse(obj):
            """Recursively traverse hierarchy using children()"""
            obj_id = id(obj)
            if obj_id in visited:
                return
            visited.add(obj_id)

            # Check if this object is a file
            if isinstance(obj, CDataFile):
                files.append(obj)

            # Traverse children using standardized children() method
            if hasattr(obj, 'children'):
                try:
                    for child in obj.children():
                        if child is not None:
                            traverse(child)
                except Exception:
                    pass

        traverse(self)
        return files

    def set_parameter(self, object_path: str, value, skip_first: bool = True):
        """
        Set a parameter with automatic database awareness.

        This method combines the functionality of set_parameter_container() with
        automatic database synchronization when running in a CPluginScript context.

        If this container's parent chain includes a CPluginScript with _dbHandler,
        the parameter update will be automatically synchronized to the database.

        Args:
            object_path: Dot-separated path to parameter (e.g., "inputData.XYZIN")
            value: New value (str, int, dict, etc.)
            skip_first: If True (default), skip first path element for legacy compatibility

        Returns:
            The updated object

        Raises:
            AttributeError: If the path is not found
            Exception: If parameter setting fails

        Example:
            >>> # Database-independent usage
            >>> container.set_parameter("controlParameters.NCYCLES", 10)

            >>> # Database-aware usage (when container is part of CPluginScript)
            >>> plugin.container.set_parameter("inputData.XYZIN", "/path/to/file.pdb")
            >>> # Automatically registers file in database via plugin._dbHandler
        """
        import logging
        logger = logging.getLogger(__name__)

        # Navigate to the target object using find_by_path
        logger.debug(
            "Setting parameter %s = %s on container %s",
            object_path, value, self.object_path()
        )

        try:
            # Use modern find_by_path to navigate
            target_obj = self.find_by_path(object_path, skip_first=skip_first)
        except AttributeError as e:
            logger.error(
                "Failed to find parameter %s on container %s: %s",
                object_path, self.object_path(), str(e)
            )
            raise

        logger.debug(
            "Found target object for %s: type=%s, has_value=%s, has_set=%s, has_update=%s",
            object_path, type(target_obj).__name__,
            hasattr(target_obj, 'value'), hasattr(target_obj, 'set'), hasattr(target_obj, 'update')
        )

        # Set the value based on object type
        if hasattr(target_obj, 'value'):
            # It's a fundamental type (CInt, CFloat, CString, CBoolean)
            logger.debug("Setting via .value attribute")
            target_obj.value = value
        elif hasattr(target_obj, 'update') and isinstance(value, dict):
            # IMPORTANT: For dict values, use update() which only modifies specified keys.
            # This MUST come before the .set() check because:
            # - All CData objects have both .set() and .update()
            # - .set() has "set these fields, unset others" semantics
            # - .update() has "only update specified fields" semantics (sparse update)
            # Using .set() with a sparse dict would unset fields not in the dict!
            logger.debug("Setting via .update() method (sparse dict update)")
            target_obj.update(value)
        elif hasattr(target_obj, 'set'):
            # It's an object with a set() method (like CDataFile)
            # Used for non-dict values where we want full replacement semantics
            logger.debug("Setting via .set() method")
            target_obj.set(value)
        else:
            # Target object is NOT a CData type - this is a problem
            # We should only allow setting values on existing CData wrappers,
            # not creating new plain Python attributes on the container.
            #
            # If the target is a plain type (int, str, etc.), the .def.xml
            # structure was not properly loaded, or the path is wrong.
            raise AttributeError(
                f"Cannot set parameter '{object_path}': target is plain type "
                f"'{type(target_obj).__name__}', not a CData wrapper. "
                f"This usually means the parameter path is incorrect or the "
                f"plugin was not properly loaded from .def.xml."
            )

        logger.debug("Successfully set parameter %s to %s", object_path, value)

        # Check if we're in a database-aware context
        plugin_parent = self._find_plugin_parent()

        if plugin_parent and hasattr(plugin_parent, '_dbHandler') and plugin_parent._dbHandler:
            logger.debug(
                "Found CPluginScript parent with dbHandler, updating database for %s",
                object_path
            )

            # Update database via dbHandler
            try:
                plugin_parent._dbHandler.updateJobStatus(
                    jobId=str(plugin_parent._dbJobId),
                    container=plugin_parent.container
                )
                logger.debug("Successfully updated database for parameter %s", object_path)
            except Exception as e:
                logger.error(
                    "Failed to update database for parameter %s: %s",
                    object_path, str(e)
                )
                # Don't raise - parameter was set successfully even if DB update failed

        return target_obj

    def _find_plugin_parent(self):
        """
        Walk up the parent chain to find a CPluginScript instance.

        This enables automatic database awareness by detecting when a container
        is part of a CPluginScript hierarchy.

        Returns:
            CPluginScript instance if found, None otherwise

        Example:
            >>> plugin = CPluginScript()
            >>> plugin.container._find_plugin_parent()  # Returns plugin
            >>> standalone_container._find_plugin_parent()  # Returns None
        """
        # Import here to avoid circular dependency at module load time
        try:
            from ccp4i2.core.CCP4PluginScript import CPluginScript
        except ImportError:
            # If CPluginScript isn't available, we can't be in that context
            return None

        current = self.get_parent()
        while current:
            if isinstance(current, CPluginScript):
                return current
            current = current.get_parent()

        return None

    def validity(self):
        """
        Validate the container and all its children recursively.

        Delegates to CData.validity() which already recursively validates
        all children that have validity() methods. No additional child
        iteration needed here - doing so would cause duplicate errors.

        Returns:
            CErrorReport containing validation errors/warnings from entire hierarchy
        """
        return super().validity()

