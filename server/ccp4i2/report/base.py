"""
Modern base classes for CCP4i2 report elements.

This module provides the foundational classes for the report system,
built on top of HierarchicalObject from base_object. These classes
are designed to be backward-compatible with existing Report subclasses
in the legacy ccp4i2 plugins while enabling modern Python patterns.

Key classes:
    - ReportElement: Base class for all report elements
    - ReportContainer: Base class for elements that contain children

The design goals are:
    a) Robust error handling with structured diagnostics
    b) Modern Python patterns (type hints, dataclasses where appropriate)
    c) Parsimonious code using base_object properties
    d) Well documented with clear API
    e) Extensible for grid layout and future features

Backward Compatibility:
    The existing Report class in CCP4ReportParser.py will be gradually
    migrated to inherit from ReportContainer. Legacy plugin code that
    subclasses Report will continue to work unchanged.
"""

from typing import Optional, Dict, Any, List, Iterator, Type, Callable, TYPE_CHECKING
from xml.etree import ElementTree as ET
import logging

from ccp4i2.core.base_object.hierarchy_system import HierarchicalObject
from ccp4i2.report.errors import DiagnosticCollector, ReportDiagnostic, DiagnosticLevel

if TYPE_CHECKING:
    from ccp4i2.report.CCP4ReportParser import Report

logger = logging.getLogger(f"ccp4i2:{__name__}")


class ReportElement(HierarchicalObject):
    """
    Base class for all report elements using modern base_object patterns.

    This class integrates with HierarchicalObject to provide:
    - Parent-child relationships with weak references
    - Object lifecycle management
    - Path-based navigation

    Attributes:
        id: Unique identifier for this element
        class_: CSS class(es) for styling
        style: Inline CSS styles
        label: Human-readable label
        title: Title/tooltip text
        diagnostics: Collector for errors and warnings

    Example:
        class MyElement(ReportElement):
            def __init__(self, value: str, **kwargs):
                super().__init__(**kwargs)
                self.value = value

            def as_data_etree(self) -> ET.Element:
                root = super().as_data_etree()
                root.set('value', self.value)
                return root
    """

    # Class-level counter for generating unique IDs
    _element_counter: int = 0

    def __init__(
        self,
        id: Optional[str] = None,
        class_: Optional[str] = None,
        style: Optional[str] = None,
        label: Optional[str] = None,
        title: Optional[str] = None,
        parent: Optional['ReportElement'] = None,
        **kwargs
    ):
        """
        Initialize a report element.

        Args:
            id: Unique identifier (auto-generated if not provided)
            class_: CSS class(es) for styling
            style: Inline CSS styles
            label: Human-readable label
            title: Title/tooltip text
            parent: Parent element in the hierarchy
            **kwargs: Additional attributes passed to HierarchicalObject
        """
        # Generate ID if not provided
        if id is None:
            id = self._generate_id()

        # Initialize HierarchicalObject with name=id
        super().__init__(name=id, parent=parent, **kwargs)

        self.id = id
        self.class_ = class_
        self.style = style
        self.label = label
        self.title = title

        # Diagnostics collector - each element can collect its own
        self._diagnostics = DiagnosticCollector()

        # For backward compatibility with legacy code
        self.internalId = id
        self.text = ''
        self.tail = ''

    @classmethod
    def _generate_id(cls) -> str:
        """Generate a unique ID for this element."""
        cls._element_counter += 1
        return f"{cls.__name__}_{cls._element_counter}"

    @property
    def diagnostics(self) -> DiagnosticCollector:
        """Get the diagnostics collector for this element."""
        return self._diagnostics

    def add_diagnostic(
        self,
        level: DiagnosticLevel,
        code: str,
        message: str,
        **kwargs
    ) -> ReportDiagnostic:
        """
        Add a diagnostic to this element.

        Convenience method that automatically sets the location.
        """
        location = kwargs.pop('location', f"{type(self).__name__}")
        return self._diagnostics.add(level, code, message, location=location, **kwargs)

    def warning(self, code: str, message: str, **kwargs) -> ReportDiagnostic:
        """Add a warning diagnostic."""
        return self.add_diagnostic(DiagnosticLevel.WARNING, code, message, **kwargs)

    def error(self, code: str, message: str, **kwargs) -> ReportDiagnostic:
        """Add an error diagnostic."""
        return self.add_diagnostic(DiagnosticLevel.ERROR, code, message, **kwargs)

    def get_report(self) -> Optional['Report']:
        """
        Get the root Report object.

        Traverses up the parent hierarchy to find the Report.
        Backward compatible with legacy getReport() method.
        """
        if hasattr(self, 'parent') and self.parent is not None:
            if hasattr(self.parent, 'get_report'):
                return self.parent.get_report()
            elif hasattr(self.parent, 'getReport'):
                return self.parent.getReport()
        return None

    # Backward compatibility alias
    def getReport(self) -> Optional['Report']:
        """Legacy alias for get_report()."""
        return self.get_report()

    def data_id(self) -> str:
        """Get the data ID for this element (used in XML output)."""
        return f'data_{self.internalId}'

    def data_url(self) -> str:
        """Get the data URL for this element."""
        return f'./tables_as_xml_files/{self.data_id()}.xml'

    def as_data_etree(self) -> ET.Element:
        """
        Generate XML element for frontend consumption.

        This is the primary method for serializing report elements
        to XML for the React frontend. Subclasses should override
        this method and call super().as_data_etree() first.

        Returns:
            ET.Element with tag 'CCP4i2Report{ClassName}'
        """
        # Tag name is CCP4i2Report + class name (e.g., CCP4i2ReportFold)
        root = ET.Element(
            f'CCP4i2Report{type(self).__name__}',
            key=self.internalId
        )

        # Standard attributes
        root.set('class', self.class_ if self.class_ is not None else '')
        root.set('style', self.style if self.style is not None else '')
        root.set('id', self.id if self.id is not None else '')

        if self.label is not None:
            root.set('label', self.label)
        if self.title is not None:
            root.set('title', self.title)

        # Text content
        root.text = self.text if (
            hasattr(self, 'text') and self.text is not None) else ''
        root.tail = self.tail if (
            hasattr(self, 'tail') and self.tail is not None) else ''

        return root

    def collect_all_diagnostics(self) -> DiagnosticCollector:
        """
        Collect diagnostics from this element and all descendants.

        Returns:
            DiagnosticCollector with all diagnostics
        """
        collector = DiagnosticCollector()
        collector.merge(self._diagnostics)

        # Collect from children if this is a container
        if hasattr(self, 'children'):
            for child in self.children:
                if hasattr(child, 'collect_all_diagnostics'):
                    child_diags = child.collect_all_diagnostics()
                    if child_diags is not None:
                        collector.merge(child_diags)
                elif hasattr(child, '_diagnostics') and child._diagnostics is not None:
                    collector.merge(child._diagnostics)

        return collector


class ReportContainer(ReportElement):
    """
    Base class for report elements that contain children.

    This provides the container functionality needed by Report, Fold,
    Results, Div, and other container elements. It maintains a list
    of child elements and provides methods for adding and iterating.

    The class is designed to be backward-compatible with the existing
    Container class in CCP4ReportParser.py.

    Attributes:
        children: List of child elements

    Example:
        container = ReportContainer(id="my_container")
        container.append(TextElement(text="Hello"))
        container.append(TableElement(data=my_data))

        for child in container:
            print(child)
    """

    def __init__(
        self,
        xrtnode: Optional[ET.Element] = None,
        xmlnode: Optional[ET.Element] = None,
        jobInfo: Optional[Dict[str, Any]] = None,
        **kwargs
    ):
        """
        Initialize a container element.

        Args:
            xrtnode: XRT template node (for legacy compatibility)
            xmlnode: XML data node (for legacy compatibility)
            jobInfo: Job information dictionary
            **kwargs: Additional arguments passed to ReportElement
        """
        super().__init__(**kwargs)

        self.children: List[ReportElement] = []
        self.xrtnode = xrtnode
        self.xmlnode = xmlnode
        self.jobInfo = jobInfo or {}

    def append(self, child: 'ReportElement'):
        """
        Add a child element to this container.

        Args:
            child: Element to add (can be ReportElement or string)
        """
        if isinstance(child, str):
            # Create a generic text element for strings
            from ccp4i2.report.CCP4ReportParser import Generic
            child = Generic(xmlnode=self.xmlnode, text=child)

        self.children.append(child)

        # Set parent relationship if child supports it
        if hasattr(child, 'parent'):
            child.parent = self

    def insert(self, index: int, child: 'ReportElement'):
        """Insert a child element at the specified index."""
        if isinstance(child, str):
            from ccp4i2.report.CCP4ReportParser import Generic
            child = Generic(xmlnode=self.xmlnode, text=child)

        self.children.insert(index, child)

        if hasattr(child, 'parent'):
            child.parent = self

    def __len__(self) -> int:
        """Return the number of children."""
        return len(self.children)

    def __iter__(self) -> Iterator[ReportElement]:
        """Iterate over children."""
        return iter(self.children)

    def __getitem__(self, index: int) -> ReportElement:
        """Get child by index."""
        return self.children[index]

    def as_data_etree(self) -> ET.Element:
        """
        Generate XML element including all children.

        Returns:
            ET.Element with children appended
        """
        root = super().as_data_etree()

        for child in self.children:
            if hasattr(child, 'as_data_etree'):
                try:
                    child_elem = child.as_data_etree()
                    if child_elem is not None:
                        root.append(child_elem)
                except Exception as e:
                    self.error(
                        "CHILD_RENDER_ERROR",
                        f"Failed to render child {type(child).__name__}: {e}",
                        exception=e
                    )
            else:
                logger.warning(f"Child {type(child).__name__} has no as_data_etree method")

        # Append diagnostics if any
        all_diagnostics = self.collect_all_diagnostics()
        if all_diagnostics.has_errors or all_diagnostics.has_warnings:
            root.append(all_diagnostics.as_data_etree())

        return root


class ElementRegistry:
    """
    Registry for report element types.

    This enables plugin extensions and provides a factory method
    for creating elements from XRT tags.

    Usage:
        # Register an element type
        ElementRegistry.register('mytag', MyElementClass)

        # Create an element from an XRT tag
        element = ElementRegistry.create('mytag', xrtnode, xmlnode, jobInfo)
    """

    _elements: Dict[str, Type[ReportElement]] = {}
    _factories: Dict[str, Callable] = {}

    @classmethod
    def register(cls, tag: str, element_class: Type[ReportElement]):
        """
        Register an element class for an XRT tag.

        Args:
            tag: XRT tag name (without namespace)
            element_class: Class to instantiate for this tag
        """
        cls._elements[tag] = element_class
        logger.debug(f"Registered element: {tag} -> {element_class.__name__}")

    @classmethod
    def register_factory(cls, tag: str, factory: Callable):
        """
        Register a factory function for creating elements.

        The factory should have signature:
            factory(xrtnode, xmlnode, jobInfo) -> ReportElement

        Args:
            tag: XRT tag name (without namespace)
            factory: Factory function
        """
        cls._factories[tag] = factory

    @classmethod
    def create(
        cls,
        tag: str,
        xrtnode: Optional[ET.Element] = None,
        xmlnode: Optional[ET.Element] = None,
        jobInfo: Optional[Dict[str, Any]] = None
    ) -> Optional[ReportElement]:
        """
        Create an element from an XRT tag.

        Args:
            tag: XRT tag name (without namespace)
            xrtnode: XRT template node
            xmlnode: XML data node
            jobInfo: Job information dictionary

        Returns:
            Created element, or None if tag not registered
        """
        # Check for factory first
        if tag in cls._factories:
            return cls._factories[tag](xrtnode, xmlnode, jobInfo)

        # Then check for element class
        if tag in cls._elements:
            return cls._elements[tag](
                xrtnode=xrtnode,
                xmlnode=xmlnode,
                jobInfo=jobInfo
            )

        return None

    @classmethod
    def is_registered(cls, tag: str) -> bool:
        """Check if a tag is registered."""
        return tag in cls._elements or tag in cls._factories

    @classmethod
    def list_tags(cls) -> List[str]:
        """Get list of all registered tags."""
        return list(set(cls._elements.keys()) | set(cls._factories.keys()))
