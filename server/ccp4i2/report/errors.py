"""
Structured error handling and diagnostics for CCP4i2 reports.

This module provides a modern, structured approach to collecting and
reporting errors and diagnostics during report generation. Diagnostics
flow through to the frontend for clean user-facing error messages.

Usage:
    from ccp4i2.report.errors import DiagnosticCollector, ReportDiagnostic

    collector = DiagnosticCollector()
    collector.warning("TABLE_EMPTY", "Table has no data rows")
    collector.error("FILE_NOT_FOUND", "Could not load data file",
                    details={"path": "/path/to/file"})

    # Get XML for frontend
    for diagnostic in collector:
        report.append_diagnostic(diagnostic)
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List, Iterator
from enum import Enum
from xml.etree import ElementTree as ET
import logging
import traceback

logger = logging.getLogger(f"ccp4x:{__name__}")


class DiagnosticLevel(Enum):
    """Severity levels for diagnostics."""
    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


@dataclass
class ReportDiagnostic:
    """
    Structured diagnostic information for frontend display.

    Diagnostics are collected during report generation and can be
    serialized to XML for transmission to the frontend, where they
    can be displayed in a consistent, user-friendly manner.

    Attributes:
        level: Severity level (debug, info, warning, error, critical)
        code: Machine-readable error code (e.g., "FILE_NOT_FOUND")
        message: Human-readable description
        location: Where the error occurred (e.g., "Table.as_data_etree")
        details: Additional context as key-value pairs
        exception: Optional exception that caused this diagnostic
    """
    level: DiagnosticLevel
    code: str
    message: str
    location: Optional[str] = None
    details: Optional[Dict[str, Any]] = None
    exception: Optional[Exception] = None

    def __post_init__(self):
        # Convert string level to enum if needed
        if isinstance(self.level, str):
            self.level = DiagnosticLevel(self.level.lower())

    def as_data_etree(self) -> ET.Element:
        """Generate XML element for frontend consumption."""
        elem = ET.Element('CCP4i2ReportDiagnostic')
        elem.set('level', self.level.value)
        elem.set('code', self.code)
        elem.set('message', self.message)

        if self.location:
            elem.set('location', self.location)

        if self.details:
            details_elem = ET.SubElement(elem, 'details')
            for key, value in self.details.items():
                detail = ET.SubElement(details_elem, 'detail')
                detail.set('key', str(key))
                detail.text = str(value)

        if self.exception:
            exc_elem = ET.SubElement(elem, 'exception')
            exc_elem.set('type', type(self.exception).__name__)
            exc_elem.text = str(self.exception)

        return elem

    def log(self):
        """Log this diagnostic using the standard logger."""
        log_message = f"[{self.code}] {self.message}"
        if self.location:
            log_message = f"{self.location}: {log_message}"
        if self.details:
            log_message += f" {self.details}"

        level_map = {
            DiagnosticLevel.DEBUG: logger.debug,
            DiagnosticLevel.INFO: logger.info,
            DiagnosticLevel.WARNING: logger.warning,
            DiagnosticLevel.ERROR: logger.error,
            DiagnosticLevel.CRITICAL: logger.critical,
        }
        log_func = level_map.get(self.level, logger.error)
        log_func(log_message, exc_info=self.exception if self.level in (
            DiagnosticLevel.ERROR, DiagnosticLevel.CRITICAL) else None)


class DiagnosticCollector:
    """
    Collects diagnostics during report generation.

    This collector gathers all diagnostics (errors, warnings, info messages)
    during report generation for later serialization and display.

    Usage:
        collector = DiagnosticCollector()

        try:
            # some operation
        except FileNotFoundError as e:
            collector.error("FILE_NOT_FOUND", f"Could not load {path}",
                          exception=e, details={"path": path})

        # Check if there were errors
        if collector.has_errors:
            # Handle error case
            pass

        # Iterate over all diagnostics
        for diagnostic in collector:
            print(diagnostic.message)
    """

    def __init__(self):
        self._diagnostics: List[ReportDiagnostic] = []

    def add(
        self,
        level: DiagnosticLevel,
        code: str,
        message: str,
        location: Optional[str] = None,
        details: Optional[Dict[str, Any]] = None,
        exception: Optional[Exception] = None,
        log: bool = True
    ) -> ReportDiagnostic:
        """
        Add a diagnostic to the collection.

        Args:
            level: Severity level
            code: Machine-readable error code
            message: Human-readable description
            location: Where the error occurred
            details: Additional context
            exception: Exception that caused this diagnostic
            log: Whether to also log this diagnostic (default True)

        Returns:
            The created ReportDiagnostic
        """
        diagnostic = ReportDiagnostic(
            level=level,
            code=code,
            message=message,
            location=location,
            details=details,
            exception=exception
        )
        self._diagnostics.append(diagnostic)

        if log:
            diagnostic.log()

        return diagnostic

    def debug(self, code: str, message: str, **kwargs) -> ReportDiagnostic:
        """Add a debug-level diagnostic."""
        return self.add(DiagnosticLevel.DEBUG, code, message, **kwargs)

    def info(self, code: str, message: str, **kwargs) -> ReportDiagnostic:
        """Add an info-level diagnostic."""
        return self.add(DiagnosticLevel.INFO, code, message, **kwargs)

    def warning(self, code: str, message: str, **kwargs) -> ReportDiagnostic:
        """Add a warning-level diagnostic."""
        return self.add(DiagnosticLevel.WARNING, code, message, **kwargs)

    def error(self, code: str, message: str, **kwargs) -> ReportDiagnostic:
        """Add an error-level diagnostic."""
        return self.add(DiagnosticLevel.ERROR, code, message, **kwargs)

    def critical(self, code: str, message: str, **kwargs) -> ReportDiagnostic:
        """Add a critical-level diagnostic."""
        return self.add(DiagnosticLevel.CRITICAL, code, message, **kwargs)

    def from_exception(
        self,
        exception: Exception,
        code: Optional[str] = None,
        location: Optional[str] = None
    ) -> ReportDiagnostic:
        """
        Create a diagnostic from an exception.

        Args:
            exception: The exception to convert
            code: Optional error code (defaults to exception class name)
            location: Where the exception occurred

        Returns:
            The created ReportDiagnostic
        """
        if code is None:
            code = type(exception).__name__.upper()

        return self.error(
            code=code,
            message=str(exception),
            location=location,
            exception=exception,
            details={"traceback": traceback.format_exc()}
        )

    def __iter__(self) -> Iterator[ReportDiagnostic]:
        """Iterate over all diagnostics."""
        return iter(self._diagnostics)

    def __len__(self) -> int:
        """Return the number of diagnostics."""
        return len(self._diagnostics)

    def __bool__(self) -> bool:
        """Return True if there are any diagnostics."""
        return len(self._diagnostics) > 0

    @property
    def has_errors(self) -> bool:
        """Check if there are any error or critical level diagnostics."""
        return any(
            d.level in (DiagnosticLevel.ERROR, DiagnosticLevel.CRITICAL)
            for d in self._diagnostics
        )

    @property
    def has_warnings(self) -> bool:
        """Check if there are any warning level or higher diagnostics."""
        return any(
            d.level in (DiagnosticLevel.WARNING, DiagnosticLevel.ERROR, DiagnosticLevel.CRITICAL)
            for d in self._diagnostics
        )

    def filter_by_level(self, *levels: DiagnosticLevel) -> List[ReportDiagnostic]:
        """Get diagnostics filtered by level(s)."""
        return [d for d in self._diagnostics if d.level in levels]

    @property
    def errors(self) -> List[ReportDiagnostic]:
        """Get all error and critical level diagnostics."""
        return self.filter_by_level(DiagnosticLevel.ERROR, DiagnosticLevel.CRITICAL)

    @property
    def warnings(self) -> List[ReportDiagnostic]:
        """Get all warning level diagnostics."""
        return self.filter_by_level(DiagnosticLevel.WARNING)

    def clear(self):
        """Clear all diagnostics."""
        self._diagnostics.clear()

    def as_data_etree(self) -> ET.Element:
        """
        Generate XML element containing all diagnostics.

        Returns:
            Element with all diagnostics as children
        """
        root = ET.Element('CCP4i2ReportDiagnostics')
        root.set('count', str(len(self._diagnostics)))
        root.set('hasErrors', str(self.has_errors).lower())
        root.set('hasWarnings', str(self.has_warnings).lower())

        for diagnostic in self._diagnostics:
            root.append(diagnostic.as_data_etree())

        return root

    def merge(self, other: 'DiagnosticCollector'):
        """Merge diagnostics from another collector."""
        if other is not None and hasattr(other, '_diagnostics') and other._diagnostics is not None:
            self._diagnostics.extend(other._diagnostics)


# Standard error codes for report generation
class ReportErrorCodes:
    """Standard error codes for report generation."""

    # File-related errors
    FILE_NOT_FOUND = "FILE_NOT_FOUND"
    FILE_READ_ERROR = "FILE_READ_ERROR"
    FILE_PARSE_ERROR = "FILE_PARSE_ERROR"

    # Data-related errors
    MISSING_DATA = "MISSING_DATA"
    INVALID_DATA = "INVALID_DATA"
    EMPTY_DATA = "EMPTY_DATA"

    # Report class errors
    REPORT_CLASS_NOT_FOUND = "REPORT_CLASS_NOT_FOUND"
    REPORT_INIT_ERROR = "REPORT_INIT_ERROR"
    REPORT_RENDER_ERROR = "REPORT_RENDER_ERROR"

    # XML-related errors
    XML_PARSE_ERROR = "XML_PARSE_ERROR"
    XML_GENERATION_ERROR = "XML_GENERATION_ERROR"
    XRT_NOT_FOUND = "XRT_NOT_FOUND"

    # Job-related errors
    JOB_NOT_FOUND = "JOB_NOT_FOUND"
    JOB_INFO_ERROR = "JOB_INFO_ERROR"

    # Element-specific warnings
    TABLE_EMPTY = "TABLE_EMPTY"
    GRAPH_NO_DATA = "GRAPH_NO_DATA"
    PICTURE_NOT_FOUND = "PICTURE_NOT_FOUND"
