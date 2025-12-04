"""Error reporting and validation system for CData objects.

This module implements the CCP4i2-compatible error handling system, providing
error tracking, severity levels, and validation reporting capabilities.
"""

from typing import List, Optional, Dict, Any
from enum import IntEnum
import xml.etree.ElementTree as ET


class Severity(IntEnum):
    """Error severity levels for validation and error reporting."""
    OK = 0
    UNDEFINED = 1
    WARNING = 2
    UNDEFINED_ERROR = 3
    ERROR = 4


# Legacy compatibility constants
SEVERITY_OK = Severity.OK
SEVERITY_UNDEFINED = Severity.UNDEFINED
SEVERITY_WARNING = Severity.WARNING
SEVERITY_UNDEFINED_ERROR = Severity.UNDEFINED_ERROR
SEVERITY_ERROR = Severity.ERROR


class CErrorReport:
    """Error report container for tracking validation errors and warnings.

    CErrorReport holds a list of errors, each containing:
    - klass: The class name where the error occurred
    - code: Error code number (from ERROR_CODES dictionary)
    - details: Additional details about the error
    - name: Object name where error occurred
    - severity: Severity level (0=OK, 1=UNDEFINED, 2=WARNING, 3=UNDEFINED_ERROR, 4=ERROR)

    Example:
        >>> report = CErrorReport()
        >>> report.append("CInt", 101, "Value exceeds maximum", "myInt", Severity.ERROR)
        >>> print(report.maxSeverity())
        4
        >>> print(report.report())
        ERROR in CInt 'myInt': Value exceeds maximum (code 101)
    """

    def __init__(self):
        """Initialize an empty error report."""
        self._errors: List[Dict[str, Any]] = []

    def append(self, klass: str, code: int, details: str, name: str = "",
               severity: int = SEVERITY_ERROR):
        """Add a single error to the report.

        Args:
            klass: Class name where error occurred
            code: Error code number
            details: Error description/details
            name: Object name (optional)
            severity: Severity level (default: SEVERITY_ERROR)
        """
        self._errors.append({
            'class': klass,
            'code': code,
            'details': details,
            'name': name,
            'severity': severity
        })

    def extend(self, other: 'CErrorReport'):
        """Merge another error report into this one.

        Args:
            other: Another CErrorReport to merge
        """
        if isinstance(other, CErrorReport):
            self._errors.extend(other._errors)

    def downgrade_to_warnings(self):
        """Downgrade all errors with severity > WARNING to WARNING.

        This is used when child validation errors should not block execution
        because the parent object is optional and not set.
        """
        for error in self._errors:
            if error['severity'] > SEVERITY_WARNING:
                error['severity'] = SEVERITY_WARNING

    def count(self, cls=None, code=None) -> int:
        """Return the number of errors in this report.

        Args:
            cls: Optional class filter (legacy Qt API compatibility)
            code: Optional error code filter

        Returns:
            Number of errors matching the filters (or all errors if no filters)
        """
        if cls is None and code is None:
            return len(self._errors)

        # Filter errors by cls and/or code
        count = 0
        for error in self._errors:
            if code is not None and error.get('code') != code:
                continue
            # cls parameter is ignored in Qt-free implementation
            # (it was used to filter by plugin class in Qt version)
            count += 1
        return count

    def maxSeverity(self) -> int:
        """Return the maximum severity level of all errors.

        Returns:
            Maximum severity (0-4), or SEVERITY_OK if no errors
        """
        if not self._errors:
            return SEVERITY_OK
        return max(error['severity'] for error in self._errors)

    def report(self, severity_threshold: int = SEVERITY_WARNING) -> str:
        """Generate a formatted report of errors at or above severity threshold.

        Args:
            severity_threshold: Minimum severity to include (default: WARNING)

        Returns:
            Formatted error report string
        """
        if not self._errors:
            return "No errors"

        lines = []
        severity_names = {
            SEVERITY_OK: "OK",
            SEVERITY_UNDEFINED: "UNDEFINED",
            SEVERITY_WARNING: "WARNING",
            SEVERITY_UNDEFINED_ERROR: "UNDEFINED_ERROR",
            SEVERITY_ERROR: "ERROR"
        }

        for error in self._errors:
            if error['severity'] >= severity_threshold:
                severity_name = severity_names.get(error['severity'], "UNKNOWN")
                name_part = f" '{error['name']}'" if error['name'] else ""
                lines.append(
                    f"{severity_name} in {error['class']}{name_part}: "
                    f"{error['details']} (code {error['code']})"
                )

        return "\n".join(lines) if lines else f"No errors at or above severity {severity_threshold}"

    def __str__(self) -> str:
        """Return a string representation of the error report.

        Returns:
            Formatted error report
        """
        return self.report(SEVERITY_WARNING)

    def __bool__(self) -> bool:
        """Return True if there are errors with severity >= WARNING.

        Returns:
            True if there are warnings or errors
        """
        return self.maxSeverity() >= SEVERITY_WARNING

    def __len__(self) -> int:
        """Return the number of errors.

        Returns:
            Number of errors
        """
        return self.count()

    def getErrors(self) -> List[Dict[str, Any]]:
        """Get the list of all errors.

        Returns:
            List of error dictionaries
        """
        return self._errors.copy()

    def clear(self):
        """Clear all errors from the report."""
        self._errors.clear()

    def getEtree(self) -> ET.Element:
        """Serialize error report to an XML ElementTree element.

        Creates an XML structure compatible with CCP4i2 diagnostic.xml format.

        Returns:
            ET.Element: Root 'errorReportList' element containing 'errorReport' children

        Example output:
            <errorReportList>
                <errorReport>
                    <class>CInt</class>
                    <code>101</code>
                    <details>Value out of range</details>
                    <name>myInt</name>
                    <severity>4</severity>
                </errorReport>
            </errorReportList>
        """
        severity_names = {
            SEVERITY_OK: "OK",
            SEVERITY_UNDEFINED: "UNDEFINED",
            SEVERITY_WARNING: "WARNING",
            SEVERITY_UNDEFINED_ERROR: "UNDEFINED_ERROR",
            SEVERITY_ERROR: "ERROR"
        }

        root = ET.Element("errorReportList")

        for error in self._errors:
            error_elem = ET.SubElement(root, "errorReport")

            # Add class element
            class_elem = ET.SubElement(error_elem, "class")
            class_elem.text = str(error.get('class', ''))

            # Add code element
            code_elem = ET.SubElement(error_elem, "code")
            code_elem.text = str(error.get('code', 0))

            # Add details element
            details_elem = ET.SubElement(error_elem, "details")
            details_elem.text = str(error.get('details', ''))

            # Add name element
            name_elem = ET.SubElement(error_elem, "name")
            name_elem.text = str(error.get('name', ''))

            # Add severity element (both as number and name)
            severity = error.get('severity', SEVERITY_OK)
            severity_elem = ET.SubElement(error_elem, "severity")
            severity_elem.text = str(severity)

            severity_name_elem = ET.SubElement(error_elem, "severityName")
            severity_name_elem.text = severity_names.get(severity, "UNKNOWN")

        return root

    # ========================================================================
    # Legacy CCP4i2 Compatibility
    # ========================================================================

    @property
    def _reports(self):
        """Legacy compatibility: Alias for _errors.

        Old CCP4i2 code expects _reports instead of _errors.
        """
        return self._errors

    def get_broken_files(self) -> List[str]:
        """Get list of file names that have errors (legacy CCP4i2 API).

        This method provides backward compatibility with legacy plugins that
        expect checkInputData() to return a list of broken file names.

        Returns:
            List of file names that have validation errors
        """
        broken_files = []
        for error in self._errors:
            if error.get('name') and error.get('name') not in broken_files:
                broken_files.append(error['name'])
        return broken_files

    def __contains__(self, item: str) -> bool:
        """Check if a file name has errors (legacy compatibility).

        Allows legacy code to use: if 'XYZIN' in error_report

        Args:
            item: File name to check

        Returns:
            True if the file name has validation errors
        """
        return item in self.get_broken_files()

    def __iter__(self):
        """Iterate over broken file names (legacy compatibility).

        Allows legacy code to use: for filename in error_report

        Returns:
            Iterator over broken file names
        """
        return iter(self.get_broken_files())

    def remove(self, item: str):
        """Remove all errors for a specific file name (legacy compatibility).

        Allows legacy code to use: error_report.remove('XYZIN')

        Args:
            item: File name to remove errors for
        """
        self._errors = [e for e in self._errors if e.get('name') != item]


class CException(CErrorReport, Exception):
    """Exception class that combines CErrorReport functionality with Python exceptions.

    This allows errors to be both raised as exceptions and tracked in error reports.

    Example:
        >>> exc = CException("CInt", 101, "Value out of range", "myInt", Severity.ERROR)
        >>> raise exc
    """

    def __init__(self, klass: str = "", code: int = 0, details: str = "",
                 name: str = "", severity: int = SEVERITY_ERROR):
        """Initialize exception with error details.

        Args:
            klass: Class name where error occurred
            code: Error code number
            details: Error description
            name: Object name
            severity: Severity level
        """
        CErrorReport.__init__(self)
        Exception.__init__(self, details)

        if klass or code or details:
            self.append(klass, code, details, name, severity)
