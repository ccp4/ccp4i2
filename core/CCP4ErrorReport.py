"""CCP4 Error Report - Re-export from base_object.error_reporting.

This module provides backward compatibility for legacy imports.
The actual implementation is in core.base_object.error_reporting.
"""

from core.base_object.error_reporting import (
    CErrorReport,
    CException,
    Severity,
    SEVERITY_OK,
    SEVERITY_UNDEFINED,
    SEVERITY_WARNING,
    SEVERITY_UNDEFINED_ERROR,
    SEVERITY_ERROR,
)

__all__ = [
    'CErrorReport',
    'CException',
    'Severity',
    'SEVERITY_OK',
    'SEVERITY_UNDEFINED',
    'SEVERITY_WARNING',
    'SEVERITY_UNDEFINED_ERROR',
    'SEVERITY_ERROR',
]
