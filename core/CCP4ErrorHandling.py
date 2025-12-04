"""CCP4 Error Handling module - compatibility wrapper.

This module provides backward compatibility by re-exporting the error handling
system from core.base_object.error_reporting.

For new code, prefer importing directly from core.base_object.error_reporting.
This module exists to support legacy CCP4i2 plugins that expect:
    from core import CCP4ErrorHandling
"""

# Export all public API from error_reporting
from core.base_object.error_reporting import (
    # Severity enum and constants
    Severity,
    SEVERITY_OK,
    SEVERITY_UNDEFINED,
    SEVERITY_WARNING,
    SEVERITY_UNDEFINED_ERROR,
    SEVERITY_ERROR,

    # Main classes
    CErrorReport,
    CException,
)

__all__ = [
    'Severity',
    'SEVERITY_OK',
    'SEVERITY_UNDEFINED',
    'SEVERITY_WARNING',
    'SEVERITY_UNDEFINED_ERROR',
    'SEVERITY_ERROR',
    'CErrorReport',
    'CException',
]
