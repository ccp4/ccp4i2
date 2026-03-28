# Copyright (C) 2025 University of York
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""CCP4 Error Handling module - compatibility wrapper.

This module provides backward compatibility by re-exporting the error handling
system from ccp4i2.core.base_object.error_reporting.

For new code, prefer importing directly from ccp4i2.core.base_object.error_reporting.
This module exists to support legacy CCP4i2 plugins that expect:
    from ccp4i2.core import CCP4ErrorHandling
"""

# Export all public API from error_reporting
from ccp4i2.core.base_object.error_reporting import (
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
