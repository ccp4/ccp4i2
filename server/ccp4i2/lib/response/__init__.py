"""
Standardized response types for CCP4i2 operations.

This module provides consistent result and error handling across
the library, API, and CLI layers.

API Response Format:
    Success: {"success": true, "data": {...}}
    Failure: {"success": false, "error": "message", "details": {...}}

Usage in ViewSets:
    from ccp4i2.lib.response import api_success, api_error, Result

    # Quick responses
    return api_success({"job_id": 123})
    return api_error("Job not found", status=404)

    # From Result objects (returned by library functions)
    result = some_operation()
    return result.to_response()
"""

from .result import Result, OperationResult, api_success, api_error
from .exceptions import (
    CCP4OperationError,
    JobNotFoundError,
    ProjectNotFoundError,
    ValidationError,
    FileOperationError,
    ParameterError,
)

__all__ = [
    'Result',
    'OperationResult',
    'api_success',
    'api_error',
    'CCP4OperationError',
    'JobNotFoundError',
    'ProjectNotFoundError',
    'ValidationError',
    'FileOperationError',
    'ParameterError',
]
