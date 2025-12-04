"""
Custom exceptions for CCP4i2 operations.

Provides well-defined exception types that can be caught and handled
consistently across the application.
"""

from typing import Dict, Any, Optional


class CCP4OperationError(Exception):
    """
    Base exception for CCP4 operations.

    All CCP4-specific exceptions inherit from this base class,
    allowing for catch-all error handling when needed.

    Attributes:
        message: Human-readable error message
        details: Optional dictionary with additional context
    """

    def __init__(self, message: str, details: Optional[Dict[str, Any]] = None):
        self.message = message
        self.details = details or {}
        super().__init__(message)

    def to_dict(self) -> Dict[str, Any]:
        """Convert exception to dictionary for API responses."""
        result = {
            "status": "Failed",
            "reason": self.message
        }
        if self.details:
            result["details"] = self.details
        return result


class JobNotFoundError(CCP4OperationError):
    """Job not found in database."""

    def __init__(self, job_id: str, details: Optional[Dict[str, Any]] = None):
        message = f"Job not found: {job_id}"
        super().__init__(message, details)
        self.job_id = job_id


class ProjectNotFoundError(CCP4OperationError):
    """Project not found in database."""

    def __init__(self, project_id: str, details: Optional[Dict[str, Any]] = None):
        message = f"Project not found: {project_id}"
        super().__init__(message, details)
        self.project_id = project_id


class ValidationError(CCP4OperationError):
    """Parameter or data validation failed."""

    def __init__(self, message: str, validation_errors: Optional[Dict[str, Any]] = None):
        details = {"validation_errors": validation_errors} if validation_errors else {}
        super().__init__(message, details)
        self.validation_errors = validation_errors


class FileOperationError(CCP4OperationError):
    """File operation (read, write, digest) failed."""

    def __init__(self, message: str, file_path: Optional[str] = None,
                 details: Optional[Dict[str, Any]] = None):
        if details is None:
            details = {}
        if file_path:
            details["file_path"] = file_path
        super().__init__(message, details)
        self.file_path = file_path


class ParameterError(CCP4OperationError):
    """Parameter setting or getting failed."""

    def __init__(self, message: str, parameter_name: Optional[str] = None,
                 details: Optional[Dict[str, Any]] = None):
        if details is None:
            details = {}
        if parameter_name:
            details["parameter_name"] = parameter_name
        super().__init__(message, details)
        self.parameter_name = parameter_name


class ExecutionError(CCP4OperationError):
    """Job execution failed."""

    def __init__(self, message: str, job_id: Optional[str] = None,
                 exit_code: Optional[int] = None,
                 details: Optional[Dict[str, Any]] = None):
        if details is None:
            details = {}
        if job_id:
            details["job_id"] = job_id
        if exit_code is not None:
            details["exit_code"] = exit_code
        super().__init__(message, details)
        self.job_id = job_id
        self.exit_code = exit_code
