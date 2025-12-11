"""
Result type for standardized operation returns.

Provides a consistent way to return success/failure from library functions
that can be easily consumed by both API endpoints and management commands.

API Response Format Convention:
    Success: {"success": true, "data": {...}}
    Failure: {"success": false, "error": "message", "details": {...}}

Usage in ViewSets:
    from ccp4i2.lib.response import api_success, api_error, Result

    # Quick responses
    return api_success({"job_id": 123})
    return api_error("Job not found", status=404)

    # From Result objects
    result = some_operation()
    return result.to_response()
"""

from typing import TypeVar, Generic, Optional, Dict, Any
from dataclasses import dataclass, field
from rest_framework.response import Response
from rest_framework import status as http_status

T = TypeVar('T')


@dataclass
class Result(Generic[T]):
    """
    Standardized result type for library operations.

    Encapsulates either a successful result with data, or a failure with
    error information. This provides a consistent interface for both
    API endpoints and CLI commands to handle operation outcomes.

    Attributes:
        success: Whether the operation succeeded
        data: Result data (if successful)
        error: Error message (if failed)
        error_details: Additional error context (if failed)

    Example:
        >>> result = Result.ok({"job_id": 123})
        >>> if result.success:
        ...     print(result.data)

        >>> result = Result.fail("Job not found", {"job_id": 999})
        >>> print(result.error)  # "Job not found"
    """

    success: bool
    data: Optional[T] = None
    error: Optional[str] = None
    error_details: Optional[Dict[str, Any]] = field(default_factory=dict)

    @classmethod
    def ok(cls, data: T) -> 'Result[T]':
        """
        Create a successful result.

        Args:
            data: The result data

        Returns:
            Result with success=True and the provided data
        """
        return cls(success=True, data=data)

    @classmethod
    def fail(cls, error: str, details: Optional[Dict[str, Any]] = None) -> 'Result[T]':
        """
        Create a failed result.

        Args:
            error: Error message
            details: Optional additional error context

        Returns:
            Result with success=False and error information
        """
        return cls(
            success=False,
            error=error,
            error_details=details or {}
        )

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary for API responses.

        Returns:
            Dictionary with success flag and data or error

        Example:
            Success: {"success": true, "data": {...}}
            Failure: {"success": false, "error": "...", "details": {...}}
        """
        if self.success:
            return {
                "success": True,
                "data": self.data
            }
        else:
            result = {
                "success": False,
                "error": self.error
            }
            if self.error_details:
                result["details"] = self.error_details
            return result

    def unwrap(self) -> T:
        """
        Unwrap the result, raising an exception if failed.

        Returns:
            The data if successful

        Raises:
            ValueError: If the result is a failure

        Example:
            >>> result = Result.ok(42)
            >>> value = result.unwrap()  # 42

            >>> result = Result.fail("Error")
            >>> value = result.unwrap()  # Raises ValueError
        """
        if self.success:
            return self.data
        else:
            raise ValueError(f"Cannot unwrap failed result: {self.error}")

    def unwrap_or(self, default: T) -> T:
        """
        Unwrap the result or return a default value.

        Args:
            default: Value to return if result is a failure

        Returns:
            The data if successful, otherwise the default
        """
        return self.data if self.success else default

    def map(self, func):
        """
        Map a function over the result data.

        If the result is successful, applies func to the data.
        If the result is a failure, returns the failure unchanged.

        Args:
            func: Function to apply to successful data

        Returns:
            New Result with transformed data or original failure
        """
        if self.success:
            return Result.ok(func(self.data))
        else:
            return self


    def to_response(self, success_status: int = 200, error_status: int = 400) -> Response:
        """
        Convert to a DRF Response object.

        Args:
            success_status: HTTP status code for successful responses (default 200)
            error_status: HTTP status code for error responses (default 400)

        Returns:
            DRF Response with appropriate status code

        Example:
            >>> result = Result.ok({"id": 123})
            >>> return result.to_response()  # Response with status 200

            >>> result = Result.fail("Not found")
            >>> return result.to_response(error_status=404)  # Response with status 404
        """
        if self.success:
            return Response(self.to_dict(), status=success_status)
        else:
            return Response(self.to_dict(), status=error_status)


# Convenience alias for clarity in some contexts
OperationResult = Result


# =============================================================================
# Standalone helper functions for quick API responses
# =============================================================================

def api_success(data: Any, status: int = 200) -> Response:
    """
    Create a successful API response.

    Args:
        data: Response data (will be wrapped in {"success": true, "data": ...})
        status: HTTP status code (default 200)

    Returns:
        DRF Response

    Example:
        return api_success({"job_id": 123})
        return api_success(serializer.data, status=201)
    """
    return Response({"success": True, "data": data}, status=status)


def api_error(
    error: str,
    status: int = 400,
    details: Optional[Dict[str, Any]] = None
) -> Response:
    """
    Create an error API response.

    Args:
        error: Error message
        status: HTTP status code (default 400)
        details: Optional additional error context

    Returns:
        DRF Response

    Example:
        return api_error("Job not found", status=404)
        return api_error("Validation failed", details={"field": "name"})
    """
    response_data = {"success": False, "error": error}
    if details:
        response_data["details"] = details
    return Response(response_data, status=status)
