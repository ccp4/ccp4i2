"""Exceptions used by the shared auth middleware contract."""


class AuthenticationFailed(Exception):
    """Raised by ``BaseAuthMiddleware.authenticate()`` to signal a 401.

    The exception's first argument is the message returned in the canonical
    401 response body (``{"success": false, "error": "..."}``).
    """
