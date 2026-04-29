"""Exceptions used by the shared auth middleware contract."""


class AuthenticationFailed(Exception):
    """Raised by ``BaseAuthMiddleware.authenticate()`` to signal a 401.

    The exception's first argument is the message returned in the canonical
    401 response body (``{"success": false, "error": "..."}``).
    """


class AuthorizationFailed(Exception):
    """Raised by ``BaseAuthMiddleware.authenticate()`` to signal a 403.

    Use when the request is *authenticated* (we know who is calling) but
    *not authorized* (e.g., not a member of an allowed group). The first
    argument is the message returned in the canonical 403 response body.
    """
