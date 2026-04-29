"""Common base class for CCP4i2 auth middleware.

Subclasses implement ``is_active()`` (whether this middleware is configured
to operate in the current process) and ``authenticate(request)`` (the
auth-specific token-validation logic). The base class handles the rest of
the request lifecycle: the canonical 401 response shape, setting
``request.user``, and setting the trust flag that the matching DRF
authentication class checks before honouring ``request.user`` (anti-
spoofing, mirroring the existing AzureADAuthMiddleware contract).
"""

from django.http import HttpRequest, HttpResponse, JsonResponse

from ..exceptions import AuthenticationFailed, AuthorizationFailed


# Attribute name set on ``request`` after a successful authentication.
# The DRF authentication class trusts ``request.user`` only when this is
# set, preventing spoofing via direct attribute writes from other code.
REQUEST_FLAG_ATTR = "_ccp4i2_auth_middleware_ran"


class BaseAuthMiddleware:
    """Abstract base for CCP4i2 auth middleware.

    Subclasses must implement:

    * ``is_active(self) -> bool`` — return True iff this middleware should
      attempt to authenticate. Typical implementation checks an env var
      or Django setting.
    * ``authenticate(self, request) -> User`` — return a Django User on
      successful auth. Raise ``AuthenticationFailed`` to signal a 401.

    When ``is_active()`` is False, the middleware is a no-op (the request
    flows to the next middleware unchanged). This lets multiple subclasses
    coexist in ``MIDDLEWARE`` without coupling between them; deployments
    activate the right one via configuration (env var presence).
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request: HttpRequest) -> HttpResponse:
        if not self.is_active():
            return self.get_response(request)
        try:
            user = self.authenticate(request)
        except AuthenticationFailed as exc:
            return self._error_response(str(exc), status=401)
        except AuthorizationFailed as exc:
            return self._error_response(str(exc), status=403, prefix="Access denied")
        request.user = user
        setattr(request, REQUEST_FLAG_ATTR, True)
        return self.get_response(request)

    def is_active(self) -> bool:
        raise NotImplementedError

    def authenticate(self, request: HttpRequest):
        raise NotImplementedError

    @staticmethod
    def _error_response(
        message: str, status: int, prefix: str = "Authentication failed"
    ) -> JsonResponse:
        return JsonResponse(
            {"success": False, "error": f"{prefix}: {message}"},
            status=status,
        )
