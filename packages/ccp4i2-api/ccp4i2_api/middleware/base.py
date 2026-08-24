"""Common base class for CCP4i2 auth middleware.

Subclasses implement ``is_active()`` (whether this middleware is configured
to operate in the current process) and ``authenticate(request)`` (the
auth-specific token-validation logic). The base class handles the rest of
the request lifecycle: the canonical 401 response shape, setting
``request.user``, and setting the trust flag that the matching DRF
authentication class checks before honouring ``request.user`` (anti-
spoofing, mirroring the existing AzureADAuthMiddleware contract).
"""

from django.contrib.auth import get_user_model
from django.http import HttpRequest, HttpResponse, JsonResponse

from ..exceptions import AuthenticationFailed, AuthorizationFailed
from ..file_grants import extract_grant, grant_user_pk


# Attribute name set on ``request`` after a successful authentication.
# The DRF authentication class trusts ``request.user`` only when this is
# set, preventing spoofing via direct attribute writes from other code.
REQUEST_FLAG_ATTR = "_ccp4i2_api_middleware_ran"


# Public endpoints that must be reachable WITHOUT authentication. These are
# unauthenticated by design (plain Django views, no DRF permissions):
#   * health  — liveness/readiness probe. Polled by container orchestrators and
#               by the desktop launch gate before any token exists.
#   * version — server version, fetched before login.
# Matched against the request path suffix so it works under any URL prefix
# (``/api/ccp4i2/health/`` in-app, ``/health/`` if mounted bare).
PUBLIC_PATH_SUFFIXES = ("/health/", "/health", "/version/", "/version")


def _is_public_path(path: str) -> bool:
    return path.endswith(PUBLIC_PATH_SUFFIXES)


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
        # Public probes (health/version) are unauthenticated by design and must
        # be reachable without a token — container liveness/readiness checks and
        # the desktop launch gate poll them before any credential exists.
        if _is_public_path(request.path):
            return self.get_response(request)
        try:
            user = self.authenticate(request)
        except AuthenticationFailed as exc:
            # Only once normal authentication has failed: a file grant, which
            # authenticates browser-issued reads that cannot carry an
            # Authorization header (a report page's own images, stylesheets
            # and relative fetches). Trying it second rather than first means
            # a request that does carry credentials is always the identity
            # they name — a grant cookie left over from a report someone
            # else opened can neither shadow a valid token nor borrow it.
            grant_user = self._grant_user(request)
            if grant_user is not None:
                request.user = grant_user
                setattr(request, REQUEST_FLAG_ATTR, True)
                return self.get_response(request)
            return self._error_response(str(exc), status=401)
        except AuthorizationFailed as exc:
            # A valid identity that is not permitted here stays refused: a
            # grant is a read capability, not a way around authorization.
            return self._error_response(str(exc), status=403, prefix="Access denied")
        request.user = user
        setattr(request, REQUEST_FLAG_ATTR, True)
        return self.get_response(request)

    @staticmethod
    def _grant_user(request: HttpRequest):
        """Return the user a valid file grant was minted by, else None."""
        token = extract_grant(request)
        if not token:
            return None
        user_pk = grant_user_pk(token, path=request.path, method=request.method)
        if user_pk is None:
            return None
        return get_user_model().objects.filter(pk=user_pk).first()

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
