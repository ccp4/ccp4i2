"""DRF authentication classes for CCP4i2.

These classes work in tandem with the middleware in
``ccp4i2_api.middleware`` — the middleware validates the bearer token
and sets ``request.user`` plus the ``REQUEST_FLAG_ATTR`` trust flag; the
DRF auth class checks the flag and surfaces ``request.user`` to DRF's
``IsAuthenticated`` permission.
"""

from rest_framework.authentication import BaseAuthentication

from .middleware.base import REQUEST_FLAG_ATTR


class AzureADAuthentication(BaseAuthentication):
    """
    DRF authentication class that uses the user set by AzureADAuthMiddleware.

    This allows DRF's IsAuthenticated permission to work with our middleware.
    The middleware does the actual JWT validation; this class just passes
    the authenticated user to DRF.

    In dev/Electron mode (CCP4I2_REQUIRE_AUTH not set), the middleware
    auto-assigns a dev_admin user, which this class also recognizes.

    Security: Only trusts users when our middleware has explicitly processed
    the request (marked by ``REQUEST_FLAG_ATTR`` attribute). This prevents
    spoofing attacks where a request might have ``request.user`` set by
    some other means.

    Note: this class is bound to the trust flag, not to the AzureAD chain
    specifically; it works equally for any middleware that inherits from
    ``BaseAuthMiddleware`` (e.g., LocalSessionAuthMiddleware in desktop
    mode), because they all set the same flag.
    """

    def authenticate(self, request):
        """
        Return the user if already authenticated by middleware, None otherwise.

        Returns:
            Tuple of (user, None) if authenticated, None if not.
        """
        # Check if middleware already validated and set user
        # The middleware sets these attributes on the underlying Django request
        django_request = getattr(request, '_request', request)

        # Security check: only trust users set by our middleware
        # This prevents spoofing where request.user might be set elsewhere
        if not getattr(django_request, REQUEST_FLAG_ATTR, False):
            return None

        # Middleware ran - trust the user it set
        user = getattr(django_request, 'user', None)
        if user and user.is_authenticated and not user.is_anonymous:
            return (user, None)

        return None
