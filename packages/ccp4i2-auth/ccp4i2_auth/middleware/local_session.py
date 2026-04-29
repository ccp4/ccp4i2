"""Per-launch token middleware for the CCP4i2 desktop app.

Validates ``Authorization: Bearer <token>`` against the secret in
``CCP4I2_LOCAL_SESSION_TOKEN``, set by the Electron main process when
spawning Django. If the env var is unset, the middleware is a no-op —
cloud deployments use ``AzureADAuthMiddleware`` instead.

Identity: the request is authenticated as the OS user who launched the
desktop app. The OS-user-derived email is computed by Electron and passed
via ``CCP4I2_LOCAL_USER_EMAIL`` (Electron sanitises the username so
domain-joined Windows boxes don't produce malformed emails). If that env
var is unset, the middleware falls back to a fixed default that is
well-formed across all platforms.
"""

import hmac
import os

from django.contrib.auth import get_user_model
from django.http import HttpRequest

from ..exceptions import AuthenticationFailed
from .base import BaseAuthMiddleware


# RFC 6761 reserves ``.invalid`` for guaranteed-non-resolvable identifiers,
# so this email is well-formed under Django's EmailValidator and cannot
# collide with a real account on any platform.
DEFAULT_DESKTOP_EMAIL = "desktop@ccp4i2.invalid"


class LocalSessionAuthMiddleware(BaseAuthMiddleware):

    def __init__(self, get_response):
        super().__init__(get_response)
        self.expected_token = os.environ.get("CCP4I2_LOCAL_SESSION_TOKEN")

    def is_active(self) -> bool:
        return self.expected_token is not None

    def authenticate(self, request: HttpRequest):
        auth = request.META.get("HTTP_AUTHORIZATION", "")
        if not auth.startswith("Bearer "):
            raise AuthenticationFailed("Missing Bearer token")
        provided = auth[len("Bearer "):]
        # Constant-time compare against length-extension / timing attacks.
        # Loopback-only is unlikely to be exploitable, but the precedent
        # is cheap to set for future cloud-side providers.
        if not hmac.compare_digest(provided, self.expected_token):
            raise AuthenticationFailed("Invalid local-session token")
        return self._desktop_user()

    @staticmethod
    def _desktop_user():
        User = get_user_model()
        email = os.environ.get("CCP4I2_LOCAL_USER_EMAIL", DEFAULT_DESKTOP_EMAIL)
        user, _ = User.objects.get_or_create(
            email=email,
            defaults={
                "username": email.split("@")[0],
                "is_staff": True,
                "is_superuser": True,
            },
        )
        return user
