"""Development-only auto-login middleware.

Auto-assigns a ``dev_admin`` Django superuser to every request — useful
for local Docker Compose development where you want CCP4i2 fully usable
without any token machinery. **Never enable this in production.**

Two defensive measures protect against accidental production exposure:

1. ``is_active()`` returns False unless ``settings.DEBUG`` is True. Even
   if this middleware is mistakenly listed in a production ``MIDDLEWARE``
   setting, it refuses to activate.
2. The CCP4i2 settings module only inserts this middleware as the
   *fallback* branch (after LocalSession + AzureAD env-var checks) and
   only when ``DEBUG`` is True. A production-shaped deploy with no auth
   env vars set falls through to *no auth middleware at all*, leaving
   requests as ``AnonymousUser`` — DRF's ``IsAuthenticated`` then 401s
   them. This is strictly safer than the previous "auto-create dev_admin
   on missing config" behaviour.

This middleware deliberately does *not* mirror the previous AzureAD
fallback's permissive default; it is opt-in via ``MIDDLEWARE`` and
DEBUG-gated.
"""

import logging

from django.conf import settings
from django.contrib.auth import get_user_model
from django.http import HttpRequest

from .base import BaseAuthMiddleware

logger = logging.getLogger(__name__)

DEV_ADMIN_USERNAME = "dev_admin"
DEV_ADMIN_EMAIL = "dev_admin@localhost"


class DevAdminMiddleware(BaseAuthMiddleware):
    """Auto-assigns a dev_admin superuser when DEBUG is True."""

    def is_active(self) -> bool:
        return getattr(settings, "DEBUG", False)

    def authenticate(self, request: HttpRequest):
        User = get_user_model()
        user, created = User.objects.get_or_create(
            username=DEV_ADMIN_USERNAME,
            defaults={
                "email": DEV_ADMIN_EMAIL,
                "first_name": "Dev",
                "last_name": "Admin",
                "is_staff": True,
                "is_superuser": True,
            },
        )
        if created:
            logger.info("Created dev_admin superuser for DEBUG mode")
        # If a previous run created the user without admin flags, top them up
        # so the dev experience stays consistent.
        if not user.is_staff or not user.is_superuser:
            user.is_staff = True
            user.is_superuser = True
            user.save(update_fields=["is_staff", "is_superuser"])
        return user
