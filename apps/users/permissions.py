"""
Permission helpers for CCP4i2/Compounds platform.

Handles the dual-mode authentication:
- Electron (desktop): No auth, everyone is effectively admin
- Web (Docker/Azure): Azure AD auth with admin checks
"""

import os
from django.conf import settings
from rest_framework.permissions import BasePermission


def require_auth():
    """Check if authentication is required (web mode)."""
    return os.environ.get("CCP4I2_REQUIRE_AUTH", "").lower() in ("true", "1", "yes")


def get_admin_emails():
    """Get list of bootstrap admin emails from environment."""
    emails_str = getattr(settings, 'PLATFORM_ADMIN_EMAILS', None)
    if emails_str is None:
        emails_str = os.environ.get('PLATFORM_ADMIN_EMAILS', '')
    if isinstance(emails_str, str):
        return [e.strip().lower() for e in emails_str.split(',') if e.strip()]
    return [e.lower() for e in emails_str]


def is_platform_admin(user):
    """
    Check if user is a platform admin.

    In Electron mode (no auth required): Always returns True
    In Web mode: Checks environment config and database flag

    Args:
        user: Django User instance or None

    Returns:
        bool: True if user has admin privileges
    """
    # Electron mode: no restrictions
    if not require_auth():
        return True

    # Web mode: check actual permissions
    if not user or not getattr(user, 'is_authenticated', False):
        return False

    # Check 1: Environment-based bootstrap admins
    user_email = getattr(user, 'email', '').lower()
    if user_email and user_email in get_admin_emails():
        return True

    # Check 2: Database flag on profile
    profile = getattr(user, 'profile', None)
    if profile and getattr(profile, 'is_platform_admin', False):
        return True

    # Check 3: Django superuser/staff (fallback)
    if getattr(user, 'is_superuser', False):
        return True

    return False


class IsPlatformAdmin(BasePermission):
    """
    DRF permission class for platform admin access.

    Usage:
        class MyAdminView(APIView):
            permission_classes = [IsAuthenticated, IsPlatformAdmin]
    """

    message = "Platform admin access required."

    def has_permission(self, request, view):
        return is_platform_admin(request.user)


class IsPlatformAdminOrReadOnly(BasePermission):
    """
    DRF permission class: admins can write, others can only read.
    """

    message = "Platform admin access required for write operations."

    def has_permission(self, request, view):
        if request.method in ('GET', 'HEAD', 'OPTIONS'):
            return True
        return is_platform_admin(request.user)
