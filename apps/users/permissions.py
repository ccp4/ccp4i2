"""
Permission helpers for CCP4i2/Compounds platform.

Handles the dual-mode authentication:
- Electron (desktop): No auth, everyone is effectively admin
- Web (Docker/Azure): Azure AD auth with admin checks

Role-based authorization:
- admin: Full access (user management, imports, all CRUD)
- contributor: Can add/edit/delete data
- user: Read-only access

Operating level:
- Users can choose to operate at any level up to their assigned role
- Stored in session, defaults to their maximum role
- Allows admins to safely browse as "user" to avoid accidental changes
"""

import os
from django.conf import settings
from rest_framework.permissions import BasePermission

from .models import UserProfile


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


# =============================================================================
# Operating Level Functions
# =============================================================================

def get_user_role(user):
    """
    Get the user's assigned role (maximum authorization level).

    In Electron mode: Returns 'admin'
    In Web mode: Returns the role from UserProfile

    Args:
        user: Django User instance or None

    Returns:
        str: One of 'admin', 'contributor', 'user'
    """
    if not require_auth():
        return UserProfile.ROLE_ADMIN

    if not user or not getattr(user, 'is_authenticated', False):
        return UserProfile.ROLE_USER

    # Check if user is a platform admin (env or superuser)
    user_email = getattr(user, 'email', '').lower()
    if user_email and user_email in get_admin_emails():
        return UserProfile.ROLE_ADMIN
    if getattr(user, 'is_superuser', False):
        return UserProfile.ROLE_ADMIN

    # Get role from profile
    profile = getattr(user, 'profile', None)
    if profile:
        # If is_platform_admin is set, they're an admin
        if profile.is_platform_admin:
            return UserProfile.ROLE_ADMIN
        return profile.role

    return UserProfile.ROLE_USER


def get_operating_level(request):
    """
    Get the user's current operating level from their profile.

    The operating level determines what actions the user can perform.
    It can be lower than their assigned role (e.g., an admin operating as 'user').

    Args:
        request: Django/DRF request object

    Returns:
        str: One of 'admin', 'contributor', 'user'
    """
    if not require_auth():
        return UserProfile.ROLE_ADMIN

    user = getattr(request, 'user', None)
    if not user or not getattr(user, 'is_authenticated', False):
        return UserProfile.ROLE_USER

    # Get from profile, default to their max role
    profile = getattr(user, 'profile', None)
    if not profile:
        return get_user_role(user)

    stored_level = profile.operating_level
    if not stored_level:
        # No operating level set - default to their max role
        return get_user_role(user)

    # Validate the stored level doesn't exceed their role
    user_role = get_user_role(user)
    user_level = UserProfile.ROLE_HIERARCHY.get(user_role, 0)
    operating_level = UserProfile.ROLE_HIERARCHY.get(stored_level, 0)

    if operating_level > user_level:
        # Stored level exceeds their role - reset to max
        return user_role

    return stored_level


def set_operating_level(request, level):
    """
    Set the user's operating level in their profile.

    Args:
        request: Django/DRF request object
        level: One of 'admin', 'contributor', 'user'

    Returns:
        str: The actual level set (may be capped at user's max role)

    Raises:
        ValueError: If level is invalid
    """
    if level not in UserProfile.ROLE_HIERARCHY:
        raise ValueError(f"Invalid operating level: {level}")

    user = getattr(request, 'user', None)
    if not user or not getattr(user, 'is_authenticated', False):
        return UserProfile.ROLE_USER

    # Cap at user's max role
    user_role = get_user_role(user)
    user_level = UserProfile.ROLE_HIERARCHY.get(user_role, 0)
    requested_level = UserProfile.ROLE_HIERARCHY.get(level, 0)

    if requested_level > user_level:
        # Can't elevate above their role
        level = user_role

    # Save to profile
    profile = getattr(user, 'profile', None)
    if profile:
        profile.operating_level = level
        profile.save(update_fields=['operating_level'])

    return level


def can_contribute(request):
    """
    Check if the current operating level allows contributing (add/edit/delete).

    Args:
        request: Django/DRF request object

    Returns:
        bool: True if operating as contributor or admin
    """
    if not require_auth():
        return True

    level = get_operating_level(request)
    return level in (UserProfile.ROLE_CONTRIBUTOR, UserProfile.ROLE_ADMIN)


def can_administer(request):
    """
    Check if the current operating level allows administration.

    Args:
        request: Django/DRF request object

    Returns:
        bool: True if operating as admin
    """
    if not require_auth():
        return True

    return get_operating_level(request) == UserProfile.ROLE_ADMIN


# =============================================================================
# New Permission Classes
# =============================================================================

class IsContributorOrAbove(BasePermission):
    """
    DRF permission class: requires contributor or admin operating level.

    Use this for endpoints that modify data (POST, PUT, PATCH, DELETE).
    """

    message = "Contributor access required. Switch to Contributor or Admin mode to make changes."

    def has_permission(self, request, view):
        return can_contribute(request)


class IsContributorOrReadOnly(BasePermission):
    """
    DRF permission class: contributors can write, users can only read.

    Use this on ViewSets to allow anyone to read but require contributor+ to modify.
    """

    message = "Contributor access required for write operations. Switch to Contributor or Admin mode."

    def has_permission(self, request, view):
        if request.method in ('GET', 'HEAD', 'OPTIONS'):
            return True
        return can_contribute(request)


class IsAdminOrReadOnly(BasePermission):
    """
    DRF permission class: admins can write, others can only read.

    Similar to IsPlatformAdminOrReadOnly but respects operating level.
    """

    message = "Admin access required for write operations. Switch to Admin mode."

    def has_permission(self, request, view):
        if request.method in ('GET', 'HEAD', 'OPTIONS'):
            return True
        return can_administer(request)


class IsContributorCreateAdminUpdate(BasePermission):
    """
    DRF permission class for resources where:
    - Contributors can CREATE new records
    - Only admins can UPDATE or DELETE existing records
    - Anyone can READ

    Use this for resources like Compounds where registration should be open
    to contributors, but editing registered data requires admin approval.
    """

    message = "Admin access required to edit or delete. Contributors can only create new records."

    def has_permission(self, request, view):
        # Read operations always allowed
        if request.method in ('GET', 'HEAD', 'OPTIONS'):
            return True

        # POST (create) allowed for contributors
        if request.method == 'POST':
            return can_contribute(request)

        # PUT, PATCH, DELETE require admin
        return can_administer(request)
