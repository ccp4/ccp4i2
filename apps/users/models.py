"""
User profile and admin management models.

This module extends Django's auth.User with additional profile information
for the CCP4i2/Compounds platform, including:
- Platform admin designation
- Legacy import tracking
- Login history
"""

from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.utils import timezone


class UserProfile(models.Model):
    """
    Extended user profile for CCP4i2/Compounds platform.

    Automatically created for each User via signal.
    """

    # Role choices for authorization levels
    ROLE_USER = 'user'
    ROLE_CONTRIBUTOR = 'contributor'
    ROLE_ADMIN = 'admin'

    ROLE_CHOICES = [
        (ROLE_USER, 'User (read-only)'),
        (ROLE_CONTRIBUTOR, 'Contributor (can add/edit/delete)'),
        (ROLE_ADMIN, 'Admin (full access)'),
    ]

    # Role hierarchy for comparison
    ROLE_HIERARCHY = {
        ROLE_USER: 0,
        ROLE_CONTRIBUTOR: 1,
        ROLE_ADMIN: 2,
    }

    user = models.OneToOneField(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        related_name='profile',
        primary_key=True,
    )

    # Role-based authorization
    role = models.CharField(
        max_length=20,
        choices=ROLE_CHOICES,
        default=ROLE_CONTRIBUTOR,
        help_text="Maximum authorization level for this user"
    )

    # Current operating level (can be <= role for self-imposed restrictions)
    operating_level = models.CharField(
        max_length=20,
        choices=ROLE_CHOICES,
        null=True,
        blank=True,
        help_text="Current operating level (defaults to role if not set)"
    )

    # Admin designation (for web mode) - DEPRECATED, use role='admin' instead
    # Kept for backwards compatibility during migration
    is_platform_admin = models.BooleanField(
        default=False,
        help_text="DEPRECATED: Use role field instead. Can manage platform settings, import data, and administer users"
    )

    # Legacy import tracking
    legacy_username = models.CharField(
        max_length=100,
        blank=True,
        help_text="Username from legacy system (for reference)"
    )
    legacy_display_name = models.CharField(
        max_length=255,
        blank=True,
        help_text="Display name from legacy system"
    )
    imported_at = models.DateTimeField(
        null=True,
        blank=True,
        help_text="When this user was imported from legacy fixtures"
    )

    # Login tracking
    first_login_at = models.DateTimeField(
        null=True,
        blank=True,
        help_text="When user first logged in via Azure AD"
    )
    last_seen_at = models.DateTimeField(
        null=True,
        blank=True,
        help_text="Last activity timestamp"
    )

    class Meta:
        verbose_name = "User Profile"
        verbose_name_plural = "User Profiles"

    def __str__(self):
        return f"Profile for {self.user.username}"

    def record_login(self):
        """Record a login event."""
        now = timezone.now()
        if self.first_login_at is None:
            self.first_login_at = now
        self.last_seen_at = now
        self.save(update_fields=['first_login_at', 'last_seen_at'])

    def record_activity(self):
        """Update last seen timestamp."""
        self.last_seen_at = timezone.now()
        self.save(update_fields=['last_seen_at'])

    def has_role(self, required_role):
        """
        Check if user has at least the required role level.

        Args:
            required_role: One of ROLE_USER, ROLE_CONTRIBUTOR, ROLE_ADMIN

        Returns:
            bool: True if user's role is >= required_role
        """
        user_level = self.ROLE_HIERARCHY.get(self.role, 0)
        required_level = self.ROLE_HIERARCHY.get(required_role, 0)
        return user_level >= required_level

    def is_admin(self):
        """Check if user has admin role."""
        return self.role == self.ROLE_ADMIN or self.is_platform_admin

    def is_contributor_or_above(self):
        """Check if user can contribute (add/edit/delete)."""
        return self.has_role(self.ROLE_CONTRIBUTOR) or self.is_platform_admin


@receiver(post_save, sender=settings.AUTH_USER_MODEL)
def create_user_profile(sender, instance, created, **kwargs):
    """Create UserProfile when User is created."""
    if created:
        UserProfile.objects.get_or_create(user=instance)


@receiver(post_save, sender=settings.AUTH_USER_MODEL)
def save_user_profile(sender, instance, **kwargs):
    """Ensure UserProfile exists and is saved."""
    UserProfile.objects.get_or_create(user=instance)
