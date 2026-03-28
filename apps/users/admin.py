# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""Django admin configuration for users app."""

from django.contrib import admin
from django.contrib.auth import get_user_model
from django.contrib.auth.admin import UserAdmin as BaseUserAdmin

from .models import UserProfile

User = get_user_model()


class UserProfileInline(admin.StackedInline):
    """Inline admin for UserProfile on User admin page."""
    model = UserProfile
    can_delete = False
    verbose_name_plural = 'Profile'
    fk_name = 'user'


class UserAdmin(BaseUserAdmin):
    """Extended User admin with profile inline."""
    inlines = (UserProfileInline,)
    list_display = (
        'username', 'email', 'first_name', 'last_name',
        'is_staff', 'get_is_platform_admin', 'date_joined'
    )
    list_filter = BaseUserAdmin.list_filter + ('profile__is_platform_admin',)

    def get_is_platform_admin(self, obj):
        return getattr(getattr(obj, 'profile', None), 'is_platform_admin', False)
    get_is_platform_admin.boolean = True
    get_is_platform_admin.short_description = 'Platform Admin'


@admin.register(UserProfile)
class UserProfileAdmin(admin.ModelAdmin):
    """Admin for UserProfile model."""
    list_display = (
        'user', 'is_platform_admin', 'legacy_username',
        'imported_at', 'first_login_at', 'last_seen_at'
    )
    list_filter = ('is_platform_admin', 'imported_at')
    search_fields = ('user__username', 'user__email', 'legacy_username')
    readonly_fields = ('imported_at', 'first_login_at', 'last_seen_at')
    raw_id_fields = ('user',)


# Re-register User with extended admin
try:
    admin.site.unregister(User)
except admin.sites.NotRegistered:
    pass
admin.site.register(User, UserAdmin)
