"""DRF serializers for users app."""

from django.contrib.auth import get_user_model
from rest_framework import serializers

from .models import UserProfile
from .permissions import is_platform_admin, get_user_role, get_operating_level

User = get_user_model()


class UserProfileSerializer(serializers.ModelSerializer):
    """Serializer for UserProfile model."""

    class Meta:
        model = UserProfile
        fields = [
            'role',
            'is_platform_admin',
            'legacy_username',
            'legacy_display_name',
            'imported_at',
            'first_login_at',
            'last_seen_at',
        ]
        read_only_fields = ['imported_at', 'first_login_at', 'last_seen_at']


class UserSerializer(serializers.ModelSerializer):
    """Serializer for User model with profile."""

    profile = UserProfileSerializer(read_only=True)
    is_admin = serializers.SerializerMethodField()

    class Meta:
        model = User
        fields = [
            'id',
            'username',
            'email',
            'first_name',
            'last_name',
            'is_active',
            'date_joined',
            'last_login',
            'profile',
            'is_admin',
        ]
        read_only_fields = ['id', 'date_joined', 'last_login']

    def get_is_admin(self, obj):
        return is_platform_admin(obj)


class UserListSerializer(serializers.ModelSerializer):
    """Serializer for user lists with profile role information."""

    is_admin = serializers.SerializerMethodField()
    display_name = serializers.SerializerMethodField()
    profile = UserProfileSerializer(read_only=True)

    class Meta:
        model = User
        fields = ['id', 'username', 'email', 'display_name', 'is_admin', 'is_active', 'profile']

    def get_is_admin(self, obj):
        return is_platform_admin(obj)

    def get_display_name(self, obj):
        if obj.first_name and obj.last_name:
            return f"{obj.first_name} {obj.last_name}"
        return obj.username


class CurrentUserSerializer(serializers.ModelSerializer):
    """Serializer for the current authenticated user."""

    profile = UserProfileSerializer(read_only=True)
    is_admin = serializers.SerializerMethodField()
    display_name = serializers.SerializerMethodField()
    role = serializers.SerializerMethodField()
    operating_level = serializers.SerializerMethodField()
    can_contribute = serializers.SerializerMethodField()
    can_administer = serializers.SerializerMethodField()

    class Meta:
        model = User
        fields = [
            'id',
            'username',
            'email',
            'first_name',
            'last_name',
            'display_name',
            'is_admin',
            'role',
            'operating_level',
            'can_contribute',
            'can_administer',
            'profile',
        ]

    def get_is_admin(self, obj):
        return is_platform_admin(obj)

    def get_display_name(self, obj):
        if obj.first_name and obj.last_name:
            return f"{obj.first_name} {obj.last_name}"
        return obj.username

    def get_role(self, obj):
        """Get the user's maximum authorized role."""
        return get_user_role(obj)

    def get_operating_level(self, obj):
        """Get the user's current operating level from session."""
        request = self.context.get('request')
        if request:
            return get_operating_level(request)
        return get_user_role(obj)

    def get_can_contribute(self, obj):
        """Check if user can currently add/edit/delete."""
        request = self.context.get('request')
        if request:
            level = get_operating_level(request)
            return level in (UserProfile.ROLE_CONTRIBUTOR, UserProfile.ROLE_ADMIN)
        return False

    def get_can_administer(self, obj):
        """Check if user can currently perform admin actions."""
        request = self.context.get('request')
        if request:
            return get_operating_level(request) == UserProfile.ROLE_ADMIN
        return False
