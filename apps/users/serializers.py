"""DRF serializers for users app."""

from django.contrib.auth import get_user_model
from rest_framework import serializers

from .models import UserProfile
from .permissions import is_platform_admin

User = get_user_model()


class UserProfileSerializer(serializers.ModelSerializer):
    """Serializer for UserProfile model."""

    class Meta:
        model = UserProfile
        fields = [
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
    """Lightweight serializer for user lists."""

    is_admin = serializers.SerializerMethodField()
    display_name = serializers.SerializerMethodField()

    class Meta:
        model = User
        fields = ['id', 'username', 'email', 'display_name', 'is_admin', 'is_active']

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
            'profile',
        ]

    def get_is_admin(self, obj):
        return is_platform_admin(obj)

    def get_display_name(self, obj):
        if obj.first_name and obj.last_name:
            return f"{obj.first_name} {obj.last_name}"
        return obj.username
