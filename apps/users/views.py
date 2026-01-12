"""API views for users app."""

from django.contrib.auth import get_user_model
from rest_framework import status
from rest_framework.decorators import action
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework.viewsets import ModelViewSet

from .models import UserProfile
from .permissions import IsPlatformAdmin, is_platform_admin, require_auth
from .serializers import (
    CurrentUserSerializer,
    UserListSerializer,
    UserSerializer,
)

User = get_user_model()


class CurrentUserView(APIView):
    """
    Get current authenticated user info.

    In Electron mode: Returns a default local user
    In Web mode: Returns the Azure AD authenticated user
    """

    def get_permissions(self):
        if require_auth():
            return [IsAuthenticated()]
        return []

    def get(self, request):
        if not require_auth():
            # Electron mode: return a synthetic local user
            return Response({
                'id': 0,
                'username': 'local_user',
                'email': '',
                'first_name': 'Local',
                'last_name': 'User',
                'display_name': 'Local User',
                'is_admin': True,
                'profile': {
                    'is_platform_admin': True,
                    'legacy_username': '',
                    'legacy_display_name': '',
                    'imported_at': None,
                    'first_login_at': None,
                    'last_seen_at': None,
                },
            })

        # Web mode: return actual user
        user = request.user
        # Record activity
        if hasattr(user, 'profile'):
            user.profile.record_activity()

        serializer = CurrentUserSerializer(user)
        return Response(serializer.data)


class UserViewSet(ModelViewSet):
    """
    ViewSet for user management.

    Only accessible by platform admins.
    """

    queryset = User.objects.select_related('profile').all()
    permission_classes = [IsAuthenticated, IsPlatformAdmin]

    def get_serializer_class(self):
        if self.action == 'list':
            return UserListSerializer
        return UserSerializer

    def get_queryset(self):
        queryset = super().get_queryset()
        # Filter by active status if specified
        is_active = self.request.query_params.get('is_active')
        if is_active is not None:
            queryset = queryset.filter(is_active=is_active.lower() == 'true')
        # Filter by admin status if specified
        is_admin = self.request.query_params.get('is_admin')
        if is_admin is not None:
            queryset = queryset.filter(profile__is_platform_admin=is_admin.lower() == 'true')
        return queryset.order_by('username')

    @action(detail=True, methods=['post'])
    def grant_admin(self, request, pk=None):
        """Grant platform admin to a user."""
        user = self.get_object()
        UserProfile.objects.update_or_create(
            user=user,
            defaults={'is_platform_admin': True}
        )
        return Response({'status': 'admin granted', 'user': user.email})

    @action(detail=True, methods=['post'])
    def revoke_admin(self, request, pk=None):
        """Revoke platform admin from a user."""
        user = self.get_object()
        # Prevent revoking your own admin
        if user == request.user:
            return Response(
                {'error': 'Cannot revoke your own admin status'},
                status=status.HTTP_400_BAD_REQUEST
            )
        UserProfile.objects.update_or_create(
            user=user,
            defaults={'is_platform_admin': False}
        )
        return Response({'status': 'admin revoked', 'user': user.email})


class AdminListView(APIView):
    """List all platform admins."""

    permission_classes = [IsAuthenticated]

    def get(self, request):
        if not require_auth():
            return Response([])

        admins = User.objects.filter(
            profile__is_platform_admin=True
        ).select_related('profile')

        serializer = UserListSerializer(admins, many=True)
        return Response(serializer.data)
