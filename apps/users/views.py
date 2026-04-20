"""API views for users app."""

from django.contrib.auth import get_user_model
from django.utils import timezone
from rest_framework import status
from rest_framework.decorators import action
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework.viewsets import ModelViewSet

from .models import PendingInvite, UserProfile
from .permissions import (
    IsPlatformAdmin,
    is_platform_admin,
    require_auth,
    get_user_role,
    get_operating_level,
    set_operating_level,
)
from .serializers import (
    CurrentUserSerializer,
    PendingInviteSerializer,
    UserListSerializer,
    UserSerializer,
)

User = get_user_model()


class CurrentUserView(APIView):
    """
    Get current authenticated user info.

    In Electron mode: Returns a default local user with full admin access
    In Web mode: Returns the Azure AD authenticated user with role/operating level

    Response includes:
    - role: The user's maximum authorized role
    - operating_level: Current session operating level (can be <= role)
    - can_contribute: Whether user can add/edit/delete at current level
    - can_administer: Whether user can perform admin actions at current level
    """

    def get_permissions(self):
        if require_auth():
            return [IsAuthenticated()]
        return []

    def get(self, request):
        if not require_auth():
            # Electron mode: return a synthetic local user with full access
            return Response({
                'id': 0,
                'username': 'local_user',
                'email': '',
                'first_name': 'Local',
                'last_name': 'User',
                'display_name': 'Local User',
                'is_admin': True,
                'role': UserProfile.ROLE_ADMIN,
                'operating_level': UserProfile.ROLE_ADMIN,
                'can_contribute': True,
                'can_administer': True,
                'profile': {
                    'role': UserProfile.ROLE_ADMIN,
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

        serializer = CurrentUserSerializer(user, context={'request': request})
        return Response(serializer.data)


class UserViewSet(ModelViewSet):
    """
    ViewSet for user management.

    Only accessible by platform admins.
    Note: IsPlatformAdmin handles the no-auth case (returns True when
    CCP4I2_REQUIRE_AUTH is not set), so we don't need IsAuthenticated here.
    """

    queryset = User.objects.select_related('profile').all()
    permission_classes = [IsPlatformAdmin]

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

    @action(detail=True, methods=['post'])
    def set_role(self, request, pk=None):
        """
        Set the role for a user.

        POST data:
        - role: One of 'admin', 'contributor', 'user'
        """
        user = self.get_object()
        role = request.data.get('role')

        if not role:
            return Response(
                {'error': 'role is required'},
                status=status.HTTP_400_BAD_REQUEST
            )

        if role not in UserProfile.ROLE_HIERARCHY:
            return Response(
                {'error': f'Invalid role. Must be one of: {list(UserProfile.ROLE_HIERARCHY.keys())}'},
                status=status.HTTP_400_BAD_REQUEST
            )

        # Prevent demoting yourself
        if user == request.user and role != UserProfile.ROLE_ADMIN:
            return Response(
                {'error': 'Cannot demote your own role'},
                status=status.HTTP_400_BAD_REQUEST
            )

        # Update role and sync is_platform_admin
        is_admin = (role == UserProfile.ROLE_ADMIN)
        UserProfile.objects.update_or_create(
            user=user,
            defaults={
                'role': role,
                'is_platform_admin': is_admin,
            }
        )

        return Response({
            'status': 'role updated',
            'user': user.email,
            'role': role,
        })


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


class PendingInviteViewSet(ModelViewSet):
    """
    ViewSet for queued B2B invites.

    Admins POST {email, note} to queue an invite, GET the list, DELETE a row
    to cancel a pending invite. A privileged operator drains the queue by
    running Docker/azure-uksouth/scripts/invite-user.sh, which also marks rows
    as sent via the `mark_sent` action.
    """

    queryset = PendingInvite.objects.select_related('requested_by', 'sent_by').all()
    serializer_class = PendingInviteSerializer
    permission_classes = [IsPlatformAdmin]

    def get_queryset(self):
        qs = super().get_queryset()
        status_filter = self.request.query_params.get('status')
        if status_filter:
            qs = qs.filter(status=status_filter)
        return qs

    def perform_create(self, serializer):
        user = self.request.user if getattr(self.request, 'user', None) and self.request.user.is_authenticated else None
        serializer.save(requested_by=user)

    @action(detail=True, methods=['post'])
    def mark_sent(self, request, pk=None):
        """Mark an invite as sent (called by invite-user.sh after Graph succeeds)."""
        invite = self.get_object()
        invite.status = PendingInvite.STATUS_SENT
        invite.sent_at = timezone.now()
        invite.sent_by = request.user if getattr(request, 'user', None) and request.user.is_authenticated else None
        invite.guest_object_id = request.data.get('guest_object_id', '')
        invite.failure_reason = ''
        invite.save()
        return Response(PendingInviteSerializer(invite).data)

    @action(detail=True, methods=['post'])
    def mark_failed(self, request, pk=None):
        """Mark an invite as failed with a reason (called by invite-user.sh)."""
        invite = self.get_object()
        invite.status = PendingInvite.STATUS_FAILED
        invite.failure_reason = request.data.get('reason', '')[:500]
        invite.save()
        return Response(PendingInviteSerializer(invite).data)


class OperatingLevelView(APIView):
    """
    Get or set the user's operating level for this session.

    The operating level controls what actions a user can perform:
    - 'admin': Full access
    - 'contributor': Can add/edit/delete data
    - 'user': Read-only access

    Users can only set their operating level to values <= their assigned role.
    """

    def get_permissions(self):
        if require_auth():
            return [IsAuthenticated()]
        return []

    def get(self, request):
        """Get current operating level and available options."""
        if not require_auth():
            return Response({
                'operating_level': UserProfile.ROLE_ADMIN,
                'role': UserProfile.ROLE_ADMIN,
                'available_levels': [UserProfile.ROLE_ADMIN],
            })

        user_role = get_user_role(request.user)
        operating_level = get_operating_level(request)

        # Get available levels (all levels up to and including user's role)
        user_level = UserProfile.ROLE_HIERARCHY.get(user_role, 0)
        available_levels = [
            role for role, level in UserProfile.ROLE_HIERARCHY.items()
            if level <= user_level
        ]

        return Response({
            'operating_level': operating_level,
            'role': user_role,
            'available_levels': available_levels,
        })

    def post(self, request):
        """Set operating level for this session."""
        if not require_auth():
            return Response({
                'operating_level': UserProfile.ROLE_ADMIN,
                'message': 'Operating level not applicable in Electron mode',
            })

        level = request.data.get('level')
        if not level:
            return Response(
                {'error': 'level is required'},
                status=status.HTTP_400_BAD_REQUEST
            )

        if level not in UserProfile.ROLE_HIERARCHY:
            return Response(
                {'error': f'Invalid level. Must be one of: {list(UserProfile.ROLE_HIERARCHY.keys())}'},
                status=status.HTTP_400_BAD_REQUEST
            )

        actual_level = set_operating_level(request, level)

        # Check if level was capped
        if actual_level != level:
            return Response({
                'operating_level': actual_level,
                'message': f'Level capped to your maximum role: {actual_level}',
            })

        return Response({
            'operating_level': actual_level,
            'message': f'Now operating as: {actual_level}',
        })
