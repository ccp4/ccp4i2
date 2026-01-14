"""
ViewSet for managing project group memberships directly.

While most membership operations should go through ProjectGroupViewSet actions
(add_member, remove_member, set_parent), this viewset provides direct CRUD
access for advanced use cases and admin operations.
"""
import logging
from rest_framework.viewsets import ModelViewSet
from rest_framework.parsers import FormParser, MultiPartParser, JSONParser
from . import serializers
from ..db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


class ProjectGroupMembershipViewSet(ModelViewSet):
    """
    ViewSet for managing project group memberships.

    Provides standard CRUD operations for ProjectGroupMembership model.
    Supports filtering by group and membership type.

    Query Parameters:
        - group: Filter by group ID
        - type: Filter by membership type ('parent' or 'member')
        - project: Filter by project ID
    """

    queryset = models.ProjectGroupMembership.objects.select_related(
        "group", "project"
    ).all()
    serializer_class = serializers.ProjectGroupMembershipSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]

    def get_queryset(self):
        """Filter by group, type, or project if provided in query params."""
        queryset = super().get_queryset()

        group_id = self.request.query_params.get("group")
        if group_id:
            queryset = queryset.filter(group_id=group_id)

        membership_type = self.request.query_params.get("type")
        if membership_type:
            queryset = queryset.filter(type=membership_type)

        project_id = self.request.query_params.get("project")
        if project_id:
            queryset = queryset.filter(project_id=project_id)

        return queryset
