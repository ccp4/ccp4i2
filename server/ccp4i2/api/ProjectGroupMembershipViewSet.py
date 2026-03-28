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
"""
ViewSet for managing project group memberships directly.

While most membership operations should go through ProjectGroupViewSet actions
(add_member, remove_member, set_parent), this viewset provides direct CRUD
access for advanced use cases and admin operations.
"""
import logging
from rest_framework.viewsets import ModelViewSet
from rest_framework.parsers import FormParser, MultiPartParser, JSONParser
from rest_framework.permissions import IsAuthenticated
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
    permission_classes = [IsAuthenticated]

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
