import logging

from rest_framework.decorators import action
from rest_framework.exceptions import ValidationError
from rest_framework.parsers import FormParser, JSONParser, MultiPartParser
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response
from rest_framework import status
from rest_framework.viewsets import ModelViewSet

from . import serializers
from ..db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


class ProjectTagViewSet(ModelViewSet):
    queryset = models.ProjectTag.objects.all()
    serializer_class = serializers.ProjectTagSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    permission_classes = [IsAuthenticated]

    def handle_exception(self, exc):
        """Give a serializer rejection a single ``error`` string as well.

        DRF reports a ``ValidationError`` per field — ``{"parent": ["..."]}`` —
        and the renderer's fetch wrapper only surfaces a body's ``error`` key,
        so "a tag cannot be moved beneath itself" reached the user as
        "HTTP 400: Bad Request". The per-field shape is kept for anything that
        wants to attribute the message to a field.
        """
        response = super().handle_exception(exc)
        if isinstance(exc, ValidationError) and isinstance(response.data, dict):
            messages = []
            for value in response.data.values():
                for message in value if isinstance(value, list) else [value]:
                    messages.append(str(message))
            response.data.setdefault("error", " ".join(messages))
        return response

    @action(detail=False, methods=["get"])
    def tree(self, request):
        """The whole tag forest, nested, with project counts.

        Two counts are reported per node: ``project_count`` is what is filed
        directly on that tag, ``total_project_count`` is the distinct projects
        anywhere in its subtree. They are distinct numbers rather than a sum,
        because a project tagged both ``Mpro`` and ``Mpro/soaks`` must not be
        counted twice in the ancestor.

        The whole forest is returned in three queries. Tag counts are desktop
        scale — a user's project list, not a warehouse — so exact roll-up in
        memory beats an approximate one in SQL.
        """
        tags = list(models.ProjectTag.objects.all().order_by("path"))
        through = models.ProjectTag.projects.through
        pairs = through.objects.values_list("projecttag_id", "project_id")

        direct = {tag.id: set() for tag in tags}
        for tag_id, project_id in pairs:
            if tag_id in direct:
                direct[tag_id].add(project_id)

        children_of = {}
        for tag in tags:
            children_of.setdefault(tag.parent_id, []).append(tag)

        # Ordering by path puts every parent before its children, so walking it
        # backwards guarantees a node's subtree is aggregated before the node.
        subtree = {}
        for tag in reversed(tags):
            collected = set(direct[tag.id])
            for child in children_of.get(tag.id, []):
                collected |= subtree[child.id]
            subtree[tag.id] = collected

        def node(tag):
            return {
                "id": tag.id,
                "text": tag.text,
                "parent": tag.parent_id,
                "path": tag.path,
                "display_path": tag.display_path,
                "depth": tag.depth,
                "project_count": len(direct[tag.id]),
                "total_project_count": len(subtree[tag.id]),
                "children": [node(child) for child in children_of.get(tag.id, [])],
            }

        forest = [node(tag) for tag in children_of.get(None, [])]
        untagged = models.Project.objects.filter(tags__isnull=True).count()
        return Response({"tags": forest, "untagged_project_count": untagged})

    @action(detail=True, methods=["post"])
    def add_projects(self, request, pk=None):
        """Tag several projects at once — what a multi-selection or a drop onto
        a tree node needs, instead of one request per project."""
        return self._bulk_membership(request, pk, add=True)

    @action(detail=True, methods=["post"])
    def remove_projects(self, request, pk=None):
        return self._bulk_membership(request, pk, add=False)

    def _bulk_membership(self, request, pk, add):
        try:
            tag = models.ProjectTag.objects.get(pk=pk)
        except models.ProjectTag.DoesNotExist:
            return Response(
                {"error": "Tag not found"}, status=status.HTTP_404_NOT_FOUND
            )

        project_ids = request.data.get("project_ids")
        if not isinstance(project_ids, list) or not project_ids:
            return Response(
                {"error": "project_ids must be a non-empty list"},
                status=status.HTTP_400_BAD_REQUEST,
            )

        projects = list(models.Project.objects.filter(pk__in=project_ids))
        found = {project.pk for project in projects}
        missing = [pk_ for pk_ in project_ids if pk_ not in found]

        if add:
            tag.projects.add(*projects)
        else:
            tag.projects.remove(*projects)

        return Response(
            {
                "status": "success",
                "tag": tag.id,
                "changed": len(projects),
                "missing": missing,
            },
            status=status.HTTP_200_OK,
        )
