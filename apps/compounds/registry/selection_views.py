"""DRF endpoints for the Selection model (slice 20).

Two endpoints:

- ``POST /api/compounds/selections/`` — create a token-addressed
  snapshot of a compound list. Used by the NLP view when the redirect
  URL would exceed the safe query-string length.
- ``GET /api/compounds/selections/<uuid>/`` — retrieve. The aggregation
  page calls this when it sees ``?selection=<uuid>`` in its URL.

V1 authorisation: only the creator can read a selection. Slice 22 adds
project-scoped sharing.
"""

from __future__ import annotations

import datetime

from django.utils import timezone
from rest_framework import status as http_status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.request import Request
from rest_framework.response import Response

from compounds.registry.models import Selection


# Auto-expire ephemeral selections after this window. Long enough for
# cross-session usefulness ("the queries I ran this morning are still
# there after lunch") without unbounded growth.
DEFAULT_EXPIRY_DAYS = 7


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def create_selection(request: Request) -> Response:
    """Create a Selection row from a compound-id list.

    Body:
        compound_ids (list[str], required): formatted compound IDs
            ("NCL-00026007"). Same shape as the existing ``?compound=``
            URL param.
        name (str, optional): display name (typically the NLP scope
            sentence — *"Top 20 CDK4 compounds by scorecard"* etc.).
        source_prompt (str, optional): original NLP prompt for audit.
    """
    body = request.data or {}

    compound_ids = body.get("compound_ids")
    if not isinstance(compound_ids, list) or not all(
        isinstance(c, str) for c in compound_ids
    ):
        return Response(
            {"error": "compound_ids must be a list of strings"},
            status=http_status.HTTP_400_BAD_REQUEST,
        )

    name = (body.get("name") or "").strip()
    source_prompt = (body.get("source_prompt") or "").strip()

    selection = Selection.objects.create(
        name=name,
        compound_ids=compound_ids,
        created_by=request.user,
        expires_at=timezone.now() + datetime.timedelta(days=DEFAULT_EXPIRY_DAYS),
        source_prompt=source_prompt,
    )
    return Response(
        {
            "id": str(selection.id),
            "name": selection.name,
            "n_compounds": len(selection.compound_ids),
            "created_at": selection.created_at.isoformat(),
            "expires_at": selection.expires_at.isoformat() if selection.expires_at else None,
        },
        status=http_status.HTTP_201_CREATED,
    )


@api_view(["GET"])
@permission_classes([IsAuthenticated])
def get_selection(request: Request, selection_id) -> Response:
    """Retrieve a Selection by token.

    Returns 404 if no such selection or if the requesting user isn't
    the creator (don't leak existence). Returns 410 Gone if expired.
    """
    selection = (
        Selection.objects
        .filter(pk=selection_id, created_by=request.user)
        .first()
    )
    if selection is None:
        return Response(
            {"error": "Selection not found"},
            status=http_status.HTTP_404_NOT_FOUND,
        )
    if selection.is_expired:
        return Response(
            {"error": "Selection has expired"},
            status=http_status.HTTP_410_GONE,
        )
    return Response(
        {
            "id": str(selection.id),
            "name": selection.name,
            "compound_ids": selection.compound_ids,
            "n_compounds": len(selection.compound_ids),
            "created_at": selection.created_at.isoformat(),
            "expires_at": selection.expires_at.isoformat() if selection.expires_at else None,
            "source_prompt": selection.source_prompt,
        },
        status=http_status.HTTP_200_OK,
    )
