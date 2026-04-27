"""DRF endpoints for the Selection model.

Endpoints:

- ``POST /api/compounds/selections/`` — create a token-addressed
  snapshot of a compound list. Used by the NLP view when the redirect
  URL would exceed the safe query-string length, and by clients that
  want to persist a selection directly.
- ``GET  /api/compounds/selections/`` — list the requesting chemist's
  selections (saved + non-expired ephemeral).
- ``GET  /api/compounds/selections/<uuid>/`` — retrieve.
- ``PATCH /api/compounds/selections/<uuid>/`` — rename / save / scope.
- ``DELETE /api/compounds/selections/<uuid>/`` — remove.

Authorisation: only the creator can read or modify a selection. Project
sharing via the optional ``target`` field is reserved for a follow-up.
"""

from __future__ import annotations

import datetime
from typing import Optional, Tuple

from django.utils import timezone
from rest_framework import status as http_status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.request import Request
from rest_framework.response import Response

from compounds.registry.models import Selection, Target


# Auto-expire ephemeral selections after this window. Long enough for
# cross-session usefulness ("the queries I ran this morning are still
# there after lunch") without unbounded growth.
DEFAULT_EXPIRY_DAYS = 7


def _default_expiry() -> datetime.datetime:
    return timezone.now() + datetime.timedelta(days=DEFAULT_EXPIRY_DAYS)


def _resolve_target_or_400(
    target_id,
) -> Tuple[Optional[Target], Optional[Response]]:
    """Look up a Target by id. Returns ``(target, None)`` on success or
    ``(None, error_response)`` if the id was provided but unknown. A
    falsy ``target_id`` is treated as "no target" — returns ``(None, None)``."""
    if not target_id:
        return None, None
    target = Target.objects.filter(pk=target_id).first()
    if target is None:
        return None, Response(
            {"error": f"Target {target_id!r} not found"},
            status=http_status.HTTP_400_BAD_REQUEST,
        )
    return target, None


def create_selection(request: Request) -> Response:
    """Create a Selection row from a compound-id list.

    Body:
        compound_ids (list[str], required): formatted compound IDs
            ("NCL-00026007"). Same shape as the existing ``?compound=``
            URL param.
        name (str, optional): display name (typically the NLP scope
            sentence — *"Top 20 CDK4 compounds by scorecard"* etc.).
        source_prompt (str, optional): original NLP prompt for audit.
        is_saved (bool, optional): create a saved selection (no expiry).
        target_id (str, optional): project context.
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

    target, error = _resolve_target_or_400(body.get("target_id"))
    if error is not None:
        return error

    is_saved = bool(body.get("is_saved", False))
    selection = Selection.objects.create(
        name=(body.get("name") or "").strip(),
        compound_ids=compound_ids,
        created_by=request.user,
        expires_at=None if is_saved else _default_expiry(),
        is_saved=is_saved,
        target=target,
        source_prompt=(body.get("source_prompt") or "").strip(),
    )
    return Response(_serialize(selection), status=http_status.HTTP_201_CREATED)


def _serialize(selection: Selection, *, include_compound_ids: bool = False) -> dict:
    """Common shape for the wire format. The list endpoint returns
    every selection the chemist has, so we keep the body lean by
    omitting ``compound_ids`` (they can be 1000s of strings) unless
    explicitly requested. The retrieve endpoint returns them."""
    body = {
        "id": str(selection.id),
        "name": selection.name,
        "n_compounds": len(selection.compound_ids),
        "created_at": selection.created_at.isoformat(),
        "expires_at": (
            selection.expires_at.isoformat() if selection.expires_at else None
        ),
        "is_saved": selection.is_saved,
        "target_id": str(selection.target_id) if selection.target_id else None,
        "target_name": selection.target.name if selection.target_id else None,
        "source_prompt": selection.source_prompt,
    }
    if include_compound_ids:
        body["compound_ids"] = selection.compound_ids
    return body


def get_selection(request: Request, selection_id) -> Response:
    """Retrieve a Selection by token.

    Returns 404 if no such selection or if the requesting user isn't
    the creator (don't leak existence). Returns 410 Gone if expired.
    """
    selection = (
        Selection.objects
        .filter(pk=selection_id, created_by=request.user)
        .select_related("target")
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
        _serialize(selection, include_compound_ids=True),
        status=http_status.HTTP_200_OK,
    )


def list_selections(request: Request) -> Response:
    """List the requesting chemist's selections.

    Returns saved selections + non-expired ephemeral ones, ordered by
    most-recent first. Lean payload: ``compound_ids`` omitted (the
    detail endpoint provides them). Frontend session-list uses this
    to populate the chip strip below the NLP prompt.

    Optional ``?is_saved=true`` filter for "show me only my saved
    queries" view.
    """
    qs = (
        Selection.objects
        .filter(created_by=request.user)
        .select_related("target")
    )

    is_saved_param = request.query_params.get("is_saved")
    if is_saved_param is not None:
        if is_saved_param.lower() in ("true", "1", "yes"):
            qs = qs.filter(is_saved=True)
        elif is_saved_param.lower() in ("false", "0", "no"):
            qs = qs.filter(is_saved=False)

    # Hide expired-and-not-saved rows; cleanup eventually deletes them
    # but the list shouldn't surface them in the meantime.
    qs = qs.exclude(is_saved=False, expires_at__lt=timezone.now())

    qs = qs.order_by("-created_at")[:200]   # generous cap

    return Response(
        {"selections": [_serialize(s) for s in qs]},
        status=http_status.HTTP_200_OK,
    )


def update_selection(request: Request, selection_id) -> Response:
    """Update a Selection.

    Body fields (all optional):
        is_saved (bool): toggles save state. ``True`` clears
            ``expires_at`` (saved selections live forever); ``False``
            re-applies the default expiry window.
        name (str): rename the chip.
        target_id (str|null): set the project context (visibility hint
            for team sharing); null clears it.
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

    body = request.data or {}
    update_fields = []

    if "name" in body:
        new_name = (body.get("name") or "").strip()
        if not new_name:
            return Response(
                {"error": "name must be a non-empty string"},
                status=http_status.HTTP_400_BAD_REQUEST,
            )
        selection.name = new_name
        update_fields.append("name")

    if "is_saved" in body:
        if not isinstance(body["is_saved"], bool):
            return Response(
                {"error": "is_saved must be a boolean"},
                status=http_status.HTTP_400_BAD_REQUEST,
            )
        selection.is_saved = body["is_saved"]
        selection.expires_at = None if selection.is_saved else _default_expiry()
        update_fields.extend(["is_saved", "expires_at"])

    if "target_id" in body:
        target_id = body["target_id"]
        if target_id is None:
            selection.target = None
        else:
            target, error = _resolve_target_or_400(target_id)
            if error is not None:
                return error
            selection.target = target
        update_fields.append("target")

    if update_fields:
        selection.save(update_fields=update_fields)

    return Response(_serialize(selection), status=http_status.HTTP_200_OK)


def delete_selection(request: Request, selection_id) -> Response:
    """Delete a Selection. Creator-scoped."""
    deleted_count, _ = (
        Selection.objects
        .filter(pk=selection_id, created_by=request.user)
        .delete()
    )
    if deleted_count == 0:
        return Response(
            {"error": "Selection not found"},
            status=http_status.HTTP_404_NOT_FOUND,
        )
    return Response(status=http_status.HTTP_204_NO_CONTENT)


# ---------------------------------------------------------------------------
# Method dispatchers — keep /selections/ and /selections/<id>/ as a single
# URL each, dispatching to the right CRUD view by HTTP method.
# ---------------------------------------------------------------------------


@api_view(["GET", "POST"])
@permission_classes([IsAuthenticated])
def selection_collection(request: Request) -> Response:
    """Method dispatcher for /selections/ — GET lists, POST creates."""
    if request.method == "POST":
        return create_selection(request)
    return list_selections(request)


@api_view(["GET", "PATCH", "DELETE"])
@permission_classes([IsAuthenticated])
def selection_detail(request: Request, selection_id) -> Response:
    """Method dispatcher for /selections/<id>/ — GET retrieves, PATCH
    updates (save/rename/scope), DELETE removes."""
    if request.method == "PATCH":
        return update_selection(request, selection_id)
    if request.method == "DELETE":
        return delete_selection(request, selection_id)
    return get_selection(request, selection_id)
