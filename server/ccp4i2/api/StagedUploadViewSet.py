"""Chunked staging transport for large-file import on web deployments.

Endpoints (owner-bound, IsAuthenticated), under the API router:

    POST   staged-uploads/                     begin  -> {upload_id, chunk_bytes, expires_at}
    GET    staged-uploads/<uuid>/              status -> {state, present_chunks, ...}
    PUT    staged-uploads/<uuid>/chunks/<i>    a raw-bytes chunk  -> 204
    POST   staged-uploads/<uuid>/finish/       assemble + verify  -> {state: ready}

Every refusal carries a distinct status code (413/429/404/410/409/422) so the
client can report honestly, rather than the a63 silent fall-through to an empty
``request.FILES``. The lifecycle logic lives in
``lib/utils/files/staged_upload.py``; this is the thin HTTP layer.
"""

from datetime import timedelta

from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.parsers import BaseParser
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from ..lib.utils.files import staged_upload as su


class RawBytesParser(BaseParser):
    """Pass an ``application/octet-stream`` body through untouched, as bytes."""
    media_type = "application/octet-stream"

    def parse(self, stream, media_type=None, parser_context=None):
        return stream.read()


_owner = su.owner_key


def _error(exc):
    return Response({"error": exc.message}, status=exc.status)


class StagedUploadViewSet(viewsets.ViewSet):
    permission_classes = [IsAuthenticated]
    # The URL id is the row's uuid; we filter by (uuid, owner) ourselves rather
    # than get_object(), so a foreign id is a 404 like an unknown one.
    lookup_field = "uuid"

    def create(self, request):
        """begin: reserve a staging directory and handle."""
        try:
            row = su.begin(
                owner=_owner(request),
                filename=request.data.get("filename"),
                size_bytes=request.data.get("size_bytes"),
                sha256=request.data.get("sha256", ""),
            )
        except su.StagedUploadError as exc:
            return _error(exc)
        expires = row.created_at + timedelta(
            hours=su.ttl_hours())
        return Response(
            {
                "upload_id": str(row.uuid),
                "chunk_bytes": su.chunk_bytes(),
                "expires_at": expires.isoformat(),
            },
            status=201,
        )

    def retrieve(self, request, uuid=None):
        """status: state + which chunks are present (for resume)."""
        try:
            row = su.get_owned(uuid, _owner(request))
        except su.StagedUploadError as exc:
            return _error(exc)
        return Response({
            "upload_id": str(row.uuid),
            "state": row.state,
            "size_bytes": row.size_bytes,
            "present_chunks": su.present_indexes(row),
        })

    @action(detail=True, methods=["put"],
            url_path="chunks/(?P<index>[0-9]+)",
            parser_classes=[RawBytesParser])
    def chunk(self, request, uuid=None, index=None):
        """Write one chunk from the raw request body."""
        try:
            row = su.get_owned(uuid, _owner(request))
            data = request.data
            if not isinstance(data, (bytes, bytearray)):
                data = request.body
            su.write_chunk(row, int(index), bytes(data))
        except su.StagedUploadError as exc:
            return _error(exc)
        return Response(status=204)

    @action(detail=True, methods=["post"], url_path="finish")
    def finish(self, request, uuid=None):
        """Assemble, verify, and mark ready."""
        try:
            row = su.finish(su.get_owned(uuid, _owner(request)))
        except su.StagedUploadError as exc:
            return _error(exc)
        return Response({"upload_id": str(row.uuid), "state": row.state})
