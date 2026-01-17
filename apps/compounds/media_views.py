"""
Protected media file serving for compounds app.

These views serve files from Django FileFields through authenticated endpoints,
ensuring that sensitive data (assay files, QC documents, etc.) are not publicly accessible.

Authentication is handled by the AzureADAuthMiddleware - these views only need to
ensure they're not added to exempt paths.
"""

import os
import mimetypes
from pathlib import Path

from django.conf import settings
from django.http import FileResponse, Http404, HttpResponseBadRequest
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated

from compounds.assays.models import Assay, ProtocolDocument
from compounds.registry.models import BatchQCFile
from compounds.constructs.models import Plasmid, CassetteUse, SequencingResult


def _serve_file(file_field, filename_override=None):
    """
    Serve a file from a FileField with proper headers.

    Args:
        file_field: Django FileField instance
        filename_override: Optional filename for Content-Disposition header

    Returns:
        FileResponse with appropriate content type and disposition
    """
    if not file_field:
        raise Http404("File not found")

    # Get the file path
    file_path = file_field.path

    if not os.path.exists(file_path):
        raise Http404("File not found on disk")

    # Determine content type
    content_type, _ = mimetypes.guess_type(file_path)
    if not content_type:
        content_type = 'application/octet-stream'

    # Get filename for Content-Disposition
    filename = filename_override or os.path.basename(file_field.name)

    # Create response
    response = FileResponse(
        open(file_path, 'rb'),
        content_type=content_type,
    )

    # Set Content-Disposition for download
    # Use inline for common viewable types, attachment for others
    viewable_types = [
        'application/pdf',
        'text/plain',
        'text/csv',
        'image/png',
        'image/jpeg',
        'image/gif',
    ]

    disposition = 'inline' if content_type in viewable_types else 'attachment'
    response['Content-Disposition'] = f'{disposition}; filename="{filename}"'

    return response


@api_view(['GET'])
@permission_classes([IsAuthenticated])
def serve_assay_data_file(request, assay_id):
    """
    Serve an assay's data file (Excel, CSV, etc.).

    URL: /api/compounds/media/assays/<assay_id>/data_file/
    """
    try:
        assay = Assay.objects.get(id=assay_id)
    except Assay.DoesNotExist:
        raise Http404("Assay not found")

    if not assay.data_file:
        raise Http404("Assay has no data file")

    return _serve_file(assay.data_file, assay.data_filename)


@api_view(['GET'])
@permission_classes([IsAuthenticated])
def serve_protocol_document(request, document_id):
    """
    Serve a protocol document file.

    URL: /api/compounds/media/protocol-documents/<document_id>/file/
    """
    try:
        doc = ProtocolDocument.objects.get(id=document_id)
    except ProtocolDocument.DoesNotExist:
        raise Http404("Protocol document not found")

    if not doc.file:
        raise Http404("Document has no file")

    return _serve_file(doc.file, doc.title + Path(doc.file.name).suffix)


@api_view(['GET'])
@permission_classes([IsAuthenticated])
def serve_batch_qc_file(request, qc_file_id):
    """
    Serve a batch QC file.

    URL: /api/compounds/media/batch-qc-files/<qc_file_id>/file/
    """
    try:
        qc_file = BatchQCFile.objects.get(id=qc_file_id)
    except BatchQCFile.DoesNotExist:
        raise Http404("QC file not found")

    if not qc_file.file:
        raise Http404("QC record has no file")

    return _serve_file(qc_file.file)


@api_view(['GET'])
@permission_classes([IsAuthenticated])
def serve_plasmid_genbank(request, plasmid_id):
    """
    Serve a plasmid's GenBank file.

    URL: /api/compounds/media/plasmids/<plasmid_id>/genbank/
    """
    try:
        plasmid = Plasmid.objects.get(id=plasmid_id)
    except Plasmid.DoesNotExist:
        raise Http404("Plasmid not found")

    if not plasmid.genbank_file:
        raise Http404("Plasmid has no GenBank file")

    return _serve_file(plasmid.genbank_file)


@api_view(['GET'])
@permission_classes([IsAuthenticated])
def serve_cassette_use_alignment(request, cassette_use_id):
    """
    Serve a cassette use's alignment file.

    URL: /api/compounds/media/cassette-uses/<cassette_use_id>/alignment/
    """
    try:
        cassette_use = CassetteUse.objects.get(id=cassette_use_id)
    except CassetteUse.DoesNotExist:
        raise Http404("Cassette use not found")

    if not cassette_use.alignment_file:
        raise Http404("Cassette use has no alignment file")

    return _serve_file(cassette_use.alignment_file)


@api_view(['GET'])
@permission_classes([IsAuthenticated])
def serve_sequencing_result(request, result_id):
    """
    Serve a sequencing result file.

    URL: /api/compounds/media/sequencing-results/<result_id>/file/
    """
    try:
        result = SequencingResult.objects.get(id=result_id)
    except SequencingResult.DoesNotExist:
        raise Http404("Sequencing result not found")

    if not result.file:
        raise Http404("Sequencing result has no file")

    return _serve_file(result.file)
