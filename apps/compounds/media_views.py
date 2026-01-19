"""
Protected media file serving for compounds app.

These views serve files from Django FileFields through authenticated endpoints,
ensuring that sensitive data (assay files, QC documents, etc.) are not publicly accessible.

Authentication is handled by the AzureADAuthMiddleware - these views only need to
ensure they're not added to exempt paths.

When Azure Blob Storage is configured (via django-storages), files are served using
time-limited SAS URLs for direct download from blob storage. This provides:
- Better performance (clients download directly from Azure CDN)
- Reduced server load (no file streaming through Django)
- Security via short-lived tokens (default 1 hour)
"""

import logging
import os
import mimetypes
from pathlib import Path

from django.conf import settings
from django.http import FileResponse, Http404, HttpResponseBadRequest, HttpResponseRedirect
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated

from compounds.assays.models import Assay, ProtocolDocument
from compounds.registry.models import BatchQCFile
from compounds.constructs.models import Plasmid, CassetteUse, SequencingResult

logger = logging.getLogger(__name__)


def _is_azure_storage_configured():
    """Check if Azure Blob Storage is configured for file storage."""
    try:
        from azure_extensions.utils.blob_sas import is_blob_storage_configured
        return is_blob_storage_configured()
    except ImportError:
        return False


def _generate_sas_download_url(blob_path: str, filename: str) -> str:
    """
    Generate a SAS URL for downloading a file from Azure Blob Storage.

    Args:
        blob_path: Path to the blob (from FileField.name)
        filename: Filename for Content-Disposition header

    Returns:
        SAS URL with time-limited read access

    Raises:
        ValueError: If blob storage is not configured or file doesn't exist
    """
    from azure_extensions.utils.blob_sas import generate_download_sas_url
    sas_url, expiry = generate_download_sas_url(blob_path, filename)
    return sas_url


def _serve_file(file_field, filename_override=None):
    """
    Serve a file from a FileField with proper headers.

    When Azure Blob Storage is configured, returns a redirect to a time-limited
    SAS URL for direct download from blob storage.

    When using local filesystem storage, streams the file through Django.

    Args:
        file_field: Django FileField instance
        filename_override: Optional filename for Content-Disposition header

    Returns:
        HttpResponseRedirect to SAS URL (Azure) or FileResponse (local filesystem)
    """
    if not file_field:
        raise Http404("File not found")

    # Get filename for Content-Disposition
    filename = filename_override or os.path.basename(file_field.name)

    # Check if using Azure Blob Storage
    if _is_azure_storage_configured():
        sas_error = None
        try:
            # Generate SAS URL for direct blob download
            # The blob path is stored in file_field.name (relative path in container)
            blob_path = file_field.name
            sas_url = _generate_sas_download_url(blob_path, filename)
            logger.info(f"[SAS] Redirecting to SAS URL for: {blob_path}")
            response = HttpResponseRedirect(sas_url)
            response['X-Download-Method'] = 'sas-redirect'
            return response
        except ValueError as e:
            sas_error = str(e)
            logger.error(f"[SAS] Failed to generate SAS URL: {e}")
        except Exception as e:
            sas_error = str(e)
            logger.exception(f"[SAS] Error generating SAS URL for {file_field.name}")

        # Stream file through Django when SAS generation fails
        # This uses the server's Managed Identity to read from Azure
        logger.info(f"[STREAM] Falling back to streaming file through Django: {file_field.name}")
        try:
            content_type, _ = mimetypes.guess_type(filename)
            if not content_type:
                content_type = 'application/octet-stream'

            # Read from Azure storage via django-storages and stream to client
            response = FileResponse(
                file_field.open('rb'),
                content_type=content_type,
            )
            viewable_types = [
                'application/pdf', 'text/plain', 'text/csv',
                'image/png', 'image/jpeg', 'image/gif',
            ]
            disposition = 'inline' if content_type in viewable_types else 'attachment'
            response['Content-Disposition'] = f'{disposition}; filename="{filename}"'
            response['X-Download-Method'] = 'django-stream'
            if sas_error:
                response['X-SAS-Error'] = sas_error[:200]  # Truncate for header safety
            return response
        except Exception as e:
            logger.exception(f"[STREAM] Failed to stream file from Azure: {e}")
            raise Http404(f"File not accessible: {e}")

    # Local filesystem access (non-Azure deployment)
    try:
        file_path = file_field.path
    except NotImplementedError:
        # django-storages backends may not support .path - but we should have handled Azure above
        raise Http404("File storage backend does not support direct file access")

    if not os.path.exists(file_path):
        raise Http404("File not found on disk")

    # Determine content type
    content_type, _ = mimetypes.guess_type(file_path)
    if not content_type:
        content_type = 'application/octet-stream'

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

    # Use the original filename from the stored path
    return _serve_file(doc.file)


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

    # Use the original filename for Content-Disposition
    return _serve_file(qc_file.file, qc_file.filename)


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
