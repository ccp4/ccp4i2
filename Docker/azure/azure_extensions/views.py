"""
ViewSet for managing staged uploads (large file uploads via SAS URL).

This enables uploading large files (up to 20GB) that would otherwise fail
due to HTTP body size limits on Azure Container Apps.

Flow:
1. POST /uploads/request - Get a SAS URL for uploading
2. Client uploads directly to Azure Blob Storage using azcopy or curl
3. POST /uploads/{id}/complete - Trigger processing of the uploaded file
"""

import logging

from django.utils import timezone
from rest_framework import status
from rest_framework.decorators import action
from rest_framework.parsers import JSONParser
from rest_framework.response import Response
from rest_framework.viewsets import ModelViewSet

from .models import StagedUpload

logger = logging.getLogger(f"azure_extensions:{__name__}")


def api_success(data, status_code=200):
    """Return a success response."""
    return Response(data, status=status_code)


def api_error(message, status=400):
    """Return an error response."""
    return Response({"error": message}, status=status)


class StagedUploadSerializer:
    """Simple serializer for StagedUpload model."""

    @staticmethod
    def to_dict(upload: StagedUpload) -> dict:
        return {
            "id": upload.id,
            "uuid": str(upload.uuid),
            "upload_type": upload.upload_type,
            "status": upload.status,
            "original_filename": upload.original_filename,
            "blob_path": upload.blob_path,
            "sas_expiry": upload.sas_expiry.isoformat() if upload.sas_expiry else None,
            "created_at": upload.created_at.isoformat() if upload.created_at else None,
            "completed_at": upload.completed_at.isoformat() if upload.completed_at else None,
            "error_message": upload.error_message,
            "is_expired": upload.is_expired,
        }


class StagedUploadViewSet(ModelViewSet):
    """
    ViewSet for staged uploads.

    Provides endpoints for:
    - Requesting an upload URL (POST /uploads/request/)
    - Completing an upload (POST /uploads/{id}/complete/)
    - Listing uploads (GET /uploads/)
    - Getting upload status (GET /uploads/{id}/)
    """

    queryset = StagedUpload.objects.all()
    parser_classes = [JSONParser]
    lookup_field = "uuid"

    def list(self, request):
        """List all staged uploads (optionally filtered by status)."""
        queryset = self.queryset

        status_filter = request.query_params.get("status")
        if status_filter:
            queryset = queryset.filter(status=status_filter)

        upload_type = request.query_params.get("type")
        if upload_type:
            queryset = queryset.filter(upload_type=upload_type)

        uploads = [StagedUploadSerializer.to_dict(u) for u in queryset.order_by("-created_at")[:50]]
        return Response(uploads)

    def retrieve(self, request, uuid=None):
        """Get details of a specific upload."""
        try:
            upload = StagedUpload.objects.get(uuid=uuid)
            return Response(StagedUploadSerializer.to_dict(upload))
        except StagedUpload.DoesNotExist:
            return api_error("Upload not found", status=404)

    @action(detail=False, methods=["post"], url_path="request")
    def request_upload(self, request):
        """
        Request a SAS URL for uploading a large file.

        Request body:
        {
            "filename": "myproject.ccp4_project.zip",
            "upload_type": "project_import" | "unmerged_data",
            "target_job_uuid": "..." (optional, for unmerged_data)
        }

        Response:
        {
            "upload_id": "uuid",
            "sas_url": "https://...",
            "blob_path": "project_import/uuid.zip",
            "expiry": "2024-01-01T12:00:00Z",
            "instructions": "..."
        }
        """
        try:
            from .utils.blob_sas import generate_upload_sas_url
        except ImportError as e:
            logger.error(f"Storage utilities not available: {e}")
            return api_error(
                "Staged uploads not configured. Azure storage libraries required.",
                status=501
            )

        # Parse request
        filename = request.data.get("filename")
        upload_type = request.data.get("upload_type")
        target_job_uuid = request.data.get("target_job_uuid")

        if not filename:
            return api_error("filename is required", status=400)

        if upload_type not in [StagedUpload.UploadType.PROJECT_IMPORT, StagedUpload.UploadType.UNMERGED_DATA]:
            return api_error(
                f"upload_type must be one of: {', '.join(StagedUpload.UploadType.values)}",
                status=400
            )

        # Validate file extension
        if upload_type == StagedUpload.UploadType.PROJECT_IMPORT:
            if not filename.endswith(".zip"):
                return api_error("Project import files must be .zip archives", status=400)

        # Get user info if available
        requested_by = ""
        if hasattr(request, "user") and request.user:
            requested_by = str(request.user)

        try:
            # Generate SAS URL
            sas_url, blob_path, expiry = generate_upload_sas_url(
                filename=filename,
                upload_type=upload_type,
                expiry_hours=2,  # 2 hours for large uploads
            )

            # Create database record
            upload = StagedUpload.objects.create(
                upload_type=upload_type,
                original_filename=filename,
                blob_path=blob_path,
                sas_expiry=expiry,
                requested_by=requested_by,
                status=StagedUpload.Status.PENDING,
            )

            # Link to target job if provided
            if target_job_uuid:
                from ccp4i2.db.models import Job
                try:
                    job = Job.objects.get(uuid=target_job_uuid)
                    upload.target_job = job
                    upload.save()
                except Job.DoesNotExist:
                    logger.warning(f"Target job not found: {target_job_uuid}")

            logger.info(f"Created staged upload {upload.uuid} for {filename}")

            return api_success({
                "upload_id": str(upload.uuid),
                "sas_url": sas_url,
                "blob_path": blob_path,
                "expiry": expiry.isoformat(),
                "instructions": (
                    f"Upload your file using:\n"
                    f"  azcopy copy '{filename}' '{sas_url}'\n"
                    f"Or:\n"
                    f"  curl -X PUT -H 'x-ms-blob-type: BlockBlob' --data-binary @'{filename}' '{sas_url}'\n"
                    f"\nThen complete the upload with:\n"
                    f"  i2remote upload-complete {upload.uuid}"
                ),
            })

        except ValueError as e:
            logger.error(f"SAS generation failed: {e}")
            return api_error(str(e), status=500)
        except Exception as e:
            logger.exception("Unexpected error generating upload URL")
            return api_error(f"Failed to generate upload URL: {e}", status=500)

    @action(detail=True, methods=["post"], url_path="complete")
    def complete_upload(self, request, uuid=None):
        """
        Mark an upload as complete and trigger processing.

        This should be called after the client has finished uploading
        the file to the SAS URL.
        """
        try:
            upload = StagedUpload.objects.get(uuid=uuid)
        except StagedUpload.DoesNotExist:
            return api_error("Upload not found", status=404)

        # Check if already processed
        if upload.status in [StagedUpload.Status.PROCESSING, StagedUpload.Status.COMPLETED]:
            return api_error(f"Upload already {upload.status}", status=400)

        # Check if expired
        if upload.is_expired:
            upload.status = StagedUpload.Status.EXPIRED
            upload.save()
            return api_error("Upload URL has expired. Please request a new one.", status=410)

        # Verify blob exists
        try:
            from .utils.blob_sas import verify_blob_exists, get_blob_local_path
        except ImportError:
            return api_error("Storage utilities not available", status=501)

        if not verify_blob_exists(upload.blob_path):
            return api_error(
                "File not found in storage. Please ensure upload completed successfully.",
                status=404
            )

        # Mark as uploaded
        upload.status = StagedUpload.Status.UPLOADED
        upload.save()

        # Process based on upload type
        try:
            if upload.upload_type == StagedUpload.UploadType.PROJECT_IMPORT:
                result = self._process_project_import(upload)
            elif upload.upload_type == StagedUpload.UploadType.UNMERGED_DATA:
                result = self._process_unmerged_data(upload)
            else:
                return api_error(f"Unknown upload type: {upload.upload_type}", status=400)

            return api_success(result)

        except Exception as e:
            upload.status = StagedUpload.Status.FAILED
            upload.error_message = str(e)
            upload.save()
            logger.exception(f"Failed to process upload {upload.uuid}")
            return api_error(f"Processing failed: {e}", status=500)

    def _process_project_import(self, upload: StagedUpload) -> dict:
        """Process a project import upload."""
        from .utils.blob_sas import get_blob_local_path, delete_blob
        from django.core.management import call_command

        upload.status = StagedUpload.Status.PROCESSING
        upload.save()

        try:
            # Download blob to local temp file
            local_path = get_blob_local_path(upload.blob_path)

            logger.info(f"Starting project import from {local_path}")

            # Run import command (detached so it doesn't block)
            call_command("import_ccp4_project_zip", local_path, "--detach")

            upload.status = StagedUpload.Status.COMPLETED
            upload.completed_at = timezone.now()
            upload.save()

            # Clean up blob (async, don't block response)
            try:
                delete_blob(upload.blob_path)
            except Exception as e:
                logger.warning(f"Failed to delete staging blob: {e}")

            return {
                "message": "Project import started",
                "upload_id": str(upload.uuid),
                "status": upload.status,
                "note": "Import runs in background. Project will appear shortly.",
            }

        except Exception as e:
            upload.status = StagedUpload.Status.FAILED
            upload.error_message = str(e)
            upload.save()
            raise

    def _process_unmerged_data(self, upload: StagedUpload) -> dict:
        """Process an unmerged data upload."""
        from .utils.blob_sas import get_blob_local_path

        upload.status = StagedUpload.Status.PROCESSING
        upload.save()

        try:
            # Download blob to local temp file
            local_path = get_blob_local_path(upload.blob_path)

            # For unmerged data, we need to move it to the appropriate job directory
            if not upload.target_job:
                raise ValueError("No target job specified for unmerged data upload")

            # TODO: Implement unmerged data handling
            # This would involve:
            # 1. Moving the file to the job's input directory
            # 2. Setting the appropriate job parameter
            # 3. Optionally triggering the job

            upload.status = StagedUpload.Status.COMPLETED
            upload.completed_at = timezone.now()
            upload.save()

            return {
                "message": "Unmerged data uploaded",
                "upload_id": str(upload.uuid),
                "local_path": local_path,
                "status": upload.status,
            }

        except Exception as e:
            upload.status = StagedUpload.Status.FAILED
            upload.error_message = str(e)
            upload.save()
            raise

    @action(detail=True, methods=["delete"], url_path="cancel")
    def cancel_upload(self, request, uuid=None):
        """Cancel a pending upload and clean up resources."""
        try:
            upload = StagedUpload.objects.get(uuid=uuid)
        except StagedUpload.DoesNotExist:
            return api_error("Upload not found", status=404)

        if upload.status not in [StagedUpload.Status.PENDING, StagedUpload.Status.UPLOADED]:
            return api_error(f"Cannot cancel upload in {upload.status} status", status=400)

        # Try to delete the blob if it exists
        try:
            from .utils.blob_sas import delete_blob
            delete_blob(upload.blob_path)
        except Exception as e:
            logger.warning(f"Failed to delete blob during cancel: {e}")

        # Delete the database record
        upload.delete()

        return api_success({"message": "Upload cancelled"})
