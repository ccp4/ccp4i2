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
        Mark an upload as complete and queue it for background processing.

        This should be called after the client has finished uploading
        the file to the SAS URL. The actual processing (download from blob,
        import/processing) happens asynchronously in the worker.
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
            from .utils.blob_sas import verify_blob_exists
        except ImportError:
            return api_error("Storage utilities not available", status=501)

        if not verify_blob_exists(upload.blob_path):
            return api_error(
                "File not found in storage. Please ensure upload completed successfully.",
                status=404
            )

        # Mark as uploaded (verified blob exists)
        upload.status = StagedUpload.Status.UPLOADED
        upload.save()

        # Queue for background processing based on upload type
        try:
            if upload.upload_type == StagedUpload.UploadType.PROJECT_IMPORT:
                result = self._queue_project_import(upload)
            elif upload.upload_type == StagedUpload.UploadType.UNMERGED_DATA:
                result = self._queue_unmerged_data(upload)
            else:
                return api_error(f"Unknown upload type: {upload.upload_type}", status=400)

            return api_success(result)

        except Exception as e:
            upload.status = StagedUpload.Status.FAILED
            upload.error_message = str(e)
            upload.save()
            logger.exception(f"Failed to queue upload {upload.uuid}")
            return api_error(f"Failed to queue processing: {e}", status=500)

    def _queue_project_import(self, upload: StagedUpload) -> dict:
        """Queue a project import for background processing."""
        from .utils.queue import send_to_queue

        # Prepare message for worker
        message = {
            "action": "import_project",
            "upload_id": str(upload.uuid),
            "blob_path": upload.blob_path,
            "original_filename": upload.original_filename,
        }

        # Send to queue
        if not send_to_queue(message):
            raise RuntimeError("Failed to send message to processing queue")

        # Mark as processing (queued for worker)
        upload.status = StagedUpload.Status.PROCESSING
        upload.save()

        logger.info(f"Queued project import for upload {upload.uuid}")

        return {
            "message": "Project import queued for processing",
            "upload_id": str(upload.uuid),
            "status": upload.status,
            "note": "Import runs in background. Project will appear when complete.",
        }

    def _queue_unmerged_data(self, upload: StagedUpload) -> dict:
        """Queue unmerged data upload for background processing."""
        from .utils.queue import send_to_queue

        if not upload.target_job:
            raise ValueError("No target job specified for unmerged data upload")

        # Prepare message for worker
        message = {
            "action": "process_unmerged_data",
            "upload_id": str(upload.uuid),
            "blob_path": upload.blob_path,
            "original_filename": upload.original_filename,
            "target_job_id": upload.target_job.id,
            "target_job_uuid": str(upload.target_job.uuid),
        }

        # Send to queue
        if not send_to_queue(message):
            raise RuntimeError("Failed to send message to processing queue")

        # Mark as processing (queued for worker)
        upload.status = StagedUpload.Status.PROCESSING
        upload.save()

        logger.info(f"Queued unmerged data processing for upload {upload.uuid}")

        return {
            "message": "Unmerged data processing queued",
            "upload_id": str(upload.uuid),
            "status": upload.status,
            "target_job": str(upload.target_job.uuid),
        }

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

    @action(detail=True, methods=["post"], url_path="reset")
    def reset_upload(self, request, uuid=None):
        """
        Reset a stuck upload back to 'uploaded' status.

        This allows re-triggering processing for uploads that got stuck
        in 'processing' state due to timeouts or errors.
        """
        logger.info(f"reset_upload called for uuid={uuid}")

        try:
            upload = StagedUpload.objects.get(uuid=uuid)
            logger.info(f"Found upload: status={upload.status}, blob_path={upload.blob_path}")
        except StagedUpload.DoesNotExist:
            logger.warning(f"Upload not found: {uuid}")
            return api_error("Upload not found", status=404)
        except Exception as e:
            logger.exception(f"Error fetching upload {uuid}: {e}")
            return api_error(f"Database error: {e}", status=500)

        if upload.status not in [StagedUpload.Status.PROCESSING, StagedUpload.Status.FAILED]:
            logger.info(f"Upload {uuid} in wrong status: {upload.status}")
            return api_error(
                f"Cannot reset upload in {upload.status} status. "
                "Only 'processing' or 'failed' uploads can be reset.",
                status=400
            )

        # Verify the blob still exists before resetting (skip on any error)
        try:
            from .utils.blob_sas import verify_blob_exists
            blob_exists = verify_blob_exists(upload.blob_path)
            logger.info(f"Blob verification: exists={blob_exists}")
            if not blob_exists:
                return api_error(
                    "Blob no longer exists. Cannot reset this upload.",
                    status=410
                )
        except ImportError as e:
            logger.info(f"Skipping blob verification (ImportError): {e}")
        except Exception as e:
            # Log warning but allow reset - blob check may fail due to permissions
            logger.warning(f"Could not verify blob exists for reset: {e}")

        # Reset to uploaded status
        try:
            upload.status = StagedUpload.Status.UPLOADED
            upload.error_message = None
            upload.save()
            logger.info(f"Reset upload {upload.uuid} to 'uploaded' status")
        except Exception as e:
            logger.exception(f"Error saving upload {uuid}: {e}")
            return api_error(f"Failed to save upload: {e}", status=500)

        return api_success({
            "message": "Upload reset to 'uploaded' status",
            "upload_id": str(upload.uuid),
            "status": upload.status,
        })

    @action(detail=True, methods=["post"], url_path="force-complete")
    def force_complete_upload(self, request, uuid=None):
        """
        Force completion of an upload, bypassing expiry check.

        Use this when the blob exists but the SAS URL has expired.
        The blob will be verified before queuing for processing.
        """
        logger.info(f"force_complete_upload called for uuid={uuid}")

        try:
            upload = StagedUpload.objects.get(uuid=uuid)
        except StagedUpload.DoesNotExist:
            return api_error("Upload not found", status=404)

        # Check if already processed
        if upload.status in [StagedUpload.Status.PROCESSING, StagedUpload.Status.COMPLETED]:
            return api_error(f"Upload already {upload.status}", status=400)

        # Verify blob exists (this is the key check - blob must exist)
        try:
            from .utils.blob_sas import verify_blob_exists
        except ImportError:
            return api_error("Storage utilities not available", status=501)

        if not verify_blob_exists(upload.blob_path):
            return api_error(
                "File not found in storage. The blob may have been deleted.",
                status=404
            )

        # Mark as uploaded (verified blob exists)
        upload.status = StagedUpload.Status.UPLOADED
        upload.save()

        # Queue for background processing based on upload type
        try:
            if upload.upload_type == StagedUpload.UploadType.PROJECT_IMPORT:
                result = self._queue_project_import(upload)
            elif upload.upload_type == StagedUpload.UploadType.UNMERGED_DATA:
                result = self._queue_unmerged_data(upload)
            else:
                return api_error(f"Unknown upload type: {upload.upload_type}", status=400)

            return api_success(result)

        except Exception as e:
            upload.status = StagedUpload.Status.FAILED
            upload.error_message = str(e)
            upload.save()
            logger.exception(f"Failed to queue upload {upload.uuid}")
            return api_error(f"Failed to queue processing: {e}", status=500)
