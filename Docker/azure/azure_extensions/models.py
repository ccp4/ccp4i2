"""
Azure-specific Django models.

These models are only used in Azure deployments.
"""

from uuid import uuid4

from django.db.models import (
    Model,
    CharField,
    DateTimeField,
    ForeignKey,
    TextField,
    TextChoices,
    UUIDField,
    SET_NULL,
)
from django.utils import timezone


class StagedUpload(Model):
    """
    Track staged uploads that use SAS URLs for large file transfers.

    This enables uploading files larger than the HTTP body size limit
    by generating time-limited SAS URLs for direct blob storage uploads.
    """

    class Status(TextChoices):
        PENDING = "pending", "Pending upload"
        UPLOADED = "uploaded", "Uploaded, awaiting processing"
        PROCESSING = "processing", "Processing"
        COMPLETED = "completed", "Completed"
        FAILED = "failed", "Failed"
        EXPIRED = "expired", "Expired"

    class UploadType(TextChoices):
        PROJECT_IMPORT = "project_import", "Project Import"
        UNMERGED_DATA = "unmerged_data", "Unmerged Data"

    uuid = UUIDField(default=uuid4, unique=True, editable=False)
    upload_type = CharField(max_length=32, choices=UploadType.choices)
    status = CharField(max_length=32, choices=Status.choices, default=Status.PENDING)
    original_filename = CharField(max_length=255)
    blob_path = CharField(max_length=512)
    sas_expiry = DateTimeField()
    # Reference to target job (for unmerged data uploads)
    target_job = ForeignKey(
        "ccp4i2.Job",
        on_delete=SET_NULL,
        blank=True,
        null=True,
        related_name="staged_uploads",
    )
    created_at = DateTimeField(default=timezone.now)
    completed_at = DateTimeField(blank=True, null=True)
    error_message = TextField(blank=True, null=True)
    requested_by = CharField(max_length=255, blank=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self):
        return f"StagedUpload {self.uuid} ({self.upload_type}: {self.status})"

    @property
    def is_expired(self) -> bool:
        """Check if the SAS URL has expired."""
        if not self.sas_expiry:
            return True
        return timezone.now() > self.sas_expiry
