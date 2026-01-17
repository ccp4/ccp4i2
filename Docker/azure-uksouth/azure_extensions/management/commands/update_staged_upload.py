"""
Management command to update staged upload status.

Used by the worker to update upload status after processing.
"""

from django.core.management.base import BaseCommand
from django.utils import timezone

from azure_extensions.models import StagedUpload


class Command(BaseCommand):
    help = "Update the status of a staged upload"

    def add_arguments(self, parser):
        parser.add_argument(
            "--upload-id",
            required=True,
            help="UUID of the staged upload",
        )
        parser.add_argument(
            "--status",
            required=True,
            choices=["completed", "failed", "uploaded"],
            help="New status for the upload",
        )
        parser.add_argument(
            "--error",
            default=None,
            help="Error message (for failed status)",
        )

    def handle(self, *args, **options):
        upload_id = options["upload_id"]
        status = options["status"]
        error_message = options.get("error")

        try:
            upload = StagedUpload.objects.get(uuid=upload_id)
        except StagedUpload.DoesNotExist:
            self.stderr.write(f"Upload not found: {upload_id}")
            return

        if status == "completed":
            upload.status = StagedUpload.Status.COMPLETED
            upload.completed_at = timezone.now()
            upload.error_message = None
        elif status == "failed":
            upload.status = StagedUpload.Status.FAILED
            upload.error_message = error_message
        elif status == "uploaded":
            upload.status = StagedUpload.Status.UPLOADED
            upload.error_message = None

        upload.save()

        self.stdout.write(
            self.style.SUCCESS(f"Updated upload {upload_id} to status {status}")
        )
