"""Reap expired and consumed staged uploads.

The deployment schedules this (e.g. from the existing maintenance job): it deletes
staging rows and directories older than ``CCP4I2_IMPORT_STAGING_TTL_HOURS`` and
any consumed row whose file is already gone. Idempotent, and a no-op when
``CCP4I2_IMPORT_STAGING_DIR`` is unset.

    ccp4-python manage.py sweep_staged_uploads
"""

from django.core.management.base import BaseCommand

from ccp4i2.lib.utils.files import staged_upload


class Command(BaseCommand):
    help = "Delete expired/consumed staged uploads and their directories."

    def handle(self, *args, **options):
        if not staged_upload.staging_enabled():
            self.stdout.write("Staged upload is not enabled (no staging dir); nothing to do.")
            return
        removed = staged_upload.sweep()
        self.stdout.write(f"Swept {removed} staged upload(s).")
