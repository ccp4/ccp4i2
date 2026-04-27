"""Delete expired Selection rows.

Selections created via the NLP redirect path auto-expire after 7 days
(see ``DEFAULT_EXPIRY_DAYS`` in ``selection_views``). This command
sweeps the expired rows. Saved selections (``is_saved=True``) are
never touched.

Run via cron / Azure scheduled job:
    ccp4-python manage.py cleanup_selections
"""

from __future__ import annotations

from django.core.management.base import BaseCommand
from django.utils import timezone

from compounds.registry.models import Selection


class Command(BaseCommand):
    help = "Delete expired Selection rows."

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="List what would be deleted without actually deleting.",
        )

    def handle(self, *args, dry_run: bool = False, **options):
        now = timezone.now()
        expired_qs = Selection.objects.filter(
            is_saved=False, expires_at__lt=now,
        )
        count = expired_qs.count()

        if count == 0:
            self.stdout.write("No expired selections.")
            return

        if dry_run:
            self.stdout.write(
                f"DRY RUN: would delete {count} expired selection(s):"
            )
            for sel in expired_qs[:20]:
                self.stdout.write(
                    f"  {sel.id}  {sel.name!r:50}  "
                    f"expired {sel.expires_at.isoformat()}"
                )
            if count > 20:
                self.stdout.write(f"  ... and {count - 20} more")
            return

        deleted, _ = expired_qs.delete()
        self.stdout.write(
            self.style.SUCCESS(f"Deleted {deleted} expired Selection(s).")
        )
