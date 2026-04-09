"""
Management command to clean up stale jobs stuck in RUNNING or RUNNING_REMOTELY state.

This handles cases where:
- Worker was OOM-killed without catching SIGTERM signal
- Container crashed unexpectedly
- Network partition prevented status update
- Worker process was killed during job execution

Usage:
    python manage.py cleanup_stale_jobs --hours 2

This is typically called at worker startup to clean up any orphaned jobs.
"""
from datetime import timedelta
from django.core.management.base import BaseCommand
from django.utils import timezone
from ccp4i2.db.models import Job


class Command(BaseCommand):
    help = "Mark stale RUNNING/RUNNING_REMOTELY jobs as FAILED"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            "--hours",
            type=float,
            default=2.0,
            help="Consider jobs stale if running longer than this many hours (default: 2)",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Show what would be cleaned up without making changes",
        )

    def handle(self, *args, **options):
        threshold_hours = options["hours"]
        dry_run = options["dry_run"]

        # Calculate the cutoff time
        cutoff_time = timezone.now() - timedelta(hours=threshold_hours)

        # Find jobs that are stuck in running state for too long
        # These jobs were likely being processed when a worker crashed
        stale_jobs = Job.objects.filter(
            status__in=[Job.Status.RUNNING, Job.Status.RUNNING_REMOTELY],
            creation_time__lt=cutoff_time,
        )

        count = stale_jobs.count()

        if count == 0:
            self.stdout.write("No stale jobs found")
            return

        if dry_run:
            self.stdout.write(
                self.style.WARNING(f"[DRY RUN] Would mark {count} stale jobs as FAILED:")
            )
            for job in stale_jobs:
                age_hours = (timezone.now() - job.creation_time).total_seconds() / 3600
                self.stdout.write(
                    f"  - Job {job.id} ({job.uuid}): {job.title} "
                    f"[status={job.get_status_display()}, age={age_hours:.1f}h]"
                )
            return

        # Mark all stale jobs as FAILED
        updated_count = 0
        for job in stale_jobs:
            age_hours = (timezone.now() - job.creation_time).total_seconds() / 3600
            old_status = job.get_status_display()

            job.status = Job.Status.FAILED
            job.save()

            self.stdout.write(
                self.style.WARNING(
                    f"Marked job {job.id} ({job.uuid}) as FAILED "
                    f"[was {old_status}, age={age_hours:.1f}h]: {job.title}"
                )
            )
            updated_count += 1

        self.stdout.write(
            self.style.SUCCESS(f"Cleaned up {updated_count} stale jobs")
        )
