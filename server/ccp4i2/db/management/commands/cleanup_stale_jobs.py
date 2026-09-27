"""
Management command to clean up jobs a worker was running when it died.

This handles cases where:
- Worker was OOM-killed without catching SIGTERM signal
- Container crashed unexpectedly
- Network partition prevented status update
- Worker process was killed during job execution

Usage:
    python manage.py cleanup_stale_jobs --hours 2

This is typically called at worker startup to clean up any orphaned jobs.

RUNNING_REMOTELY is deliberately NOT swept by age. That status means no
worker holds the job: its program runs on a program run target and the job
waits for a reconcile, so a worker dying tells us nothing about it and age
tells us nothing either. A PanDDA campaign routinely runs for longer than any
threshold worth setting. This command used to include it, and on 2026-09-27 a
deployment restarted the worker 2.7 hours into a healthy Azure Batch run: the
job was marked FAILED while its task kept running, and nothing could harvest
it afterwards because the reconcile only acts on RUNNING_REMOTELY.

Whether such a run is alive is knowable exactly rather than by a timer, by
asking the target -- which is what `reconcile_dispatch --all` does, and what a
worker should run at startup beside this command. The one dispatched job this
command can still judge is one with no dispatch record at all: nothing can
ever reconcile it, so nothing will ever finish it.
"""
from datetime import timedelta
from django.core.management.base import BaseCommand
from django.utils import timezone
from ccp4i2.db.models import Job


class Command(BaseCommand):
    help = "Mark jobs FAILED that a dead worker left RUNNING, and dispatched jobs with no record"
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

    @staticmethod
    def _has_dispatch_record(job) -> bool:
        """Whether anything could still reconcile this job.

        Unreadable is treated as present: a record that cannot be parsed right
        now (a partial write, a filesystem blip) is a reason to leave the job
        alone and let the reconcile report it, not to declare the run dead.
        """
        from ccp4i2.lib.utils.jobs.dispatch_record import record_path

        try:
            return record_path(job.directory).is_file()
        except (OSError, TypeError, ValueError):
            return True

    def handle(self, *args, **options):
        threshold_hours = options["hours"]
        dry_run = options["dry_run"]

        # Calculate the cutoff time
        cutoff_time = timezone.now() - timedelta(hours=threshold_hours)

        # Jobs a worker was executing when it died. An interactive job (a
        # recorded Moorhen session) is RUNNING for as long as its window is
        # open, which is routinely longer than the threshold; while its
        # session is open it is not stale.
        stale_jobs = list(
            Job.objects.filter(
                status=Job.Status.RUNNING,
                creation_time__lt=cutoff_time,
            ).exclude(interactive_session__finished=False)
        )

        # Dispatched jobs are judged by whether anything CAN finish them, not
        # by age (see the module docstring). A job in RUNNING_REMOTELY with no
        # dispatch record has no handle to poll and no target to ask, so no
        # reconcile will ever move it; it is stuck in the only sense this
        # command can establish. One with a record is left alone however old
        # it is, and `reconcile_dispatch --all` is what speaks for it.
        stale_jobs += [
            job
            for job in Job.objects.filter(
                status=Job.Status.RUNNING_REMOTELY,
                creation_time__lt=cutoff_time,
            )
            if not self._has_dispatch_record(job)
        ]

        count = len(stale_jobs)

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
