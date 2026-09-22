"""
Django management command to clear non-finite KPI values out of the database.

A KPI that came out NaN or infinite is a failed measurement, not a measurement.
Nothing should write one any more (see `ccp4i2.lib.kpi_values`), but rows
written before that gate existed are still there, and PostgreSQL stores them
quite happily. While such a row exists, every endpoint serving that job leans
on the renderer backstop, and its KPI is silently missing from the UI.

Reports by default and changes nothing without --apply.

Usage:
    python manage.py prune_nonfinite_kpis
    python manage.py prune_nonfinite_kpis --project toxd
    python manage.py prune_nonfinite_kpis --apply
"""

import math

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from ccp4i2.db.models import JobFloatValue, Project


class Command(BaseCommand):
    """Find, report and optionally delete non-finite JobFloatValue rows."""

    help = "Report (or with --apply, delete) KPI rows holding NaN or infinity"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            help="Limit to one project, by name (default: every project)",
            type=str,
        )
        parser.add_argument(
            "--apply",
            action="store_true",
            help="Delete the rows found. Without this, only reports.",
        )

    def handle(self, *args, **options):
        rows = JobFloatValue.objects.select_related("job", "job__project", "key")

        project_name = options.get("project")
        if project_name:
            try:
                project = Project.objects.get(name=project_name)
            except Project.DoesNotExist:
                raise CommandError(f"No project named {project_name}")
            rows = rows.filter(job__project=project)

        # No portable ORM predicate for NaN (and SQLite cannot even hold one),
        # so the scan happens in Python. .iterator() keeps a large database off
        # the heap; these tables are one row per KPI per job.
        offenders = [
            row for row in rows.iterator()
            if not math.isfinite(row.value)
        ]

        if not offenders:
            self.stdout.write(self.style.SUCCESS("No non-finite KPI values found"))
            return

        self.stdout.write(
            f"Found {len(offenders)} non-finite KPI value(s):"
        )
        for row in offenders:
            self.stdout.write(
                f"  {row.job.project.name} job {row.job.number} "
                f"{row.key_id} = {row.value}"
            )

        if not options.get("apply"):
            self.stdout.write(
                self.style.WARNING(
                    "Nothing changed. Re-run with --apply to delete these rows."
                )
            )
            return

        with transaction.atomic():
            deleted, _ = JobFloatValue.objects.filter(
                pk__in=[row.pk for row in offenders]
            ).delete()

        self.stdout.write(self.style.SUCCESS(f"Deleted {deleted} row(s)"))
