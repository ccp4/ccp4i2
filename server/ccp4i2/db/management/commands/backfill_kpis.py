"""
Django management command to backfill KPIs from params.xml for finished jobs.

Parses the PERFORMANCE or PERFORMANCEINDICATOR element from each job's
params.xml and registers any numeric/string values as JobFloatValue/JobCharValue.

Usage:
    # Dry run (show what would be backfilled, no DB writes)
    python manage.py backfill_kpis

    # Actually write to the database
    python manage.py backfill_kpis --commit

    # Limit to a specific project
    python manage.py backfill_kpis --commit --projectname toxd

    # Only backfill jobs that currently have no KPIs
    python manage.py backfill_kpis --commit --missing-only

    # Limit to specific task names
    python manage.py backfill_kpis --commit --tasks aimless_pipe xia2_dials
"""

from lxml import etree
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from ccp4i2.db.models import Job, Project, JobValueKey, JobFloatValue, JobCharValue


# KPI element names to look for in outputData, in priority order
KPI_ELEMENT_NAMES = ('PERFORMANCEINDICATOR', 'PERFORMANCE')

# Elements that are not KPI values (structural/container elements)
SKIP_ELEMENTS = {'value', 'annotation'}

# Elements whose text should be stored as char (string) values
CHAR_ELEMENTS = {'spaceGroup'}


class Command(BaseCommand):
    """Backfill KPIs from params.xml for finished jobs."""

    help = "Backfill job KPIs from params.xml files"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            "--commit",
            action="store_true",
            help="Actually write to the database (default is dry run)",
        )
        parser.add_argument(
            "-pn", "--projectname",
            type=str,
            help="Limit to a specific project",
        )
        parser.add_argument(
            "--missing-only",
            action="store_true",
            help="Only backfill jobs that currently have no KPIs",
        )
        parser.add_argument(
            "--tasks",
            nargs="+",
            type=str,
            help="Limit to specific task names (e.g. aimless_pipe xia2_dials)",
        )

    def handle(self, *args, **options):
        commit = options["commit"]
        missing_only = options["missing_only"]

        if not commit:
            self.stdout.write(self.style.WARNING(
                "DRY RUN — pass --commit to write to the database\n"
            ))

        # Build queryset
        jobs = Job.objects.filter(status=Job.Status.FINISHED).select_related("project")

        if options["projectname"]:
            try:
                project = Project.objects.get(name=options["projectname"])
            except Project.DoesNotExist:
                raise CommandError(f"Project '{options['projectname']}' not found")
            jobs = jobs.filter(project=project)

        if options["tasks"]:
            jobs = jobs.filter(task_name__in=options["tasks"])

        if missing_only:
            # Exclude jobs that already have float or char values
            jobs = jobs.exclude(float_values__isnull=False).exclude(char_values__isnull=False)

        jobs = jobs.order_by("project__name", "number")

        total_jobs = 0
        total_kpis = 0
        skipped_no_params = 0
        skipped_no_perf = 0

        for job in jobs.iterator():
            params_path = job.directory / "params.xml"
            if not params_path.exists():
                skipped_no_params += 1
                continue

            kpis = self._extract_kpis_from_xml(params_path)
            if not kpis:
                skipped_no_perf += 1
                continue

            total_jobs += 1
            total_kpis += len(kpis)

            self.stdout.write(
                f"  {job.project.name}/{job.number} ({job.task_name}): "
                f"{len(kpis)} KPIs — {', '.join(f'{k}={v}' for k, v in kpis.items())}"
            )

            if commit:
                self._register_kpis(job, kpis)

        # Summary
        self.stdout.write("")
        self.stdout.write(self.style.SUCCESS(f"{'='*60}"))
        action = "Backfilled" if commit else "Would backfill"
        self.stdout.write(self.style.SUCCESS(
            f"{action} {total_kpis} KPIs across {total_jobs} jobs"
        ))
        self.stdout.write(f"  Skipped (no params.xml): {skipped_no_params}")
        self.stdout.write(f"  Skipped (no PERFORMANCE element): {skipped_no_perf}")
        self.stdout.write(self.style.SUCCESS(f"{'='*60}\n"))

    def _extract_kpis_from_xml(self, params_path: Path) -> dict:
        """Parse params.xml and extract KPI values from PERFORMANCE* element."""
        try:
            tree = etree.parse(str(params_path))
        except Exception:
            return {}

        # Find the outputData element
        output_data = tree.find('.//outputData')
        if output_data is None:
            return {}

        # Find the performance element
        perf_elem = None
        for name in KPI_ELEMENT_NAMES:
            perf_elem = output_data.find(name)
            if perf_elem is not None:
                break

        if perf_elem is None:
            return {}

        # Extract KPI values from child elements
        kpis = {}
        for child in perf_elem:
            tag = child.tag
            if tag in SKIP_ELEMENTS:
                continue

            text = (child.text or '').strip()
            if not text:
                continue

            if tag in CHAR_ELEMENTS:
                if len(text) > 0:
                    kpis[tag] = text
            else:
                try:
                    val = float(text)
                    # Skip zero values — they typically mean "not set"
                    if val != 0.0:
                        kpis[tag] = val
                except ValueError:
                    # Non-numeric, store as string
                    if len(text) > 0:
                        kpis[tag] = text

        return kpis

    def _register_kpis(self, job: Job, kpis: dict):
        """Write KPI values to the database."""
        with transaction.atomic():
            for key, value in kpis.items():
                job_value_key, _ = JobValueKey.objects.get_or_create(
                    name=key,
                    defaults={"description": key},
                )

                if isinstance(value, float):
                    JobFloatValue.objects.update_or_create(
                        job=job,
                        key=job_value_key,
                        defaults={"value": value},
                    )
                elif isinstance(value, str):
                    JobCharValue.objects.update_or_create(
                        job=job,
                        key=job_value_key,
                        defaults={"value": value},
                    )
