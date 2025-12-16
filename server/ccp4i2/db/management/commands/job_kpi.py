"""
Django management command to display job KPIs (key performance indicators).

Shows JobFloatValue and JobCharValue entries associated with a job.

Usage:
    python manage.py job_kpi --projectname toxd --jobnumber 5
    python manage.py job_kpi --jobid 42
    python manage.py job_kpi --jobuuid <uuid>
    python manage.py job_kpi --projectname toxd --jobnumber 5 --json
"""

import json
import uuid as uuid_module
from django.core.management.base import BaseCommand, CommandError
from ccp4i2.db.models import Job, Project, JobFloatValue, JobCharValue


class Command(BaseCommand):
    """Display KPIs (JobFloatValue and JobCharValue) for a job."""

    help = "Display job KPIs (float and char values)"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument(
            "--json",
            action="store_true",
            help="Output as JSON"
        )
        parser.add_argument(
            "--floats-only",
            action="store_true",
            help="Show only float values"
        )
        parser.add_argument(
            "--chars-only",
            action="store_true",
            help="Show only char values"
        )

    def handle(self, *args, **options):
        # Find job
        try:
            the_job = self.get_job(options)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            raise CommandError(str(e))

        json_output = options.get("json", False)
        floats_only = options.get("floats_only", False)
        chars_only = options.get("chars_only", False)

        # Get KPI values
        float_values = []
        char_values = []

        if not chars_only:
            float_values = list(
                JobFloatValue.objects.filter(job=the_job)
                .select_related('key')
                .order_by('key__name')
            )

        if not floats_only:
            char_values = list(
                JobCharValue.objects.filter(job=the_job)
                .select_related('key')
                .order_by('key__name')
            )

        if json_output:
            self._output_json(the_job, float_values, char_values)
        else:
            self._output_table(the_job, float_values, char_values)

    def _output_json(self, job, float_values, char_values):
        """Output KPIs as JSON."""
        output = {
            "job": {
                "id": job.id,
                "uuid": str(job.uuid),
                "number": job.number,
                "title": job.title,
                "task_name": job.task_name,
            },
            "float_values": [
                {
                    "key": fv.key.name,
                    "description": fv.key.description,
                    "value": fv.value,
                }
                for fv in float_values
            ],
            "char_values": [
                {
                    "key": cv.key.name,
                    "description": cv.key.description,
                    "value": cv.value,
                }
                for cv in char_values
            ],
        }
        self.stdout.write(json.dumps(output, indent=2))

    def _output_table(self, job, float_values, char_values):
        """Output KPIs as formatted table."""
        self.stdout.write(self.style.SUCCESS(f"\n{'='*80}"))
        self.stdout.write(self.style.SUCCESS(
            f"KPIs for Job #{job.number}: {job.title} ({job.task_name})"
        ))
        self.stdout.write(self.style.SUCCESS(f"{'='*80}\n"))

        if not float_values and not char_values:
            self.stdout.write("No KPI values found for this job.")
            return

        # Float values
        if float_values:
            self.stdout.write(self.style.MIGRATE_HEADING("Float Values:"))
            self.stdout.write(f"{'Key':<30} {'Value':<20} {'Description'}")
            self.stdout.write("-" * 80)
            for fv in float_values:
                key = fv.key.name[:29] if len(fv.key.name) > 29 else fv.key.name
                desc = fv.key.description[:30] if len(fv.key.description) > 30 else fv.key.description
                self.stdout.write(f"{key:<30} {fv.value:<20.6g} {desc}")
            self.stdout.write("")

        # Char values
        if char_values:
            self.stdout.write(self.style.MIGRATE_HEADING("Char Values:"))
            self.stdout.write(f"{'Key':<30} {'Value':<30} {'Description'}")
            self.stdout.write("-" * 80)
            for cv in char_values:
                key = cv.key.name[:29] if len(cv.key.name) > 29 else cv.key.name
                value = cv.value[:29] if len(cv.value) > 29 else cv.value
                desc = cv.key.description[:20] if len(cv.key.description) > 20 else cv.key.description
                self.stdout.write(f"{key:<30} {value:<30} {desc}")
            self.stdout.write("")

        # Summary
        self.stdout.write(f"{'='*80}")
        self.stdout.write(f"Total: {len(float_values)} float values, {len(char_values)} char values")
        self.stdout.write(f"{'='*80}\n")

    def get_job(self, options):
        """Find job by various identifiers."""
        if options["jobid"] is not None:
            return Job.objects.get(id=options["jobid"])
        if options["jobuuid"] is not None:
            return Job.objects.get(uuid=uuid_module.UUID(options["jobuuid"]))
        if options["projectname"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        raise CommandError(
            "Must specify job using --projectname and --jobnumber, --jobid, or --jobuuid"
        )
