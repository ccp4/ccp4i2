"""
Django management command to get job reports (params, report, diagnostic XML).

Usage:
    python manage.py get_job_report --jobuuid <uuid> --type params
    python manage.py get_job_report --jobuuid <uuid> --type report
    python manage.py get_job_report --jobuuid <uuid> --type diagnostic
"""

import uuid
from django.core.management.base import BaseCommand, CommandError
from ccp4x.db.models import Job, Project
from ccp4x.lib.utils.jobs.reports import (
    get_job_params_xml,
    get_job_report_xml,
    get_job_diagnostic_xml
)


class Command(BaseCommand):
    """Get job reports in XML format."""

    help = "Get job reports (params/report/diagnostic XML)"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument(
            "--type",
            choices=["params", "report", "diagnostic"],
            default="params",
            help="Type of report to retrieve"
        )
        parser.add_argument("-o", "--output", help="Output file path", type=str)
        parser.add_argument("--regenerate", action="store_true", help="Regenerate cached report")

    def handle(self, *args, **options):
        # Find job
        try:
            the_job = self.get_job(options)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            raise CommandError(str(e))

        report_type = options["type"]
        output_path = options.get("output")
        regenerate = options.get("regenerate", False)

        # Get appropriate report
        if report_type == "params":
            result = get_job_params_xml(the_job)
        elif report_type == "report":
            result = get_job_report_xml(the_job, regenerate=regenerate)
        elif report_type == "diagnostic":
            result = get_job_diagnostic_xml(the_job)

        if result.success:
            content = result.data
            if output_path:
                mode = 'wb' if isinstance(content, bytes) else 'w'
                with open(output_path, mode) as f:
                    f.write(content)
                self.stdout.write(self.style.SUCCESS(
                    f"✓ {report_type.capitalize()} report written to {output_path}"
                ))
            else:
                # Output to stdout
                if isinstance(content, bytes):
                    self.stdout.write(content.decode('utf-8'))
                else:
                    self.stdout.write(content)
        else:
            raise CommandError(f"Failed to get {report_type} report: {result.error}")

    def get_job(self, options):
        if options["jobid"] is not None:
            return Job.objects.get(id=options["jobid"])
        if options["jobuuid"] is not None:
            return Job.objects.get(uuid=uuid.UUID(options["jobuuid"]))
        if options["projectname"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        raise Job.DoesNotExist("No job found with the provided criteria.")
