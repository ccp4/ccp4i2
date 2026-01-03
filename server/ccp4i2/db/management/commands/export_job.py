"""
Django management command to export a job as ZIP archive.

Usage:
    python manage.py export_job --jobuuid <uuid> --output /path/to/export.zip
    python manage.py export_job --jobid <id> -o job_export.zip
"""

import uuid
from pathlib import Path
from django.core.management.base import BaseCommand, CommandError
from ccp4i2.db.models import Job, Project
from ccp4i2.lib.utils.jobs.export import export_job


class Command(BaseCommand):
    """Export a job as a ZIP archive."""

    help = "Export a job as ZIP archive"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument(
            "-o", "--output",
            required=True,
            help="Output ZIP file path"
        )
        parser.add_argument("--json", action="store_true", help="Output as JSON")

    def handle(self, *args, **options):
        import json

        json_output = options.get("json", False)

        # Find job
        try:
            the_job = self.get_job(options)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            if json_output:
                self.stdout.write(json.dumps({"status": "Failed", "reason": str(e)}))
            else:
                raise CommandError(str(e))
            return

        output_path = Path(options["output"])

        # Export job
        result = export_job(the_job, output_path)

        if result.success:
            zip_path = result.data

            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Success",
                    "job_uuid": str(the_job.uuid),
                    "job_number": the_job.number,
                    "task_name": the_job.task_name,
                    "export_path": str(zip_path),
                    "file_size_bytes": zip_path.stat().st_size
                }, indent=2))
            else:
                size_mb = zip_path.stat().st_size / (1024 * 1024)
                self.stdout.write(self.style.SUCCESS(
                    f"✓ Job {the_job.number} exported successfully"
                ))
                self.stdout.write(f"  Task: {the_job.task_name}")
                self.stdout.write(f"  Output: {zip_path}")
                self.stdout.write(f"  Size: {size_mb:.2f} MB")
        else:
            if json_output:
                self.stdout.write(json.dumps(result.to_dict(), indent=2))
            else:
                raise CommandError(f"Failed to export job: {result.error}")

    def get_job(self, options):
        if options["jobid"] is not None:
            return Job.objects.get(id=options["jobid"])
        if options["jobuuid"] is not None:
            return Job.objects.get(uuid=uuid.UUID(options["jobuuid"]))
        if options["projectname"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        raise Job.DoesNotExist("No job found with the provided criteria.")
