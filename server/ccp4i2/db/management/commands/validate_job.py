"""
Django management command to validate a CCP4i2 job's parameters.

Usage:
    python manage.py validate_job --jobuuid <uuid>
    python manage.py validate_job --jobid <id>
    python manage.py validate_job --projectname <name> --jobnumber <num>
"""

import uuid
from django.core.management.base import BaseCommand, CommandError
from ccp4i2.db.models import Job, Project
from ccp4i2.lib.utils.jobs.validate import validate_job
from xml.etree import ElementTree as ET


class Command(BaseCommand):
    """Validate job parameters and show validation errors."""

    help = "Validate job parameters"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument("--json", action="store_true", help="Output as JSON")
        parser.add_argument("-o", "--output", help="Output file path", type=str)

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

        # Validate
        result = validate_job(the_job)

        if result.success:
            error_tree = result.data
            errors = error_tree.findall('.//error')
            warnings = error_tree.findall('.//warning')

            output_path = options.get("output")
            if output_path:
                with open(output_path, 'wb') as f:
                    f.write(ET.tostring(error_tree))
                self.stdout.write(self.style.SUCCESS(f"Validation XML written to {output_path}"))

            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Success",
                    "error_count": len(errors),
                    "warning_count": len(warnings),
                    "has_errors": len(errors) > 0
                }, indent=2))
            else:
                if len(errors) == 0 and len(warnings) == 0:
                    self.stdout.write(self.style.SUCCESS(
                        f"✓ Job {the_job.number} validation passed - no errors or warnings"
                    ))
                else:
                    if len(errors) > 0:
                        self.stdout.write(self.style.ERROR(
                            f"✗ Found {len(errors)} validation error(s)"
                        ))
                        for err in errors[:5]:  # Show first 5
                            msg = err.findtext('message', 'No message')
                            self.stdout.write(f"  - {msg}")
                        if len(errors) > 5:
                            self.stdout.write(f"  ... and {len(errors) - 5} more")

                    if len(warnings) > 0:
                        self.stdout.write(self.style.WARNING(
                            f"! Found {len(warnings)} warning(s)"
                        ))

                    if output_path:
                        self.stdout.write(f"\nFull report: {output_path}")
        else:
            if json_output:
                self.stdout.write(json.dumps(result.to_dict(), indent=2))
            else:
                raise CommandError(f"Validation failed: {result.error}")

    def get_job(self, options):
        if options["jobid"] is not None:
            return Job.objects.get(id=options["jobid"])
        if options["jobuuid"] is not None:
            return Job.objects.get(uuid=uuid.UUID(options["jobuuid"]))
        if options["projectname"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        raise Job.DoesNotExist("No job found with the provided criteria.")
