"""
Django management command to set job parameters.

Usage:
    python manage.py set_job_parameter --jobuuid <uuid> --path "inputData.XYZIN" --value "/path/to/file.pdb"
    python manage.py set_job_parameter --jobid <id> --path "container.parameter" --value "some_value"
"""

import uuid
import json
from django.core.management.base import BaseCommand, CommandError
from ccp4i2.db.models import Job, Project
from ccp4i2.lib.utils.parameters.set_param import set_parameter


class Command(BaseCommand):
    """Set job parameter values."""

    help = "Set job parameter"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument(
            "--path",
            required=True,
            help="Object path (e.g., 'inputData.XYZIN.fileName')"
        )
        parser.add_argument(
            "--value",
            required=True,
            help="New value (string, number, or JSON object)"
        )
        parser.add_argument(
            "--type",
            choices=["string", "int", "float", "json", "auto"],
            default="auto",
            help="Value type (auto-detect by default)"
        )
        parser.add_argument("--json-output", action="store_true", help="Output as JSON")

    def handle(self, *args, **options):
        json_output = options.get("json_output", False)

        # Find job
        try:
            the_job = self.get_job(options)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            if json_output:
                self.stdout.write(json.dumps({"status": "Failed", "reason": str(e)}))
            else:
                raise CommandError(str(e))
            return

        # Parse value
        path = options["path"]
        value_str = options["value"]
        value_type = options["type"]

        try:
            if value_type == "int":
                value = int(value_str)
            elif value_type == "float":
                value = float(value_str)
            elif value_type == "json":
                value = json.loads(value_str)
            elif value_type == "auto":
                # Try to auto-detect
                try:
                    value = json.loads(value_str)
                except (json.JSONDecodeError, ValueError):
                    # It's a string
                    value = value_str
            else:  # string
                value = value_str
        except Exception as e:
            raise CommandError(f"Failed to parse value: {e}")

        # Set parameter
        result = set_parameter(the_job, path, value)

        if result.success:
            updated_obj = result.data

            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Success",
                    "job_uuid": str(the_job.uuid),
                    "job_number": the_job.number,
                    "parameter_path": path,
                    "updated_object": updated_obj
                }, indent=2))
            else:
                self.stdout.write(self.style.SUCCESS(
                    f"✓ Set parameter '{path}' = {value}"
                ))
                self.stdout.write(f"  Job: {the_job.number} ({the_job.task_name})")
                self.stdout.write(f"  Updated object: {json.dumps(updated_obj, indent=2)}")
        else:
            if json_output:
                self.stdout.write(json.dumps(result.to_dict(), indent=2))
            else:
                raise CommandError(f"Failed to set parameter: {result.error}")

    def get_job(self, options):
        if options["jobid"] is not None:
            return Job.objects.get(id=options["jobid"])
        if options["jobuuid"] is not None:
            return Job.objects.get(uuid=uuid.UUID(options["jobuuid"]))
        if options["projectname"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        raise Job.DoesNotExist("No job found with the provided criteria.")
