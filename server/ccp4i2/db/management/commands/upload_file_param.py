"""
Django management command to upload a file and set it as a job parameter.

Usage:
    python manage.py upload_file_param -pn myproject -jn 5 --path "inputData.XYZIN" --file /path/to/file.pdb
    python manage.py upload_file_param -pn myproject -jn 5 --path "inputData.F_SIGF" --file /path/to/data.mtz --column-selector "/*/*/[FP,SIGFP]"
"""

import uuid
import json
import pathlib
from django.core.management.base import BaseCommand, CommandError
from django.http import QueryDict
from django.core.files.uploadedfile import SimpleUploadedFile
from ccp4i2.db.models import Job, Project
from ccp4i2.lib.utils.files.upload_param import upload_file_param


class MockRequest:
    """Mock HTTP request for upload_file_param utility."""

    def __init__(self, file_path: pathlib.Path, object_path: str, column_selector: str = None):
        self.POST = QueryDict(mutable=True)
        self.POST["objectPath"] = object_path
        if column_selector:
            self.POST["column_selector"] = column_selector

        # Read file and create uploaded file object
        with open(file_path, "rb") as f:
            file_content = f.read()

        self.FILES = {
            "file": [SimpleUploadedFile(file_path.name, file_content)]
        }


class Command(BaseCommand):
    """Upload a file and set it as a job parameter."""

    help = "Upload file and set job parameter"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument(
            "--path",
            required=True,
            help="Object path (e.g., 'inputData.XYZIN')"
        )
        parser.add_argument(
            "--file",
            required=True,
            help="Path to file to upload"
        )
        parser.add_argument(
            "--column-selector",
            help="MTZ column selector for reflection files (e.g., '/*/*/[FP,SIGFP]')"
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

        # Validate file exists
        file_path = pathlib.Path(options["file"])
        if not file_path.exists():
            error_msg = f"File not found: {file_path}"
            if json_output:
                self.stdout.write(json.dumps({"status": "Failed", "reason": error_msg}))
            else:
                raise CommandError(error_msg)
            return

        object_path = options["path"]
        column_selector = options.get("column_selector")

        if not json_output:
            self.stdout.write(f"Uploading {file_path} to {object_path}...")
            if column_selector:
                self.stdout.write(f"  Column selector: {column_selector}")

        # Create mock request
        request = MockRequest(file_path, object_path, column_selector)

        try:
            result = upload_file_param(the_job, request)

            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Success",
                    "job_uuid": str(the_job.uuid),
                    "job_number": the_job.number,
                    "parameter_path": object_path,
                    "updated_item": result
                }, indent=2))
            else:
                self.stdout.write(self.style.SUCCESS("Upload successful:"))
                self.stdout.write(f"  Base name: {result.get('baseName', 'N/A')}")
                self.stdout.write(f"  Annotation: {result.get('annotation', 'N/A')}")
                if result.get('dbFileId'):
                    self.stdout.write(f"  DB ID: {result.get('dbFileId')}")

        except Exception as e:
            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Failed",
                    "reason": str(e)
                }, indent=2))
            else:
                raise CommandError(f"Upload failed: {e}")

    def get_job(self, options):
        if options["jobid"] is not None:
            return Job.objects.get(id=options["jobid"])
        if options["jobuuid"] is not None:
            return Job.objects.get(uuid=uuid.UUID(options["jobuuid"]))
        if options["projectname"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        raise Job.DoesNotExist("No job found with the provided criteria.")
