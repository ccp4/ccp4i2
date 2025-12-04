"""
Django management command to clone a CCP4i2 job.

Usage:
    python manage.py clone_job --jobuuid <uuid>
    python manage.py clone_job --jobid <id>
    python manage.py clone_job --projectname <name> --jobnumber <num>
"""

import uuid
from django.core.management.base import BaseCommand, CommandError
from ccp4x.db.models import Job, Project
from ccp4x.lib.utils.jobs.clone import clone_job


class Command(BaseCommand):
    """
    Clone an existing job with identical parameters.

    Creates a new job as a copy of an existing job, preserving all input
    parameters but resetting status to PENDING and generating new UUIDs.
    """

    help = "Clone a job with identical parameters"
    requires_system_checks = []

    def add_arguments(self, parser):
        """Add command-line arguments."""
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-pi", "--projectid", help="Integer project id", type=int)
        parser.add_argument("-pu", "--projectuuid", help="Project uuid", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument(
            "--json",
            action="store_true",
            help="Output result as JSON"
        )

    def handle(self, *args, **options):
        """Execute the clone operation."""
        import json

        json_output = options.get("json", False)

        # Find the job to clone
        try:
            the_job = self.get_job(options)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Failed",
                    "reason": str(e)
                }))
            else:
                raise CommandError(str(e))
            return

        # Clone the job
        result = clone_job(str(the_job.uuid))

        # Handle result
        if result.success:
            new_job = result.data
            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Success",
                    "original_job_uuid": str(the_job.uuid),
                    "original_job_number": the_job.number,
                    "new_job_uuid": str(new_job.uuid),
                    "new_job_number": new_job.number,
                    "new_job_id": new_job.id,
                    "task_name": new_job.task_name,
                    "project_name": new_job.project.name
                }, indent=2))
            else:
                self.stdout.write(self.style.SUCCESS(
                    f"Successfully cloned job {the_job.number} to new job {new_job.number}"
                ))
                self.stdout.write(f"  Original UUID: {the_job.uuid}")
                self.stdout.write(f"  New UUID:      {new_job.uuid}")
                self.stdout.write(f"  New Job ID:    {new_job.id}")
                self.stdout.write(f"  Task:          {new_job.task_name}")
                self.stdout.write(f"  Project:       {new_job.project.name}")
        else:
            if json_output:
                self.stdout.write(json.dumps(result.to_dict(), indent=2))
            else:
                raise CommandError(f"Clone failed: {result.error}")

    def get_job(self, options):
        """
        Retrieve job based on provided options.

        Returns:
            Job: The job instance

        Raises:
            Job.DoesNotExist: If no job found with provided criteria
            Project.DoesNotExist: If project not found (when searching by project)
        """
        if options["jobid"] is not None:
            return Job.objects.get(id=options["jobid"])

        if options["jobuuid"] is not None:
            return Job.objects.get(uuid=uuid.UUID(options["jobuuid"]))

        if options["projectname"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)

        if options["projectid"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(id=options["projectid"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)

        if options["projectuuid"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(uuid=uuid.UUID(options["projectuuid"]))
            return Job.objects.get(number=options["jobnumber"], project=the_project)

        raise Job.DoesNotExist("No job found with the provided criteria.")
