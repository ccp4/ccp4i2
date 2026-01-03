"""
Django management command to execute a CCP4i2 job.

Usage:
    python manage.py execute_job --jobuuid <uuid>
    python manage.py execute_job --jobid <id>
    python manage.py execute_job --projectname <name> --jobnumber <num>
    python manage.py execute_job --jobuuid <uuid> --force-local
"""

import uuid
from django.core.management.base import BaseCommand, CommandError
from ccp4i2.db.models import Job, Project
from ccp4i2.lib.utils.jobs.execute import execute_job


class Command(BaseCommand):
    """Execute a CCP4i2 job (context-aware: local or remote)."""

    help = "Execute a CCP4i2 job"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument(
            "--force-local",
            action="store_true",
            help="Force local execution (bypass Azure queue)"
        )
        parser.add_argument("--json", action="store_true", help="Output as JSON")

    def handle(self, *args, **options):
        import json

        json_output = options.get("json", False)
        force_local = options.get("force_local", False)

        # Find job
        try:
            the_job = self.get_job(options)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            if json_output:
                self.stdout.write(json.dumps({"status": "Failed", "reason": str(e)}))
            else:
                raise CommandError(str(e))
            return

        # Execute job
        result = execute_job(the_job, force_local=force_local)

        if result.success:
            updated_job = result.data
            mode = "local" if force_local else "auto"

            if json_output:
                self.stdout.write(json.dumps({
                    "status": "Success",
                    "job_uuid": str(updated_job.uuid),
                    "job_number": updated_job.number,
                    "job_status": updated_job.get_status_display(),
                    "task_name": updated_job.task_name,
                    "execution_mode": mode
                }, indent=2))
            else:
                self.stdout.write(self.style.SUCCESS(
                    f"✓ Job {updated_job.number} started successfully ({mode} mode)"
                ))
                self.stdout.write(f"  UUID: {updated_job.uuid}")
                self.stdout.write(f"  Task: {updated_job.task_name}")
                self.stdout.write(f"  Status: {updated_job.get_status_display()}")
                self.stdout.write(f"  Directory: {updated_job.directory}")
        else:
            if json_output:
                self.stdout.write(json.dumps(result.to_dict(), indent=2))
            else:
                raise CommandError(f"Failed to run job: {result.error}")

    def get_job(self, options):
        if options["jobid"] is not None:
            return Job.objects.get(id=options["jobid"])
        if options["jobuuid"] is not None:
            return Job.objects.get(uuid=uuid.UUID(options["jobuuid"]))
        if options["projectname"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        raise Job.DoesNotExist("No job found with the provided criteria.")
