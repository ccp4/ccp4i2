import uuid
import logging
import asyncio
from django.core.management.base import BaseCommand
from ccp4i2.db.models import Job, Project
from ccp4i2.lib.async_create_job import create_job_async

logger = logging.getLogger(f"ccp4i2:{__name__}")


class Command(BaseCommand):
    """
    A Django management command to create a new job with a given task name and optional context job to define inputs.

    Attributes:
        help (str): Description of the command.
        requires_system_checks (list): List of system checks required before running the command.

    Methods:
        add_arguments(parser):
            Adds command-line arguments to the parser.

        handle(*args, **options):
            Handles the command execution. Retrieves the job based on provided options and runs it.
            If the detach option is specified, the job is run in a detached subprocess.

        get_job(options):
            Retrieves the job based on the provided options. Raises Job.DoesNotExist if no job is found.
    """

    help = "Clone a job"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            "-cn", "--contextjobnumber", help="Context job number", type=str
        )
        parser.add_argument(
            "-ci", "--contextjobid", help="Context integer job id", type=int
        )
        parser.add_argument(
            "-cu", "--contextjobuuid", help="Context job UUID", type=str
        )
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-pi", "--projectid", help="Integer project id", type=int)
        parser.add_argument("-pu", "--projectuuid", help="Project uuid", type=str)
        parser.add_argument(
            "-tn", "--taskname", help="Plugin name for new job", type=str
        )

    def handle(self, *args, **options):
        the_context_job = None
        try:
            the_context_job = self.get_context_job(options)
        except Project.DoesNotExist as e:
            self.stderr.write(self.style.ERROR(str(e)))
            return
        except Job.DoesNotExist as e:
            logger.info("Unable to identify context job %s", e)
            the_context_job = None

        try:
            the_project = self.get_project(options)
        except Project.DoesNotExist as e:
            self.stderr.write(self.style.ERROR(str(e)))
            return

        # Use async create_job_async function
        parent_job_uuid = the_context_job.uuid if the_context_job is not None else None

        try:
            result = asyncio.run(
                create_job_async(
                    project_uuid=the_project.uuid,
                    task_name=options["taskname"],
                    parent_job_uuid=parent_job_uuid,
                    save_params=True
                )
            )

            # Result contains job info directly
            self.stdout.write(
                f"Created job with number {result['job_number']}, "
                f"uuid {result['job_uuid']}"
            )
        except Exception as e:
            self.stderr.write(
                self.style.ERROR(f"Failed to create job: {e}")
            )
            import traceback
            traceback.print_exc()

    def get_context_job(self, options):
        if options["contextjobid"] is not None:
            return Job.objects.get(id=options["contextjobid"])
        if options["contextjobuuid"] is not None:
            return Job.objects.get(uuid=uuid.UUID(options["contextjobuuid"]))
        if (
            options["projectname"] is not None
            and options["contextjobnumber"] is not None
        ):
            return Job.objects.get(
                number=options["contextjobnumber"], project__name=options["projectname"]
            )
        if options["projectid"] is not None and options["contextjobnumber"] is not None:
            return Job.objects.get(
                number=options["contextjobnumber"], project__id=options["projectid"]
            )
        if (
            options["projectuuid"] is not None
            and options["contextjobnumber"] is not None
        ):
            return Job.objects.get(
                number=options["contextjobnumber"], project__uuid=options["projectuuid"]
            )
        raise Job.DoesNotExist("No job found with the provided criteria.")

    def get_project(self, options):
        if options["projectname"] is not None:
            return Project.objects.get(name=options["projectname"])
        if options["projectid"] is not None:
            return Project.objects.get(id=options["projectid"])
        if options["projectuuid"] is not None:
            return Project.objects.get(uuid=options["projectuuid"])
        raise Project.DoesNotExist("No project found with the provided criteria.")
