import uuid
from django.core.management.base import BaseCommand
from ccp4x.db.models import Job, Project


class Command(BaseCommand):
    """
    A Django management command to set the status of a job.

    Attributes:
        help (str): Description of the command.
        requires_system_checks (list): List of system checks required before running the command.

    Methods:
        add_arguments(parser):
            Adds command-line arguments to the parser.

        handle(*args, **options):
            Handles the command execution. Retrieves the job based on provided options and sets its status.

        get_job(options):
            Retrieves the job based on the provided options. Raises Job.DoesNotExist if no job is found.
    """

    help = "Set the status of a job"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-pi", "--projectid", help="Integer project id", type=int)
        parser.add_argument("-pu", "--projectuuid", help="Project uuid", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument(
            "-s",
            "--status",
            help="Job status",
            choices=[
                "UNKNOWN",
                "PENDING",
                "QUEUED",
                "RUNNING",
                "INTERRUPTED",
                "FAILED",
                "FINISHED",
                "RUNNING_REMOTELY",
                "FILE_HOLDER",
                "TO_DELETE",
                "UNSATISFACTORY",
            ],
            required=True,
        )

    def handle(self, *args, **options):
        try:
            the_job = self.get_job(options)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            self.stderr.write(self.style.ERROR(str(e)))
            return

        # Get the status value from the Status enum
        status_value = options["status"]
        try:
            new_status = getattr(Job.Status, status_value)
        except AttributeError:
            self.stderr.write(self.style.ERROR(f"Invalid status: {status_value}"))
            return

        # Set the job status
        old_status = the_job.status
        the_job.status = new_status
        the_job.save()

        self.stdout.write(
            self.style.SUCCESS(
                f"Job {the_job.id} ({the_job.uuid}) status changed from "
                f"{old_status} to {new_status}"
            )
        )

    def get_job(self, options):
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
