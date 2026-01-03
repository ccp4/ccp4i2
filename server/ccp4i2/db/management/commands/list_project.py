import json
from django.core.management.base import BaseCommand
from ccp4i2.db.models import Project
from ccp4i2.lib.utils.navigation.list_project import list_project


class Command(BaseCommand):
    """
    A Django management command to run a job.

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

    help = "Import a project"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-pi", "--projectid", help="Integer project id", type=int)
        parser.add_argument("-pu", "--projectuuid", help="Project uuid", type=str)

    def handle(self, *args, **options):
        try:
            the_project = self.get_project(options)
        except Project.DoesNotExist as e:
            self.stderr.write(self.style.ERROR(str(e)))
            return

        project_files = list_project(str(the_project.uuid))
        self.stdout.write(json.dumps(project_files, indent=2))

    def get_project(self, options):
        if options["projectname"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return the_project
        if options["projectid"] is not None:
            the_project = Project.objects.get(id=options["projectid"])
            return the_project
        if options["projectuuid"]:
            the_project = Project.objects.get(uuid=uuid.UUID(options["projectuuid"]))
            return the_project
        raise Project.DoesNotExist("No project found with the provided criteria.")
