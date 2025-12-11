import uuid
import subprocess
import pathlib
from django.core.management.base import BaseCommand
from ccp4i2.db.models import Project, File
from ccp4i2.lib.utils.files.preview import preview_file


class Command(BaseCommand):
    """
    Command class to preview a file.

    This command provides functionality to preview a file based on various criteria such as
    project name, project ID, project UUID, file path, file ID, and file UUID. It also supports
    detaching the job to run in the background.

    Attributes:
        help (str): Description of the command.
        requires_system_checks (list): List of system checks required before running the command.

    Methods:
        add_arguments(parser):
            Adds command-line arguments to the parser.

        handle(*args, **options):
            Handles the command execution. Retrieves the file path based on the provided options
            and either runs the preview in the foreground or detaches it to run in the background.

        get_path(options):
            Retrieves the file path based on the provided options. Raises File.DoesNotExist if no
            file is found with the provided criteria.
    """

    help = "Preview a file"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-pi", "--projectid", help="Integer project id", type=int)
        parser.add_argument("-pu", "--projectuuid", help="Project uuid", type=str)
        parser.add_argument(
            "-r", "--rel_path", help="File path relative to project root", type=str
        )
        parser.add_argument(
            "-p", "--path", help="Full file path of file to preview", type=str
        )
        parser.add_argument("-fi", "--file_id", help="File id", type=int)
        parser.add_argument("-fu", "--file_uuid", help="FIle uuid", type=str)
        parser.add_argument(
            "-e", "--executable", help="Executable of viewer to use", type=str
        )

    def handle(self, *args, **options):
        try:
            file_path = self.get_path(options)
        except (File.DoesNotExist, Project.DoesNotExist) as e:
            self.stderr.write(self.style.ERROR(str(e)))
            return
        preview_file(options["executable"], str(file_path))

    def get_path(self, options):
        if options["path"] is not None:
            return pathlib.Path(options["path"])
        if options["file_id"] is not None:
            return File.objects.get(id=options["file_id"]).path
        if options["file_uuid"] is not None:
            return File.objects.get(uuid=uuid.UUID(options["file_uuid"])).path
        if options["projectname"] is not None and options["rel_path"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return pathlib.Path(the_project.directory) / options["rel_path"]
        if options["projectid"] is not None and options["rel_path"] is not None:
            the_project = Project.objects.get(id=options["projectid"])
            return pathlib.Path(the_project.directory) / options["rel_path"]
        if options["projectuuid"] is not None and options["rel_path"] is not None:
            the_project = Project.objects.get(uuid=uuid.UUID(options["projectuuid"]))
            return pathlib.Path(the_project.directory) / options["rel_path"]
        raise File.DoesNotExist("No file found with the provided criteria.")
