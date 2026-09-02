import platform
import subprocess
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from ccp4i2.db.import_i2xml import ProjectArchiveError, import_ccp4_project_zip


class Command(BaseCommand):

    help = "Import a project"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("zip_file", nargs="*")
        parser.add_argument("-d", "--detach", help="Detach job", action="store_true")

    def handle(self, *args, **options):
        zip_files = options["zip_file"]
        if not zip_files:
            raise CommandError("No zip file to import")

        if options["detach"]:
            # Determine the program name based on the OS
            ccp4_python_program = "ccp4-python"
            if platform.system() == "Windows":
                ccp4_python_program += ".bat"
            # One child per archive, each argument passed separately: joining
            # them into a single argv element made every import after the first
            # arrive as part of one unopenable path.
            for zip_file in zip_files:
                _ = subprocess.Popen(
                    [
                        ccp4_python_program,
                        "manage.py",
                        "import_ccp4_project_zip",
                        zip_file,
                    ],
                    start_new_session=True,
                )
        else:
            for zip_file in zip_files:
                try:
                    import_ccp4_project_zip(
                        zip_file, relocate_path=settings.CCP4I2_PROJECTS_DIR
                    )
                except ProjectArchiveError as err:
                    # A bad archive is the user's problem to fix, not a
                    # traceback: say what is wrong with it and stop.
                    raise CommandError(f"{zip_file}: {err}") from err
