import platform
import subprocess
from django.conf import settings
from django.core.management.base import BaseCommand

from ccp4x.db.import_i2xml import import_ccp4_project_zip


class Command(BaseCommand):

    help = "Import a project"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("zip_file", nargs="*")
        parser.add_argument("-d", "--detach", help="Detach job", action="store_true")

    def handle(self, *args, **options):
        if options["detach"]:
            # Determine the program name based on the OS
            ccp4_python_program = "ccp4-python"
            if platform.system() == "Windows":
                ccp4_python_program += ".bat"
            _ = subprocess.Popen(
                [
                    ccp4_python_program,
                    "manage.py",
                    "import_ccp4_project_zip",
                    f"{' '.join(options['zip_file'])}",
                ],
                start_new_session=True,
            )
        else:
            for zip_file in options["zip_file"]:
                import_ccp4_project_zip(
                    zip_file, relocate_path=settings.CCP4I2_PROJECTS_DIR
                )
