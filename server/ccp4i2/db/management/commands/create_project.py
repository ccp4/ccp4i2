"""
Django management command to create a new CCP4i2 project.

This command uses the ProjectSerializer to create a project with all
necessary directory structure (CCP4_JOBS, CCP4_IMPORTED_FILES, etc.).

Usage:
    python manage.py create_project <project_name>
    python manage.py create_project my_project --description "My crystallography project"
    python manage.py create_project test_project --directory /custom/path/test_project

Example:
    python manage.py create_project toxd_project
    python manage.py create_project toxd_project --description "Toxd structure determination"
"""

from django.core.management.base import BaseCommand, CommandError
from ccp4i2.api.serializers import ProjectSerializer


class Command(BaseCommand):
    """
    Django management command to create a new CCP4i2 project.

    The serializer automatically:
    - Creates the project directory (defaults to CCP4I2_PROJECTS_DIR/project_name)
    - Creates subdirectories: CCP4_JOBS, CCP4_IMPORTED_FILES, CCP4_COOT, CCP4_TMP, CCP4_PROJECT_FILES
    - Validates project name (alphanumeric, underscore, hyphen only)
    - Checks for duplicate project names
    - Verifies write permissions

    Attributes:
        help (str): Description of the command displayed with --help
        requires_system_checks (list): System checks required before running
    """

    help = "Create a new CCP4i2 project with directory structure"
    requires_system_checks = []

    def add_arguments(self, parser):
        """
        Add command-line arguments.

        Args:
            parser: Django argument parser
        """
        # Required argument: project name
        parser.add_argument(
            "name",
            type=str,
            help="Project name (alphanumeric, underscore, hyphen only)",
        )

        # Optional arguments
        parser.add_argument(
            "-d",
            "--description",
            type=str,
            default="",
            help="Project description (optional)",
        )

        parser.add_argument(
            "--directory",
            type=str,
            default="__default__",
            help="Custom project directory path (optional, defaults to CCP4I2_PROJECTS_DIR/name)",
        )

        parser.add_argument(
            "-q",
            "--quiet",
            action="store_true",
            help="Suppress output (only print project UUID)",
        )

        parser.add_argument(
            "--json",
            action="store_true",
            help="Output project details as JSON",
        )

    def handle(self, *args, **options):
        """
        Handle the create_project command.

        Creates a new project using ProjectSerializer, which handles:
        - Directory creation with proper structure
        - Name validation
        - Duplicate checking

        Args:
            *args: Positional arguments (unused)
            **options: Command options including project name and description

        Raises:
            CommandError: If project creation fails (validation, permissions, etc.)
        """
        import json

        project_name = options["name"]
        description = options["description"]
        directory = options["directory"]
        quiet = options["quiet"]
        json_output = options["json"]

        # Prepare data for serializer
        data = {
            "name": project_name,
            "description": description,
            "directory": directory,
        }

        # Create project using serializer
        serializer = ProjectSerializer(data=data)

        try:
            # Validate
            serializer.is_valid(raise_exception=True)

            # Create project (this also creates the directory structure)
            project = serializer.save()

            # Output results
            if json_output:
                # JSON output for programmatic use
                output = {
                    "uuid": str(project.uuid),
                    "name": project.name,
                    "description": project.description,
                    "directory": project.directory,
                    "creation_time": project.creation_time.isoformat() if project.creation_time else None,
                }
                self.stdout.write(json.dumps(output, indent=2))

            elif quiet:
                # Quiet mode - just output UUID for scripting
                self.stdout.write(str(project.uuid))

            else:
                # Normal verbose output
                self.stdout.write(self.style.SUCCESS(f"\n{'='*70}"))
                self.stdout.write(self.style.SUCCESS("Project Created Successfully"))
                self.stdout.write(self.style.SUCCESS(f"{'='*70}\n"))

                self.stdout.write(f"  Project UUID:  {project.uuid}")
                self.stdout.write(f"  Name:          {project.name}")
                if project.description:
                    self.stdout.write(f"  Description:   {project.description}")
                self.stdout.write(f"  Directory:     {project.directory}")
                self.stdout.write(f"  Created:       {project.creation_time}")

                self.stdout.write(f"\n  Subdirectories created:")
                for subdir in [
                    "CCP4_JOBS",
                    "CCP4_IMPORTED_FILES",
                    "CCP4_COOT",
                    "CCP4_TMP",
                    "CCP4_PROJECT_FILES",
                ]:
                    self.stdout.write(f"    - {subdir}/")

                self.stdout.write(f"\n{'='*70}\n")

        except Exception as e:
            # Handle validation or creation errors
            if hasattr(e, "detail"):
                # DRF validation error
                error_msg = f"Project creation failed:\n"
                if isinstance(e.detail, dict):
                    for field, errors in e.detail.items():
                        if isinstance(errors, list):
                            for error in errors:
                                error_msg += f"  {field}: {error}\n"
                        else:
                            error_msg += f"  {field}: {errors}\n"
                else:
                    error_msg += f"  {e.detail}"
                raise CommandError(error_msg)
            else:
                # Other errors
                raise CommandError(f"Project creation failed: {str(e)}")
