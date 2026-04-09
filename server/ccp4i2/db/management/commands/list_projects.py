"""
Django management command to list CCP4i2 projects.

Usage:
    python manage.py list_projects
    python manage.py list_projects --json
    python manage.py list_projects --names-only
    python manage.py list_projects --sort-by creation_time
    python manage.py list_projects --filter "toxd*"

Example:
    python manage.py list_projects
    python manage.py list_projects --sort-by jobs --reverse
    python manage.py list_projects --json
"""

from django.core.management.base import BaseCommand, CommandError
from ccp4i2.db import models
from django.db.models import Count
import fnmatch


class Command(BaseCommand):
    """
    Django management command to list all CCP4i2 projects.

    Displays project information including:
    - UUID
    - Name
    - Description
    - Directory
    - Creation time
    - Job count
    - Last access time

    Supports multiple output formats and sorting/filtering options.

    Attributes:
        help (str): Description of the command displayed with --help
        requires_system_checks (list): System checks required before running
    """

    help = "List all CCP4i2 projects"
    requires_system_checks = []

    def add_arguments(self, parser):
        """
        Add command-line arguments.

        Args:
            parser: Django argument parser
        """
        # Output format options
        parser.add_argument(
            "--json",
            action="store_true",
            help="Output as JSON array",
        )

        parser.add_argument(
            "--names-only",
            action="store_true",
            help="Output only project names (one per line)",
        )

        # Sorting options
        parser.add_argument(
            "--sort-by",
            type=str,
            choices=["name", "creation_time", "last_access", "jobs"],
            default="creation_time",
            help="Sort projects by field (default: creation_time)",
        )

        parser.add_argument(
            "--reverse",
            action="store_true",
            help="Reverse sort order",
        )

        # Filtering options
        parser.add_argument(
            "--filter",
            type=str,
            help="Filter projects by name pattern (supports wildcards like 'toxd*')",
        )

        parser.add_argument(
            "--limit",
            type=int,
            help="Limit number of results",
        )

    def handle(self, *args, **options):
        """
        Handle the list_projects command.

        Args:
            *args: Positional arguments (unused)
            **options: Command options

        Raises:
            CommandError: If listing fails
        """
        import json
        from datetime import datetime

        try:
            # Get all projects with job counts
            projects = models.Project.objects.annotate(
                job_count=Count('jobs')
            )

            # Apply name filter if specified
            name_filter = options.get("filter")
            if name_filter:
                # Get all projects and filter with fnmatch
                all_projects = list(projects)
                projects = [p for p in all_projects if fnmatch.fnmatch(p.name, name_filter)]
            else:
                projects = list(projects)

            # Sort projects
            sort_by = options["sort_by"]
            reverse = options["reverse"]

            if sort_by == "name":
                projects.sort(key=lambda p: p.name, reverse=reverse)
            elif sort_by == "creation_time":
                projects.sort(key=lambda p: p.creation_time, reverse=reverse)
            elif sort_by == "last_access":
                projects.sort(key=lambda p: p.last_access, reverse=reverse)
            elif sort_by == "jobs":
                projects.sort(key=lambda p: p.job_count, reverse=reverse)

            # Apply limit if specified
            limit = options.get("limit")
            if limit:
                projects = projects[:limit]

            # Handle empty results
            if not projects:
                if name_filter:
                    self.stdout.write(f"No projects found matching '{name_filter}'")
                else:
                    self.stdout.write("No projects found")
                return

            # Output based on format
            json_output = options["json"]
            names_only = options["names_only"]

            if json_output:
                # JSON output
                output = []
                for project in projects:
                    output.append({
                        "uuid": str(project.uuid),
                        "name": project.name,
                        "description": project.description,
                        "directory": project.directory,
                        "creation_time": project.creation_time.isoformat() if project.creation_time else None,
                        "creation_user": project.creation_user,
                        "creation_host": project.creation_host,
                        "last_access": project.last_access.isoformat() if project.last_access else None,
                        "job_count": project.job_count,
                        "last_job_number": project.last_job_number,
                    })
                self.stdout.write(json.dumps(output, indent=2))

            elif names_only:
                # Names only output (one per line)
                for project in projects:
                    self.stdout.write(project.name)

            else:
                # Verbose table output
                self.stdout.write(self.style.SUCCESS(f"\n{'='*130}"))
                self.stdout.write(self.style.SUCCESS(f"CCP4i2 Projects ({len(projects)} found)"))
                self.stdout.write(self.style.SUCCESS(f"{'='*130}\n"))

                # Table header
                self.stdout.write(
                    f"{'ID':<6} {'UUID':<38} {'Name':<20} {'Jobs':<6} {'Created':<20} {'Description':<36}"
                )
                self.stdout.write("-" * 130)

                # Table rows
                for project in projects:
                    # Truncate/format fields
                    project_id = str(project.id)
                    uuid_full = str(project.uuid)
                    name = project.name[:19] if len(project.name) > 19 else project.name
                    description = project.description[:35] if len(project.description) > 35 else project.description
                    created = project.creation_time.strftime("%Y-%m-%d %H:%M:%S") if project.creation_time else "Unknown"
                    job_count = str(project.job_count)

                    self.stdout.write(
                        f"{project_id:<6} {uuid_full:<38} {name:<20} {job_count:<6} {created:<20} {description:<36}"
                    )

                # Summary footer
                self.stdout.write(f"\n{'='*130}")
                total_jobs = sum(p.job_count for p in projects)
                self.stdout.write(f"Total: {len(projects)} projects, {total_jobs} jobs")
                self.stdout.write(f"{'='*130}\n")

                # Show sort/filter info if applied
                if name_filter:
                    self.stdout.write(f"Filtered by: {name_filter}")
                if sort_by != "creation_time" or reverse:
                    sort_desc = f"Sorted by: {sort_by}"
                    if reverse:
                        sort_desc += " (reversed)"
                    self.stdout.write(sort_desc)

        except Exception as e:
            raise CommandError(f"Failed to list projects: {str(e)}")
