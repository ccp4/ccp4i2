"""
Django management command to show directory tree for a project.

Usage:
    python manage.py tree_project <project_name_or_uuid>
    python manage.py tree_project toxd --depth 3
    python manage.py tree_project toxd --no-sizes
"""

import uuid
from django.core.management.base import BaseCommand, CommandError
from ccp4i2.db.models import Project
from ccp4i2.lib.utils.directory_tree import visualize_project_directory


class Command(BaseCommand):
    """Show directory tree structure for a project."""

    help = "Display directory tree structure for a CCP4i2 project"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            'project',
            help='Project name or UUID'
        )
        parser.add_argument(
            '--depth',
            type=int,
            default=3,
            help='Maximum depth to display (default: 3, 0 = unlimited)'
        )
        parser.add_argument(
            '--no-sizes',
            action='store_true',
            help='Do not show file sizes'
        )
        parser.add_argument(
            '--show-hidden',
            action='store_true',
            help='Show hidden files and directories'
        )

    def handle(self, *args, **options):
        project_identifier = options['project']

        # Find project
        try:
            project = self.get_project(project_identifier)
        except Project.DoesNotExist:
            raise CommandError(
                f"Project not found: '{project_identifier}'. "
                "Use 'ccp4i2 projects list' to see available projects."
            )

        # Get options
        max_depth = options['depth'] if options['depth'] > 0 else None
        show_sizes = not options['no_sizes']
        show_hidden = options['show_hidden']

        # Import here to avoid circular dependency
        from ccp4i2.lib.utils.directory_tree import DirectoryTree
        from pathlib import Path

        # Generate tree
        tree = DirectoryTree(
            max_depth=max_depth,
            show_sizes=show_sizes,
            show_hidden=show_hidden
        )

        title = f"Project: {project.name} (UUID: {str(project.uuid)[:8]}...)"
        output = tree.visualize(Path(project.directory), title=title)

        # Print tree
        self.stdout.write(output)

        # Additional project info
        self.stdout.write(f"\nProject Directory: {project.directory}")
        self.stdout.write(f"Created: {project.creation_time}")
        if hasattr(project, 'job_count'):
            from django.db.models import Count
            job_count = Project.objects.filter(pk=project.pk).annotate(
                jobs=Count('jobs')
            ).first().jobs
            self.stdout.write(f"Jobs: {job_count}")

    def get_project(self, identifier: str) -> Project:
        """Get project by name or UUID."""
        # Try UUID first
        try:
            uuid.UUID(identifier)
            return Project.objects.get(uuid=identifier)
        except (ValueError, Project.DoesNotExist):
            pass

        # Try name
        try:
            return Project.objects.get(name=identifier)
        except Project.DoesNotExist:
            pass

        raise Project.DoesNotExist()
