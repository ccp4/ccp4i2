"""
Django management command to show directory tree for a job.

Usage:
    python manage.py tree_job <project> <job>
    python manage.py tree_job toxd 1 --depth 2
    python manage.py tree_job toxd 1 --no-sizes
"""

import uuid
from django.core.management.base import BaseCommand, CommandError
from ccp4i2.db.models import Job, Project
from ccp4i2.lib.utils.directory_tree import visualize_job_directory


class Command(BaseCommand):
    """Show directory tree structure for a job."""

    help = "Display directory tree structure for a CCP4i2 job"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            'project',
            help='Project name or UUID'
        )
        parser.add_argument(
            'job',
            help='Job number or UUID'
        )
        parser.add_argument(
            '--depth',
            type=int,
            default=2,
            help='Maximum depth to display (default: 2, 0 = unlimited)'
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
        job_identifier = options['job']

        # Find job
        try:
            job = self.get_job(project_identifier, job_identifier)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            raise CommandError(str(e))

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

        title = (
            f"Job #{job.number}: {job.title}\n"
            f"Task: {job.task_name} | Status: {job.get_status_display()} | "
            f"Project: {job.project.name}"
        )
        output = tree.visualize(Path(job.directory), title=title)

        # Print tree
        self.stdout.write(output)

        # Additional job info
        self.stdout.write(f"\nJob Directory: {job.directory}")
        self.stdout.write(f"Created: {job.creation_time}")
        if job.finish_time:
            self.stdout.write(f"Finished: {job.finish_time}")

    def get_job(self, project_identifier: str, job_identifier: str) -> Job:
        """Get job by project and job identifiers."""
        # First, get the project
        project = None

        # Try project UUID
        try:
            uuid.UUID(project_identifier)
            project = Project.objects.get(uuid=project_identifier)
        except (ValueError, Project.DoesNotExist):
            # Try project name
            try:
                project = Project.objects.get(name=project_identifier)
            except Project.DoesNotExist:
                raise Project.DoesNotExist(
                    f"Project not found: '{project_identifier}'"
                )

        # Now get the job
        # Try job UUID
        try:
            uuid.UUID(job_identifier)
            job = Job.objects.get(uuid=job_identifier, project=project)
            return job
        except (ValueError, Job.DoesNotExist):
            pass

        # Try job number
        try:
            job = Job.objects.get(number=job_identifier, project=project)
            return job
        except Job.DoesNotExist:
            raise Job.DoesNotExist(
                f"Job not found: '{job_identifier}' in project '{project.name}'"
            )
