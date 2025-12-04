"""
Django management command to list files in the database.

Usage:
    python manage.py list_files --project <name>
    python manage.py list_files --job <uuid>
    python manage.py list_files --all
"""

import json
from django.core.management.base import BaseCommand
from ccp4x.db.models import File, Project, Job


class Command(BaseCommand):
    """List files registered in the database."""

    help = "List files in the database"

    def add_arguments(self, parser):
        """Add command-line arguments."""
        parser.add_argument(
            '--project',
            type=str,
            help='Project name or UUID'
        )
        parser.add_argument(
            '--job',
            type=str,
            help='Job number or UUID'
        )
        parser.add_argument(
            '--all',
            action='store_true',
            help='List all files across all projects'
        )
        parser.add_argument(
            '--format',
            type=str,
            choices=['table', 'csv', 'json'],
            default='table',
            help='Output format'
        )
        parser.add_argument(
            '--filter',
            type=str,
            help='Filter files by name (substring match)'
        )

    def handle(self, *args, **options):
        """Execute the list operation."""
        format_type = options.get('format', 'table')
        filter_name = options.get('filter')

        # Determine scope
        if options.get('job'):
            files = self._get_job_files(options['job'])
            scope = f"job {options['job']}"
        elif options.get('project'):
            files = self._get_project_files(options['project'])
            scope = f"project {options['project']}"
        elif options.get('all'):
            files = File.objects.all()
            scope = "all projects"
        else:
            self.stdout.write(self.style.ERROR(
                "Error: Must specify --project, --job, or --all"
            ))
            return

        # Apply filter
        if filter_name:
            files = files.filter(name__icontains=filter_name)

        # Order by name (File model doesn't have created_at)
        files = files.order_by('name')

        # Output
        if format_type == 'table':
            self._output_table(files, scope)
        elif format_type == 'csv':
            self._output_csv(files)
        elif format_type == 'json':
            self._output_json(files)

    def _get_project_files(self, project_identifier):
        """Get files for a specific project."""
        try:
            # Try UUID first
            import uuid
            project_uuid = uuid.UUID(project_identifier)
            project = Project.objects.get(uuid=project_uuid)
        except (ValueError, Project.DoesNotExist):
            # Try name
            project = Project.objects.get(name=project_identifier)

        # Get all jobs in project, then get files for those jobs
        jobs = Job.objects.filter(project=project)
        return File.objects.filter(job__in=jobs)

    def _get_job_files(self, job_identifier):
        """Get files used by a specific job."""
        from ccp4x.db.models import FileUse

        try:
            # Try UUID first
            import uuid
            job_uuid = uuid.UUID(job_identifier)
            job = Job.objects.get(uuid=job_uuid)
        except (ValueError, Job.DoesNotExist):
            # Try job number (need project context)
            self.stdout.write(self.style.ERROR(
                "Error: Job must be specified by UUID"
            ))
            return File.objects.none()

        # Get all files used by this job
        file_uses = FileUse.objects.filter(job=job)
        file_ids = [fu.file.id for fu in file_uses]
        return File.objects.filter(id__in=file_ids)

    def _output_table(self, files, scope):
        """Output files as a formatted table."""
        if files.count() == 0:
            self.stdout.write(self.style.WARNING(f"No files found for {scope}"))
            return

        # Column widths
        name_width = min(max(len(f.name) for f in files), 60)
        name_width = max(name_width, len('File Name'))

        uuid_width = 36  # Standard UUID width
        project_width = 20

        # Header
        header = (
            f"{'File Name':<{name_width}}  "
            f"{'UUID':<{uuid_width}}  "
            f"{'Project':<{project_width}}  "
            f"{'Type':<20}  "
            f"{'Job':<10}"
        )
        separator = "=" * len(header)

        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(f"Files in {scope}"))
        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(header))
        self.stdout.write(self.style.SUCCESS(separator))

        # Rows
        for f in files:
            name = f.name[:name_width] if len(f.name) > name_width else f.name
            # Files are associated with jobs, not directly with projects
            project_name = f.job.project.name if f.job and f.job.project else "N/A"
            project_name = project_name[:project_width] if len(project_name) > project_width else project_name

            # Get file type name safely
            file_type = "Unknown"
            if f.type:
                file_type = f.type.name[:20]

            # Job number
            job_num = f.job.number if f.job else "N/A"

            row = (
                f"{name:<{name_width}}  "
                f"{str(f.uuid):<{uuid_width}}  "
                f"{project_name:<{project_width}}  "
                f"{file_type:<20}  "
                f"{job_num:<10}"
            )
            self.stdout.write(row)

        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(f"Total: {files.count()} files"))

    def _output_csv(self, files):
        """Output files as CSV."""
        import csv
        import sys

        writer = csv.writer(sys.stdout)
        writer.writerow(['Name', 'UUID', 'Job', 'Project', 'Type'])

        for f in files:
            job_num = f.job.number if f.job else "N/A"
            project_name = f.job.project.name if f.job and f.job.project else "N/A"
            file_type = f.type.name if f.type else "Unknown"

            writer.writerow([
                f.name,
                str(f.uuid),
                job_num,
                project_name,
                file_type
            ])

    def _output_json(self, files):
        """Output files as JSON."""
        data = [
            {
                'name': f.name,
                'uuid': str(f.uuid),
                'job': {
                    'number': f.job.number if f.job else None,
                    'uuid': str(f.job.uuid) if f.job else None,
                },
                'project': {
                    'name': f.job.project.name if f.job and f.job.project else None,
                    'uuid': str(f.job.project.uuid) if f.job and f.job.project else None,
                },
                'type': f.type.name if f.type else None,
            }
            for f in files
        ]

        self.stdout.write(json.dumps(data, indent=2))
