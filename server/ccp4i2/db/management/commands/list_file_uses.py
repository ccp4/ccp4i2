"""
Django management command to list file uses (file-job relationships).

Usage:
    python manage.py list_file_uses --job <uuid>
    python manage.py list_file_uses --file <uuid>
    python manage.py list_file_uses --project <name>
"""

import json
from django.core.management.base import BaseCommand
from ccp4i2.db.models import FileUse, File, Job, Project


class Command(BaseCommand):
    """List file uses (which files are used by which jobs)."""

    help = "List file uses in the database"

    def add_arguments(self, parser):
        """Add command-line arguments."""
        parser.add_argument(
            '--job',
            type=str,
            help='Job UUID'
        )
        parser.add_argument(
            '--file',
            type=str,
            help='File UUID'
        )
        parser.add_argument(
            '--project',
            type=str,
            help='Project name or UUID'
        )
        parser.add_argument(
            '--format',
            type=str,
            choices=['table', 'csv', 'json'],
            default='table',
            help='Output format'
        )

    def handle(self, *args, **options):
        """Execute the list operation."""
        format_type = options.get('format', 'table')

        # Determine scope
        if options.get('job'):
            file_uses = self._get_job_file_uses(options['job'])
            scope = f"job {options['job']}"
        elif options.get('file'):
            file_uses = self._get_file_uses(options['file'])
            scope = f"file {options['file']}"
        elif options.get('project'):
            file_uses = self._get_project_file_uses(options['project'])
            scope = f"project {options['project']}"
        else:
            self.stdout.write(self.style.ERROR(
                "Error: Must specify --job, --file, or --project"
            ))
            return

        # Order by id (FileUse doesn't have created_at)
        file_uses = file_uses.order_by('-id')

        # Output
        if format_type == 'table':
            self._output_table(file_uses, scope)
        elif format_type == 'csv':
            self._output_csv(file_uses)
        elif format_type == 'json':
            self._output_json(file_uses)

    def _get_job_file_uses(self, job_identifier):
        """Get file uses for a specific job."""
        import uuid
        try:
            job_uuid = uuid.UUID(job_identifier)
            job = Job.objects.get(uuid=job_uuid)
        except (ValueError, Job.DoesNotExist) as e:
            self.stdout.write(self.style.ERROR(f"Job not found: {job_identifier}"))
            return FileUse.objects.none()

        return FileUse.objects.filter(job=job)

    def _get_file_uses(self, file_identifier):
        """Get uses of a specific file."""
        import uuid
        try:
            file_uuid = uuid.UUID(file_identifier)
            file_obj = File.objects.get(uuid=file_uuid)
        except (ValueError, File.DoesNotExist):
            self.stdout.write(self.style.ERROR(f"File not found: {file_identifier}"))
            return FileUse.objects.none()

        return FileUse.objects.filter(file=file_obj)

    def _get_project_file_uses(self, project_identifier):
        """Get all file uses in a project."""
        import uuid
        try:
            project_uuid = uuid.UUID(project_identifier)
            project = Project.objects.get(uuid=project_uuid)
        except (ValueError, Project.DoesNotExist):
            project = Project.objects.get(name=project_identifier)

        jobs = Job.objects.filter(project=project)
        return FileUse.objects.filter(job__in=jobs)

    def _output_table(self, file_uses, scope):
        """Output file uses as a formatted table."""
        if file_uses.count() == 0:
            self.stdout.write(self.style.WARNING(f"No file uses found for {scope}"))
            return

        # Column widths
        file_width = 40
        job_width = 20
        param_width = 20

        # Header
        header = (
            f"{'File Name':<{file_width}}  "
            f"{'Job':<{job_width}}  "
            f"{'Parameter':<{param_width}}  "
            f"{'Role':<15}  "
            f"{'Created':<20}"
        )
        separator = "=" * len(header)

        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(f"File Uses in {scope}"))
        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(header))
        self.stdout.write(self.style.SUCCESS(separator))

        # Rows
        for fu in file_uses:
            file_name = fu.file.name[:file_width] if len(fu.file.name) > file_width else fu.file.name
            task_name = str(fu.job.task_name)[:10] if fu.job.task_name else "N/A"
            job_desc = f"{fu.job.number} ({task_name})"
            job_desc = job_desc[:job_width]
            param_name = str(fu.job_param_name or "N/A")[:param_width]
            role = fu.get_role_display() if fu.role is not None else "Unknown"
            role = str(role)[:15]
            created = "N/A"  # FileUse doesn't have created_at

            row = (
                f"{file_name:<{file_width}}  "
                f"{job_desc:<{job_width}}  "
                f"{param_name:<{param_width}}  "
                f"{role:<15}  "
                f"{created:<20}"
            )
            self.stdout.write(row)

        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(f"Total: {file_uses.count()} file uses"))

    def _output_csv(self, file_uses):
        """Output file uses as CSV."""
        import csv
        import sys

        writer = csv.writer(sys.stdout)
        writer.writerow(['File UUID', 'File Name', 'Job UUID', 'Job Number', 'Parameter', 'Role'])

        for fu in file_uses:
            writer.writerow([
                str(fu.file.uuid),
                fu.file.name,
                str(fu.job.uuid),
                fu.job.number,
                fu.job_param_name or '',
                fu.get_role_display() if fu.role is not None else ''
            ])

    def _output_json(self, file_uses):
        """Output file uses as JSON."""
        data = [
            {
                'file': {
                    'uuid': str(fu.file.uuid),
                    'name': fu.file.name,
                },
                'job': {
                    'uuid': str(fu.job.uuid),
                    'number': fu.job.number,
                    'task_name': fu.job.task_name,
                },
                'param_name': fu.job_param_name,
                'role': fu.get_role_display() if fu.role is not None else None,
                'role_code': fu.role,
            }
            for fu in file_uses
        ]

        self.stdout.write(json.dumps(data, indent=2))
