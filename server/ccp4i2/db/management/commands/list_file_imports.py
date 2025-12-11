"""
Django management command to list file imports.

Usage:
    python manage.py list_file_imports --job <uuid>
    python manage.py list_file_imports --file <uuid>
    python manage.py list_file_imports --project <name>
"""

import json
from django.core.management.base import BaseCommand
from ccp4x.db.models import FileImport, File, Job, Project


class Command(BaseCommand):
    """List file imports (imported external files)."""

    help = "List file imports in the database"

    def add_arguments(self, parser):
        """Add command-line arguments."""
        parser.add_argument(
            '--job',
            type=str,
            help='Job UUID that imported files'
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
            imports = self._get_job_imports(options['job'])
            scope = f"job {options['job']}"
        elif options.get('file'):
            imports = self._get_file_imports(options['file'])
            scope = f"file {options['file']}"
        elif options.get('project'):
            imports = self._get_project_imports(options['project'])
            scope = f"project {options['project']}"
        else:
            self.stdout.write(self.style.ERROR(
                "Error: Must specify --job, --file, or --project"
            ))
            return

        # Order by time
        imports = imports.order_by('-time')

        # Output
        if format_type == 'table':
            self._output_table(imports, scope)
        elif format_type == 'csv':
            self._output_csv(imports)
        elif format_type == 'json':
            self._output_json(imports)

    def _get_job_imports(self, job_identifier):
        """Get imports for a specific job."""
        import uuid
        try:
            job_uuid = uuid.UUID(job_identifier)
            job = Job.objects.get(uuid=job_uuid)
        except (ValueError, Job.DoesNotExist):
            self.stdout.write(self.style.ERROR(f"Job not found: {job_identifier}"))
            return FileImport.objects.none()

        # Get files for this job, then get imports for those files
        files = File.objects.filter(job=job)
        return FileImport.objects.filter(file__in=files)

    def _get_file_imports(self, file_identifier):
        """Get imports of a specific file."""
        import uuid
        try:
            file_uuid = uuid.UUID(file_identifier)
            file_obj = File.objects.get(uuid=file_uuid)
        except (ValueError, File.DoesNotExist):
            self.stdout.write(self.style.ERROR(f"File not found: {file_identifier}"))
            return FileImport.objects.none()

        return FileImport.objects.filter(file=file_obj)

    def _get_project_imports(self, project_identifier):
        """Get all imports in a project."""
        import uuid
        try:
            project_uuid = uuid.UUID(project_identifier)
            project = Project.objects.get(uuid=project_uuid)
        except (ValueError, Project.DoesNotExist):
            project = Project.objects.get(name=project_identifier)

        # Get all jobs in project, then files for those jobs, then imports for those files
        jobs = Job.objects.filter(project=project)
        files = File.objects.filter(job__in=jobs)
        return FileImport.objects.filter(file__in=files)

    def _output_table(self, imports, scope):
        """Output imports as a formatted table."""
        if imports.count() == 0:
            self.stdout.write(self.style.WARNING(f"No file imports found for {scope}"))
            return

        # Column widths
        source_width = 50
        dest_width = 30
        checksum_width = 16  # First 8 chars + ...

        # Header
        header = (
            f"{'Source Path':<{source_width}}  "
            f"{'Imported As':<{dest_width}}  "
            f"{'Checksum':<{checksum_width}}  "
            f"{'Job':<10}  "
            f"{'Time':<20}"
        )
        separator = "=" * len(header)

        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(f"File Imports in {scope}"))
        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(header))
        self.stdout.write(self.style.SUCCESS(separator))

        # Rows
        for imp in imports:
            source = imp.name[:source_width] if len(imp.name) > source_width else imp.name
            dest = imp.file.name[:dest_width] if len(imp.file.name) > dest_width else imp.file.name
            checksum = (imp.checksum[:8] + "...") if imp.checksum else "N/A"
            # FileImport -> File -> Job
            job_num = imp.file.job.number if imp.file and imp.file.job else "N/A"
            time_str = imp.time.strftime("%Y-%m-%d %H:%M") if imp.time else "N/A"

            row = (
                f"{source:<{source_width}}  "
                f"{dest:<{dest_width}}  "
                f"{checksum:<{checksum_width}}  "
                f"{job_num:<10}  "
                f"{time_str:<20}"
            )
            self.stdout.write(row)

        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(f"Total: {imports.count()} imports"))

    def _output_csv(self, imports):
        """Output imports as CSV."""
        import csv
        import sys

        writer = csv.writer(sys.stdout)
        writer.writerow(['Source Path', 'File UUID', 'File Name', 'Job', 'Checksum', 'Time'])

        for imp in imports:
            job_num = imp.file.job.number if imp.file and imp.file.job else ''
            writer.writerow([
                imp.name,
                str(imp.file.uuid),
                imp.file.name,
                job_num,
                imp.checksum or '',
                imp.time.isoformat() if imp.time else ''
            ])

    def _output_json(self, imports):
        """Output imports as JSON."""
        data = [
            {
                'source_path': imp.name,
                'file': {
                    'uuid': str(imp.file.uuid),
                    'name': imp.file.name,
                },
                'job': {
                    'uuid': str(imp.file.job.uuid) if imp.file and imp.file.job else None,
                    'number': imp.file.job.number if imp.file and imp.file.job else None,
                },
                'checksum': imp.checksum,
                'time': imp.time.isoformat() if imp.time else None,
            }
            for imp in imports
        ]

        self.stdout.write(json.dumps(data, indent=2))
