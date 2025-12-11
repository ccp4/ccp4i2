"""
Django management command to list all file types in the database.

Usage:
    python manage.py list_filetypes
    ccp4i2 filetypes list
"""

from django.core.management.base import BaseCommand
from ccp4i2.db.models import FileType


class Command(BaseCommand):
    """List all file types registered in the database."""

    help = "List all file types with their IDs, names, and descriptions"

    def add_arguments(self, parser):
        parser.add_argument(
            '--format',
            type=str,
            choices=['table', 'csv', 'json'],
            default='table',
            help='Output format (default: table)'
        )
        parser.add_argument(
            '--filter',
            type=str,
            help='Filter file types by name (case-insensitive substring match)'
        )

    def handle(self, *args, **options):
        output_format = options['format']
        name_filter = options.get('filter')

        # Query file types
        filetypes = FileType.objects.all().order_by('name')

        if name_filter:
            filetypes = filetypes.filter(name__icontains=name_filter)

        if not filetypes.exists():
            if name_filter:
                self.stdout.write(self.style.WARNING(
                    f"No file types found matching '{name_filter}'"
                ))
            else:
                self.stdout.write(self.style.WARNING(
                    "No file types found in database"
                ))
            return

        if output_format == 'table':
            self._output_table(filetypes)
        elif output_format == 'csv':
            self._output_csv(filetypes)
        elif output_format == 'json':
            self._output_json(filetypes)

    def _output_table(self, filetypes):
        """Output file types as a formatted table."""
        # Calculate column widths
        name_width = max(len(ft.name) for ft in filetypes)
        name_width = max(name_width, len('MIME Type'))
        name_width = min(name_width, 50)  # Cap at 50 chars

        desc_width = 60  # Fixed width for description

        # Header
        header = f"{'MIME Type':<{name_width}}  {'Description':<{desc_width}}"
        separator = "=" * len(header)

        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(header))
        self.stdout.write(self.style.SUCCESS(separator))

        # Rows
        for ft in filetypes:
            ft_name = ft.name[:name_width] if len(ft.name) > name_width else ft.name
            ft_desc = ft.description[:desc_width] if ft.description and len(ft.description) > desc_width else (ft.description or '')

            row = f"{ft_name:<{name_width}}  {ft_desc:<{desc_width}}"
            self.stdout.write(row)

        self.stdout.write(self.style.SUCCESS(separator))
        self.stdout.write(self.style.SUCCESS(f"Total: {filetypes.count()} file types"))

    def _output_csv(self, filetypes):
        """Output file types as CSV."""
        import csv
        import sys

        writer = csv.writer(sys.stdout)
        writer.writerow(['MIME Type', 'Description'])

        for ft in filetypes:
            writer.writerow([
                ft.name,
                ft.description or ''
            ])

    def _output_json(self, filetypes):
        """Output file types as JSON."""
        import json

        data = [
            {
                'name': ft.name,
                'description': ft.description or ''
            }
            for ft in filetypes
        ]

        self.stdout.write(json.dumps(data, indent=2))
