"""
Management command to clear all construct data for reimport.

This command deletes all data from the constructs app tables, allowing
a fresh reimport from the legacy database.

Usage:
    python manage.py clear_constructs

    # Confirm without prompt:
    python manage.py clear_constructs --yes

    # Also delete media files:
    python manage.py clear_constructs --yes --delete-files
"""

import shutil
from pathlib import Path

from django.conf import settings
from django.core.management.base import BaseCommand
from django.db import connection


class Command(BaseCommand):
    help = 'Clear all construct data for reimport'

    def add_arguments(self, parser):
        parser.add_argument(
            '--yes',
            action='store_true',
            help='Skip confirmation prompt'
        )
        parser.add_argument(
            '--delete-files',
            action='store_true',
            help='Also delete associated media files'
        )

    def handle(self, *args, **options):
        confirm = options['yes']
        delete_files = options['delete_files']

        # Tables to clear in dependency order (children before parents)
        tables = [
            'constructs_expressiontag',
            'constructs_sequencingresult',
            'constructs_cassetteuse',
            'constructs_cassette',
            'constructs_proteinsynonym',
            'constructs_proteinuse',
            'constructs_protein',
            'constructs_plasmid',
            'constructs_constructproject',
            # Reference data - uncomment if you want to clear these too
            # 'constructs_expressiontaglocation',
            # 'constructs_protease',
            # 'constructs_expressiontagtype',
        ]

        if not confirm:
            self.stdout.write(self.style.WARNING(
                '\nThis will DELETE all data from the following tables:'
            ))
            for table in tables:
                self.stdout.write(f'  - {table}')

            if delete_files:
                self.stdout.write(self.style.WARNING(
                    '\nAND delete files from: compounds/constructs/'
                ))

            response = input('\nAre you sure? (yes/no): ')
            if response.lower() != 'yes':
                self.stdout.write(self.style.ERROR('Aborted.'))
                return

        # Get counts before deletion
        counts = {}
        with connection.cursor() as cursor:
            for table in tables:
                try:
                    cursor.execute(f'SELECT COUNT(*) FROM {table}')
                    counts[table] = cursor.fetchone()[0]
                except Exception:
                    counts[table] = 0

        self.stdout.write('Clearing construct tables...')

        # Delete in dependency order
        with connection.cursor() as cursor:
            for table in tables:
                try:
                    cursor.execute(f'DELETE FROM {table}')
                    self.stdout.write(f'  Deleted {counts[table]} rows from {table}')
                except Exception as e:
                    self.stdout.write(
                        self.style.WARNING(f'  Skipped {table}: {e}')
                    )

        # Optionally delete media files
        if delete_files:
            media_dir = Path(settings.MEDIA_ROOT) / 'compounds' / 'constructs'
            if media_dir.exists():
                file_count = sum(1 for _ in media_dir.rglob('*') if _.is_file())
                shutil.rmtree(media_dir)
                self.stdout.write(f'  Deleted {file_count} files from {media_dir}')
            else:
                self.stdout.write(f'  No files to delete (directory does not exist)')

        self.stdout.write(self.style.SUCCESS('\nConstruct data cleared successfully.'))
        self.stdout.write('You can now reimport with: python manage.py migrate_legacy_constructs')
