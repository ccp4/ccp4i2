"""
Management command to create database backups and upload to Azure Blob Storage.

This command creates JSON fixture dumps of all app data and optionally
uploads them to Azure Blob Storage for archival.

Usage:
    # Create backups locally
    python manage.py backup_database --output-dir /tmp/backups

    # Create backups and upload to Azure Blob Storage
    python manage.py backup_database --upload

    # Backup specific apps only
    python manage.py backup_database --apps registry assays constructs

    # List available backups in blob storage
    python manage.py backup_database --list
"""

import json
import os
import tempfile
from datetime import datetime
from pathlib import Path

from django.apps import apps
from django.core.management import call_command
from django.core.management.base import BaseCommand, CommandError

# Try to import Azure SDK (optional - only needed for upload)
try:
    from azure.storage.blob import BlobServiceClient
    from azure.identity import DefaultAzureCredential
    AZURE_SDK_AVAILABLE = True
except ImportError:
    AZURE_SDK_AVAILABLE = False


# Apps to backup and their Django app labels
BACKUP_APPS = {
    'ccp4i2': {
        'apps': ['ccp4i2.db'],
        'description': 'CCP4i2 projects, jobs, and files',
    },
    'registry': {
        'apps': ['compounds.registry'],
        'description': 'Compound registry data',
    },
    'assays': {
        'apps': ['compounds.assays'],
        'description': 'Assay experiments and results',
    },
    'constructs': {
        'apps': ['compounds.constructs'],
        'description': 'Construct/plasmid database',
    },
    'users': {
        'apps': ['users', 'auth'],
        'description': 'User accounts and profiles',
    },
}


class Command(BaseCommand):
    """Create database backups and upload to Azure Blob Storage."""

    help = "Create database backups and optionally upload to Azure Blob Storage"

    def add_arguments(self, parser):
        parser.add_argument(
            '--apps',
            nargs='+',
            choices=list(BACKUP_APPS.keys()),
            default=list(BACKUP_APPS.keys()),
            help='Apps to backup (default: all)',
        )
        parser.add_argument(
            '--output-dir',
            type=str,
            help='Directory to save backup files (default: temp directory)',
        )
        parser.add_argument(
            '--upload',
            action='store_true',
            help='Upload backups to Azure Blob Storage',
        )
        parser.add_argument(
            '--container',
            type=str,
            default='django-uploads',
            help='Azure Blob container name (default: django-uploads)',
        )
        parser.add_argument(
            '--prefix',
            type=str,
            default='fixtures',
            help='Blob prefix/folder for backups (default: fixtures)',
        )
        parser.add_argument(
            '--list',
            action='store_true',
            dest='list_backups',
            help='List available backups in blob storage',
        )
        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Show detailed progress',
        )

    def handle(self, *args, **options):
        if options['list_backups']:
            return self._list_backups(options)

        apps_to_backup = options['apps']
        output_dir = options.get('output_dir')
        upload = options['upload']
        verbose = options['verbose']

        self.verbose = verbose
        self.stdout.write("\nDatabase Backup")
        self.stdout.write("=" * 60)

        # Create timestamp for this backup run
        timestamp = datetime.utcnow().strftime('%Y%m%d-%H-%M')
        self.stdout.write(f"Timestamp: {timestamp}")

        # Setup output directory
        if output_dir:
            backup_dir = Path(output_dir)
            backup_dir.mkdir(parents=True, exist_ok=True)
            cleanup = False
        else:
            backup_dir = Path(tempfile.mkdtemp(prefix='ccp4i2_backup_'))
            cleanup = not upload  # Keep temp files if uploading

        self.stdout.write(f"Output directory: {backup_dir}")
        self.stdout.write("")

        backup_files = []

        # Create backups for each app
        for app_name in apps_to_backup:
            app_config = BACKUP_APPS[app_name]
            self.stdout.write(f"Backing up: {app_name} ({app_config['description']})")

            filename = f"{timestamp}-{app_name}.json"
            filepath = backup_dir / filename

            try:
                # Get the Django app labels
                django_apps = app_config['apps']

                # Filter to only installed apps
                installed_apps = []
                for app_label in django_apps:
                    try:
                        apps.get_app_config(app_label.split('.')[-1])
                        installed_apps.append(app_label)
                    except LookupError:
                        if verbose:
                            self.stdout.write(f"  Skipping {app_label} (not installed)")

                if not installed_apps:
                    self.stdout.write(self.style.WARNING(f"  No apps to backup for {app_name}"))
                    continue

                # Run dumpdata
                with open(filepath, 'w') as f:
                    call_command(
                        'dumpdata',
                        *installed_apps,
                        indent=2,
                        stdout=f,
                        verbosity=0,
                    )

                # Get file size
                size = filepath.stat().st_size
                size_str = self._format_size(size)

                self.stdout.write(self.style.SUCCESS(f"  Created: {filename} ({size_str})"))
                backup_files.append((filename, filepath, size))

            except Exception as e:
                self.stderr.write(self.style.ERROR(f"  Error: {e}"))

        self.stdout.write("")

        # Upload to Azure if requested
        if upload and backup_files:
            self._upload_backups(backup_files, options)

        # Cleanup temp directory
        if cleanup and not output_dir:
            import shutil
            shutil.rmtree(backup_dir, ignore_errors=True)
            self.stdout.write(f"Cleaned up temp directory")

        # Summary
        self.stdout.write("")
        self.stdout.write("=" * 60)
        total_size = sum(f[2] for f in backup_files)
        self.stdout.write(f"Total: {len(backup_files)} backups, {self._format_size(total_size)}")
        self.stdout.write("=" * 60)

    def _upload_backups(self, backup_files, options):
        """Upload backup files to Azure Blob Storage."""
        if not AZURE_SDK_AVAILABLE:
            self.stderr.write(self.style.ERROR(
                "Azure SDK not available. Install with: pip install azure-storage-blob azure-identity"
            ))
            return

        container = options['container']
        prefix = options['prefix']

        # Get storage account from environment
        account_name = os.environ.get('AZURE_STORAGE_ACCOUNT_NAME')
        if not account_name:
            self.stderr.write(self.style.ERROR(
                "AZURE_STORAGE_ACCOUNT_NAME environment variable not set"
            ))
            return

        self.stdout.write(f"Uploading to Azure Blob Storage...")
        self.stdout.write(f"  Account: {account_name}")
        self.stdout.write(f"  Container: {container}")
        self.stdout.write(f"  Prefix: {prefix}/")

        try:
            # Use DefaultAzureCredential (Managed Identity in production)
            credential = DefaultAzureCredential()
            blob_service = BlobServiceClient(
                account_url=f"https://{account_name}.blob.core.windows.net",
                credential=credential
            )
            container_client = blob_service.get_container_client(container)

            for filename, filepath, size in backup_files:
                blob_name = f"{prefix}/{filename}"
                blob_client = container_client.get_blob_client(blob_name)

                with open(filepath, 'rb') as f:
                    blob_client.upload_blob(f, overwrite=True)

                self.stdout.write(self.style.SUCCESS(f"  Uploaded: {blob_name}"))

        except Exception as e:
            self.stderr.write(self.style.ERROR(f"Upload failed: {e}"))

    def _list_backups(self, options):
        """List available backups in Azure Blob Storage."""
        if not AZURE_SDK_AVAILABLE:
            self.stderr.write(self.style.ERROR(
                "Azure SDK not available. Install with: pip install azure-storage-blob azure-identity"
            ))
            return

        container = options['container']
        prefix = options['prefix']

        account_name = os.environ.get('AZURE_STORAGE_ACCOUNT_NAME')
        if not account_name:
            self.stderr.write(self.style.ERROR(
                "AZURE_STORAGE_ACCOUNT_NAME environment variable not set"
            ))
            return

        self.stdout.write(f"\nBackups in Azure Blob Storage")
        self.stdout.write("=" * 60)
        self.stdout.write(f"Account: {account_name}")
        self.stdout.write(f"Container: {container}")
        self.stdout.write(f"Prefix: {prefix}/")
        self.stdout.write("")

        try:
            credential = DefaultAzureCredential()
            blob_service = BlobServiceClient(
                account_url=f"https://{account_name}.blob.core.windows.net",
                credential=credential
            )
            container_client = blob_service.get_container_client(container)

            blobs = list(container_client.list_blobs(name_starts_with=f"{prefix}/"))

            if not blobs:
                self.stdout.write("No backups found.")
                return

            # Group by date
            by_date = {}
            for blob in blobs:
                name = blob.name.replace(f"{prefix}/", "")
                # Extract date from filename (YYYYMMDD-HH-MM-app.json)
                parts = name.split('-')
                if len(parts) >= 3:
                    date_key = parts[0]  # YYYYMMDD
                    if date_key not in by_date:
                        by_date[date_key] = []
                    by_date[date_key].append((name, blob.size, blob.last_modified))

            # Display sorted by date (newest first)
            for date_key in sorted(by_date.keys(), reverse=True)[:10]:
                files = by_date[date_key]
                date_str = f"{date_key[:4]}-{date_key[4:6]}-{date_key[6:8]}"
                self.stdout.write(f"\n{date_str}:")
                for name, size, modified in sorted(files):
                    self.stdout.write(f"  {name} ({self._format_size(size)})")

            self.stdout.write(f"\nTotal: {len(blobs)} backup files")

        except Exception as e:
            self.stderr.write(self.style.ERROR(f"Failed to list backups: {e}"))

    def _format_size(self, size):
        """Format file size in human-readable format."""
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size < 1024:
                return f"{size:.1f} {unit}"
            size /= 1024
        return f"{size:.1f} TB"
