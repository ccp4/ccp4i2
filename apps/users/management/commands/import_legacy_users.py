"""
Management command to import users from legacy auth fixtures.

This command imports auth.user records from legacy CCP4i2Docker fixtures,
creating Django users with associated UserProfile records.

Usage:
    python manage.py import_legacy_users auth.json
    python manage.py import_legacy_users auth.json --dry-run --verbose
"""

import json
import re
from datetime import datetime

from django.contrib.auth import get_user_model
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.utils import timezone

from users.models import UserProfile

User = get_user_model()


class Command(BaseCommand):
    """Import users from legacy auth fixtures."""

    help = "Import users from legacy auth.json fixtures"

    def add_arguments(self, parser):
        parser.add_argument(
            'fixture_file',
            type=str,
            help='Path to auth.json fixture file',
        )
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Validate without committing changes',
        )
        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Show detailed progress',
        )
        parser.add_argument(
            '--update-existing',
            action='store_true',
            help='Update existing users instead of skipping',
        )

    def handle(self, *args, **options):
        fixture_file = options['fixture_file']
        dry_run = options['dry_run']
        verbose = options['verbose']
        update_existing = options['update_existing']

        try:
            with open(fixture_file, 'r') as f:
                content = f.read()
                # Handle Django debug output in fixture
                start = content.find('[')
                if start > 0:
                    content = content[start:]
                data = json.loads(content)
        except FileNotFoundError:
            raise CommandError(f"Fixture file not found: {fixture_file}")
        except json.JSONDecodeError as e:
            raise CommandError(f"Invalid JSON: {e}")

        self.stdout.write(f"\nImporting users from: {fixture_file}")
        self.stdout.write("-" * 60)

        if dry_run:
            self.stdout.write(self.style.WARNING("[DRY RUN - no changes will be saved]"))

        # Filter to auth.user records
        user_records = [r for r in data if r.get('model') == 'auth.user']
        self.stdout.write(f"Found {len(user_records)} user records")

        stats = {'created': 0, 'updated': 0, 'skipped': 0, 'errors': 0}

        try:
            with transaction.atomic():
                for record in user_records:
                    try:
                        result = self._import_user(
                            record, update_existing, verbose, dry_run
                        )
                        stats[result] += 1
                    except Exception as e:
                        self.stderr.write(
                            self.style.ERROR(f"  Error: {e}")
                        )
                        stats['errors'] += 1

                if dry_run:
                    raise DryRunRollback()

        except DryRunRollback:
            self.stdout.write(self.style.WARNING("\n[DRY RUN - all changes rolled back]"))

        # Summary
        self.stdout.write("\n" + "-" * 60)
        self.stdout.write(f"Created: {stats['created']}")
        self.stdout.write(f"Updated: {stats['updated']}")
        self.stdout.write(f"Skipped: {stats['skipped']}")
        if stats['errors']:
            self.stdout.write(self.style.ERROR(f"Errors: {stats['errors']}"))
        self.stdout.write("-" * 60)

    def _import_user(self, record, update_existing, verbose, dry_run):
        """Import a single user record."""
        fields = record.get('fields', {})
        username = fields.get('username', '')
        email = fields.get('email', '')

        if not username:
            if verbose:
                self.stdout.write(f"  Skipping record with no username")
            return 'skipped'

        # Normalize email
        if not email and '@' in username:
            email = username

        # Check for existing user
        existing = User.objects.filter(username=username).first()
        if not existing and email:
            existing = User.objects.filter(email__iexact=email).first()

        if existing and not update_existing:
            if verbose:
                self.stdout.write(f"  Skipping existing: {username}")
            return 'skipped'

        # Prepare user data
        user_data = {
            'email': email,
            'first_name': fields.get('first_name', '')[:30],
            'last_name': fields.get('last_name', '')[:150],
            'is_active': fields.get('is_active', True),
            'is_staff': fields.get('is_staff', False),
            'is_superuser': fields.get('is_superuser', False),
        }

        # Handle date_joined
        date_joined = fields.get('date_joined')
        if date_joined:
            user_data['date_joined'] = self._parse_datetime(date_joined)

        if existing:
            # Update existing user
            for key, value in user_data.items():
                setattr(existing, key, value)
            existing.save()
            user = existing
            result = 'updated'
            if verbose:
                self.stdout.write(f"  Updated: {username}")
        else:
            # Create new user
            user = User.objects.create(
                username=username,
                **user_data
            )
            result = 'created'
            if verbose:
                self.stdout.write(self.style.SUCCESS(f"  Created: {username}"))

        # Update profile
        profile, _ = UserProfile.objects.get_or_create(user=user)
        profile.legacy_username = username
        profile.legacy_display_name = f"{fields.get('first_name', '')} {fields.get('last_name', '')}".strip()
        profile.imported_at = timezone.now()
        profile.save()

        return result

    def _parse_datetime(self, value):
        """Parse datetime string to timezone-aware datetime."""
        if not value:
            return None

        # Add timezone if missing
        if not value.endswith('Z') and '+' not in value[-6:]:
            value = value + 'Z'

        try:
            # Try ISO format
            if value.endswith('Z'):
                value = value[:-1] + '+00:00'
            dt = datetime.fromisoformat(value)
            if dt.tzinfo is None:
                dt = timezone.make_aware(dt)
            return dt
        except (ValueError, TypeError):
            return None


class DryRunRollback(Exception):
    """Exception to trigger rollback in dry-run mode."""
    pass
