"""
Management command to import aliases from legacy RegisterCompounds fixture.

The legacy database stored aliases in a nested format:
    {"data": ["alias1", "alias2"]}

This command extracts those aliases and maps them to current compounds
by reg_number (reg_id in legacy).

Usage:
    # Dry run (preview what would be imported)
    python manage.py import_legacy_aliases /path/to/fixture.json --dry-run

    # Actually import
    python manage.py import_legacy_aliases /path/to/fixture.json

    # Merge with existing aliases (don't overwrite)
    python manage.py import_legacy_aliases /path/to/fixture.json --merge
"""

import json
import logging

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from compounds.registry.models import Compound

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'Import compound aliases from legacy RegisterCompounds fixture file'

    def add_arguments(self, parser):
        parser.add_argument(
            'fixture_file',
            help='Path to the legacy RegisterCompounds JSON fixture file'
        )
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Preview what would be imported without making changes'
        )
        parser.add_argument(
            '--merge',
            action='store_true',
            help='Merge with existing aliases instead of overwriting'
        )
        parser.add_argument(
            '--skip-lines',
            type=int,
            default=3,
            help='Number of header lines to skip in fixture file (default: 3)'
        )

    def handle(self, *args, **options):
        fixture_path = options['fixture_file']
        dry_run = options['dry_run']
        merge = options['merge']
        skip_lines = options['skip_lines']

        self.stdout.write(f"Reading fixture file: {fixture_path}")

        # Read and parse the fixture file
        try:
            with open(fixture_path, 'r') as f:
                # Skip header lines (e.g., "In CCP4i2AppConfig.ready\nccp4i2 version...")
                for _ in range(skip_lines):
                    f.readline()
                data = json.load(f)
        except FileNotFoundError:
            raise CommandError(f"Fixture file not found: {fixture_path}")
        except json.JSONDecodeError as e:
            raise CommandError(f"Invalid JSON in fixture file: {e}")

        # Extract compounds with aliases from the fixture
        legacy_aliases = {}
        for record in data:
            if record.get('model') != 'RegisterCompounds.regdata':
                continue

            fields = record.get('fields', {})
            reg_id = fields.get('reg_id')
            aliases_data = fields.get('aliases')

            if not reg_id:
                continue

            # Extract aliases from nested format {"data": [...]}
            if aliases_data and isinstance(aliases_data, dict):
                alias_list = aliases_data.get('data', [])
                if alias_list and isinstance(alias_list, list) and len(alias_list) > 0:
                    legacy_aliases[reg_id] = alias_list

        self.stdout.write(f"Found {len(legacy_aliases)} compounds with aliases in fixture")

        if not legacy_aliases:
            self.stdout.write(self.style.WARNING("No aliases to import"))
            return

        # Match with current compounds and prepare updates
        updates = []
        not_found = []
        already_set = []

        for reg_number, aliases in legacy_aliases.items():
            try:
                compound = Compound.objects.get(reg_number=reg_number)
            except Compound.DoesNotExist:
                not_found.append(reg_number)
                continue

            # Check if compound already has aliases
            current_aliases = compound.aliases or []

            if merge:
                # Merge: add new aliases without duplicates
                merged = list(current_aliases)
                for alias in aliases:
                    if alias not in merged:
                        merged.append(alias)
                if merged != current_aliases:
                    updates.append((compound, merged, current_aliases))
                else:
                    already_set.append(reg_number)
            else:
                # Overwrite mode
                if aliases != current_aliases:
                    updates.append((compound, aliases, current_aliases))
                else:
                    already_set.append(reg_number)

        # Report findings
        self.stdout.write(f"\nSummary:")
        self.stdout.write(f"  - Compounds to update: {len(updates)}")
        self.stdout.write(f"  - Already up-to-date: {len(already_set)}")
        self.stdout.write(f"  - Not found in database: {len(not_found)}")

        if not_found:
            self.stdout.write(f"\nCompounds not found (first 10): {not_found[:10]}")

        if dry_run:
            self.stdout.write(self.style.WARNING("\n=== DRY RUN - No changes will be saved ===\n"))

        # Show what will be updated
        if updates:
            self.stdout.write("\nUpdates to apply:")
            for compound, new_aliases, old_aliases in updates[:20]:
                old_str = str(old_aliases) if old_aliases else '[]'
                self.stdout.write(
                    f"  {compound.formatted_id}: {old_str} -> {new_aliases}"
                )
            if len(updates) > 20:
                self.stdout.write(f"  ... and {len(updates) - 20} more")

        # Apply updates
        if not dry_run and updates:
            with transaction.atomic():
                updated_count = 0
                for compound, new_aliases, _ in updates:
                    compound.aliases = new_aliases
                    compound.save(update_fields=['aliases'])
                    updated_count += 1

                self.stdout.write(
                    self.style.SUCCESS(f"\nSuccessfully updated {updated_count} compounds")
                )
        elif dry_run:
            self.stdout.write("\nRun without --dry-run to apply changes.")
