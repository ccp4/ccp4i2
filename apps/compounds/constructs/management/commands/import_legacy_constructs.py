"""
Management command to import legacy ConstructDatabase fixtures.

This command transforms fixtures from the old ConstructDatabase
Django app and loads them into the new compounds.constructs schema.

Usage:
    python manage.py import_legacy_constructs \\
        --auth-fixture auth.json \\
        --constructs-fixture ConstructDatabase.json

    python manage.py import_legacy_constructs \\
        --auth-fixture auth.json \\
        --constructs-fixture ConstructDatabase.json \\
        --dry-run --verbose
"""

import json
import re
from datetime import datetime

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.utils import timezone


# Datetime fields that need timezone info
DATETIME_FIELDS = {
    'created_at', 'updated_at', 'createdAt', 'updatedAt',
}

# Model name mappings: old -> new
MODEL_MAPPINGS = {
    # ConstructDatabase -> compounds.constructs
    'ConstructDatabase.project': 'constructs.constructproject',
    'ConstructDatabase.plasmid': 'constructs.plasmid',
    'ConstructDatabase.protein': 'constructs.protein',
    'ConstructDatabase.proteinsynonym': 'constructs.proteinsynonym',
    'ConstructDatabase.proteinuse': 'constructs.proteinuse',
    'ConstructDatabase.cassette': 'constructs.cassette',
    'ConstructDatabase.cassetteuse': 'constructs.cassetteuse',
    'ConstructDatabase.sequencingresult': 'constructs.sequencingresult',
    'ConstructDatabase.expressiontagtype': 'constructs.expressiontagtype',
    'ConstructDatabase.protease': 'constructs.protease',
    'ConstructDatabase.expressiontaglocation': 'constructs.expressiontaglocation',
    'ConstructDatabase.expressiontag': 'constructs.expressiontag',
}


class Command(BaseCommand):
    """Import legacy ConstructDatabase fixtures."""

    help = "Import legacy ConstructDatabase fixtures"

    def add_arguments(self, parser):
        parser.add_argument(
            '--auth-fixture',
            type=str,
            help='Path to auth.json fixture (optional, for user import)',
        )
        parser.add_argument(
            '--constructs-fixture',
            type=str,
            required=True,
            help='Path to ConstructDatabase.json fixture',
        )
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Transform and validate without loading',
        )
        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Show detailed progress',
        )

    def handle(self, *args, **options):
        auth_fixture = options.get('auth_fixture')
        constructs_fixture = options['constructs_fixture']
        dry_run = options['dry_run']
        verbose = options['verbose']

        self.verbose = verbose
        self.stdout.write("\nImporting legacy constructs fixtures")
        self.stdout.write("-" * 60)

        if dry_run:
            self.stdout.write(self.style.WARNING("[DRY RUN - no changes will be saved]"))

        # Load and build user lookup
        user_lookup = {}
        if auth_fixture:
            self.log(f"Loading auth fixture: {auth_fixture}")
            auth_data = self._load_fixture(auth_fixture)
            user_lookup = self._build_user_lookup(auth_data)
            self.stdout.write(f"  Built user lookup with {len(user_lookup)} users")

        # Load and group fixture data by model
        self.stdout.write(f"\nProcessing constructs fixture: {constructs_fixture}")
        constructs_data = self._load_fixture(constructs_fixture)

        # Group records by model type
        records_by_model = {}
        for record in constructs_data:
            old_model = record['model']
            if old_model in MODEL_MAPPINGS:
                new_model = MODEL_MAPPINGS[old_model]
                if new_model not in records_by_model:
                    records_by_model[new_model] = []
                records_by_model[new_model].append(record)

        # Show counts
        for model, records in sorted(records_by_model.items()):
            self.stdout.write(f"    {model}: {len(records)}")

        if dry_run:
            self.stdout.write(self.style.WARNING("\n[DRY RUN] No changes made"))
            return

        # Import in dependency order, preserving dates
        try:
            with transaction.atomic():
                self._import_reference_data(records_by_model)
                self._import_construct_projects(records_by_model)
                self._import_proteins(records_by_model)
                self._import_plasmids(records_by_model)
                self._import_cassettes(records_by_model)
                self._import_cassette_uses(records_by_model)
                self._import_protein_uses(records_by_model)
                self._import_sequencing_results(records_by_model)
                self._import_expression_tags(records_by_model)

                self.stdout.write(self.style.SUCCESS("\nImport completed successfully!"))
        except Exception as e:
            self.stdout.write(self.style.ERROR(f"\nImport failed: {e}"))
            raise

    def log(self, message):
        """Log message if verbose mode is enabled."""
        if self.verbose:
            self.stdout.write(f"  {message}")

    def _load_fixture(self, path):
        """Load and parse a fixture file."""
        try:
            with open(path, 'r') as f:
                content = f.read()
                # Handle Django debug output in fixture
                start = content.find('[')
                if start > 0:
                    content = content[start:]
                return json.loads(content)
        except FileNotFoundError:
            raise CommandError(f"Fixture file not found: {path}")
        except json.JSONDecodeError as e:
            raise CommandError(f"Invalid JSON in {path}: {e}")

    def _build_user_lookup(self, auth_data):
        """Build a lookup from username to user info."""
        lookup = {}
        for record in auth_data:
            if record['model'] == 'auth.user':
                fields = record['fields']
                username = fields.get('username', '')
                lookup[username] = {
                    'email': fields.get('email', ''),
                    'first_name': fields.get('first_name', ''),
                    'last_name': fields.get('last_name', ''),
                }
        return lookup

    def _make_timezone_aware(self, value):
        """Add UTC timezone to naive datetime strings."""
        if not value or not isinstance(value, str):
            return value
        if value.endswith('Z') or '+' in value[-6:]:
            return value
        if re.match(r'^\d{4}-\d{2}-\d{2}', value):
            return value + 'Z'
        return value

    def _parse_datetime(self, value):
        """Parse datetime string from fixture, making it timezone-aware."""
        if not value:
            return None
        if isinstance(value, datetime):
            if value.tzinfo is None:
                return timezone.make_aware(value, timezone.utc)
            return value
        # Parse string datetime
        try:
            # Handle Z suffix
            value_str = value.replace('Z', '+00:00')
            dt = datetime.fromisoformat(value_str)
            if dt.tzinfo is None:
                dt = timezone.make_aware(dt, timezone.utc)
            return dt
        except (ValueError, AttributeError):
            return None

    def _get_field(self, fields, *keys):
        """Get field value trying multiple key names."""
        for key in keys:
            if key in fields:
                return fields[key]
        return None

    def _import_reference_data(self, records_by_model):
        """Import reference data (expression tag types, proteases, locations)."""
        from compounds.constructs.models import (
            ExpressionTagType, Protease, ExpressionTagLocation
        )

        # ExpressionTagType
        records = records_by_model.get('constructs.expressiontagtype', [])
        for record in records:
            fields = record['fields']
            ExpressionTagType.objects.update_or_create(
                id=record['pk'],
                defaults={'name': fields.get('name', '')}
            )
        self.stdout.write(f"  Imported {len(records)} ExpressionTagTypes")

        # Protease
        records = records_by_model.get('constructs.protease', [])
        for record in records:
            fields = record['fields']
            Protease.objects.update_or_create(
                id=record['pk'],
                defaults={'name': fields.get('name', '')}
            )
        self.stdout.write(f"  Imported {len(records)} Proteases")

        # ExpressionTagLocation
        records = records_by_model.get('constructs.expressiontaglocation', [])
        for record in records:
            fields = record['fields']
            ExpressionTagLocation.objects.update_or_create(
                id=record['pk'],
                defaults={'name': fields.get('name', '')}
            )
        self.stdout.write(f"  Imported {len(records)} ExpressionTagLocations")

    def _import_construct_projects(self, records_by_model):
        """Import ConstructProject records preserving dates."""
        from compounds.constructs.models import ConstructProject

        records = records_by_model.get('constructs.constructproject', [])
        for record in records:
            fields = record['fields']
            created_at = self._get_field(fields, 'createdAt', 'created_at')

            ConstructProject.objects.update_or_create(
                id=record['pk'],
                defaults={
                    'name': fields.get('name', ''),
                    'parent_id': self._get_field(fields, 'parent_id', 'parent'),
                }
            )
            # Preserve original created_at
            legacy_created_at = self._parse_datetime(created_at)
            if legacy_created_at:
                ConstructProject.objects.filter(id=record['pk']).update(
                    created_at=legacy_created_at
                )
        self.stdout.write(f"  Imported {len(records)} ConstructProjects")

    def _import_proteins(self, records_by_model):
        """Import Protein and ProteinSynonym records preserving dates."""
        from compounds.constructs.models import Protein, ProteinSynonym

        # Proteins
        records = records_by_model.get('constructs.protein', [])
        for record in records:
            fields = record['fields']
            created_at = self._get_field(fields, 'createdAt', 'created_at')

            Protein.objects.update_or_create(
                id=record['pk'],
                defaults={
                    'uniprot_id': self._get_field(fields, 'uniprotId', 'uniprot_id', 'name') or '',
                }
            )
            legacy_created_at = self._parse_datetime(created_at)
            if legacy_created_at:
                Protein.objects.filter(id=record['pk']).update(
                    created_at=legacy_created_at
                )
        self.stdout.write(f"  Imported {len(records)} Proteins")

        # ProteinSynonyms
        records = records_by_model.get('constructs.proteinsynonym', [])
        for record in records:
            fields = record['fields']
            created_at = self._get_field(fields, 'createdAt', 'created_at')

            ProteinSynonym.objects.update_or_create(
                id=record['pk'],
                defaults={
                    'name': fields.get('name', ''),
                    'protein_id': self._get_field(fields, 'protein_id', 'protein'),
                }
            )
            legacy_created_at = self._parse_datetime(created_at)
            if legacy_created_at:
                ProteinSynonym.objects.filter(id=record['pk']).update(
                    created_at=legacy_created_at
                )
        self.stdout.write(f"  Imported {len(records)} ProteinSynonyms")

    def _import_plasmids(self, records_by_model):
        """Import Plasmid records preserving dates."""
        from compounds.constructs.models import Plasmid

        records = records_by_model.get('constructs.plasmid', [])
        for record in records:
            fields = record['fields']
            created_at = self._get_field(fields, 'createdAt', 'created_at')

            Plasmid.objects.update_or_create(
                id=record['pk'],
                defaults={
                    'ncn_id': self._get_field(fields, 'ncnId', 'ncn_id'),
                    'name': fields.get('name', ''),
                    'parent_id': self._get_field(fields, 'parent_id', 'parent'),
                    'project_id': self._get_field(fields, 'project_id', 'project'),
                    'genbank_file': self._get_field(fields, 'snapgeneFile', 'genbank_file') or '',
                }
            )
            legacy_created_at = self._parse_datetime(created_at)
            if legacy_created_at:
                Plasmid.objects.filter(id=record['pk']).update(
                    created_at=legacy_created_at
                )
        self.stdout.write(f"  Imported {len(records)} Plasmids")

    def _import_cassettes(self, records_by_model):
        """Import Cassette records preserving dates."""
        from compounds.constructs.models import Cassette

        records = records_by_model.get('constructs.cassette', [])
        for record in records:
            fields = record['fields']
            created_at = self._get_field(fields, 'createdAt', 'created_at')

            Cassette.objects.update_or_create(
                id=record['pk'],
                defaults={
                    'protein_id': self._get_field(fields, 'protein_id', 'protein'),
                    'start': fields.get('start', 0),
                    'end': fields.get('end', 0),
                }
            )
            legacy_created_at = self._parse_datetime(created_at)
            if legacy_created_at:
                Cassette.objects.filter(id=record['pk']).update(
                    created_at=legacy_created_at
                )
        self.stdout.write(f"  Imported {len(records)} Cassettes")

    def _import_cassette_uses(self, records_by_model):
        """Import CassetteUse records preserving dates."""
        from compounds.constructs.models import CassetteUse

        records = records_by_model.get('constructs.cassetteuse', [])
        for record in records:
            fields = record['fields']
            created_at = self._get_field(fields, 'createdAt', 'created_at')

            CassetteUse.objects.update_or_create(
                id=record['pk'],
                defaults={
                    'cassette_id': self._get_field(fields, 'cassette_id', 'cassette'),
                    'plasmid_id': self._get_field(fields, 'plasmid_id', 'plasmid'),
                    'alignment_file': self._get_field(fields, 'alignmentFile', 'alignment_file') or '',
                }
            )
            legacy_created_at = self._parse_datetime(created_at)
            if legacy_created_at:
                CassetteUse.objects.filter(id=record['pk']).update(
                    created_at=legacy_created_at
                )
        self.stdout.write(f"  Imported {len(records)} CassetteUses")

    def _import_protein_uses(self, records_by_model):
        """Import ProteinUse records preserving dates."""
        from compounds.constructs.models import ProteinUse

        records = records_by_model.get('constructs.proteinuse', [])
        for record in records:
            fields = record['fields']
            created_at = self._get_field(fields, 'createdAt', 'created_at')

            ProteinUse.objects.update_or_create(
                id=record['pk'],
                defaults={
                    'protein_id': self._get_field(fields, 'protein_id', 'protein'),
                    'project_id': self._get_field(fields, 'project_id', 'project'),
                }
            )
            legacy_created_at = self._parse_datetime(created_at)
            if legacy_created_at:
                ProteinUse.objects.filter(id=record['pk']).update(
                    created_at=legacy_created_at
                )
        self.stdout.write(f"  Imported {len(records)} ProteinUses")

    def _import_sequencing_results(self, records_by_model):
        """Import SequencingResult records preserving dates."""
        from compounds.constructs.models import SequencingResult

        records = records_by_model.get('constructs.sequencingresult', [])
        for record in records:
            fields = record['fields']
            created_at = self._get_field(fields, 'createdAt', 'created_at')

            SequencingResult.objects.update_or_create(
                id=record['pk'],
                defaults={
                    'cassette_use_id': self._get_field(fields, 'cassetteUse_id', 'cassette_use'),
                    'plasmid_id': self._get_field(fields, 'plasmid_id', 'plasmid'),
                    'file': fields.get('file', ''),
                }
            )
            legacy_created_at = self._parse_datetime(created_at)
            if legacy_created_at:
                SequencingResult.objects.filter(id=record['pk']).update(
                    created_at=legacy_created_at
                )
        self.stdout.write(f"  Imported {len(records)} SequencingResults")

    def _import_expression_tags(self, records_by_model):
        """Import ExpressionTag records preserving dates."""
        from compounds.constructs.models import (
            ExpressionTag, ExpressionTagType, Protease, ExpressionTagLocation, CassetteUse
        )

        records = records_by_model.get('constructs.expressiontag', [])
        count = 0
        skipped = 0

        for record in records:
            fields = record['fields']
            created_at = self._get_field(fields, 'createdAt', 'created_at')

            tag_type_id = self._get_field(fields, 'expressionTagType_id', 'expression_tag_type')
            protease_id = self._get_field(fields, 'protease_id', 'protease')
            location_id = self._get_field(fields, 'location_id', 'location')
            cassette_use_id = self._get_field(fields, 'cassetteUse_id', 'cassette_use')

            # Validate references exist
            tag_type = ExpressionTagType.objects.filter(id=tag_type_id).first() if tag_type_id else None
            location = ExpressionTagLocation.objects.filter(id=location_id).first() if location_id else None
            cassette_use = CassetteUse.objects.filter(id=cassette_use_id).first() if cassette_use_id else None

            if not tag_type or not location or not cassette_use:
                skipped += 1
                continue

            protease = Protease.objects.filter(id=protease_id).first() if protease_id else None

            ExpressionTag.objects.update_or_create(
                id=record['pk'],
                defaults={
                    'expression_tag_type': tag_type,
                    'protease': protease,
                    'location': location,
                    'cassette_use': cassette_use,
                }
            )
            legacy_created_at = self._parse_datetime(created_at)
            if legacy_created_at:
                ExpressionTag.objects.filter(id=record['pk']).update(
                    created_at=legacy_created_at
                )
            count += 1

        self.stdout.write(f"  Imported {count} ExpressionTags, skipped {skipped}")
