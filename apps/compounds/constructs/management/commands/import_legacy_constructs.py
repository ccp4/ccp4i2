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
import tempfile
from pathlib import Path

from django.core.management import call_command
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction


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
        parser.add_argument(
            '--output-dir',
            type=str,
            help='Directory to save transformed fixtures (for inspection)',
        )

    def handle(self, *args, **options):
        auth_fixture = options.get('auth_fixture')
        constructs_fixture = options['constructs_fixture']
        dry_run = options['dry_run']
        verbose = options['verbose']
        output_dir = options.get('output_dir')

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

        # Use temp dir or specified output dir
        if output_dir:
            temp_dir = Path(output_dir)
            temp_dir.mkdir(parents=True, exist_ok=True)
        else:
            temp_dir = Path(tempfile.mkdtemp(prefix='constructs_import_'))

        self.stdout.write(f"  Working directory: {temp_dir}")

        try:
            transformed_files = []

            # Transform constructs fixture
            self.stdout.write(f"\nProcessing constructs fixture: {constructs_fixture}")
            constructs_data = self._load_fixture(constructs_fixture)
            constructs_records = self._transform_records(constructs_data, user_lookup)

            # Sort by model for FK dependencies
            # Order matters for foreign key relationships:
            # 1. Reference data (expression tag types, proteases, locations)
            # 2. Projects
            # 3. Proteins and synonyms
            # 4. Plasmids
            # 5. Cassettes
            # 6. CassetteUses and ProteinUses
            # 7. SequencingResults and ExpressionTags
            model_order = [
                'constructs.expressiontagtype',
                'constructs.protease',
                'constructs.expressiontaglocation',
                'constructs.constructproject',
                'constructs.protein',
                'constructs.proteinsynonym',
                'constructs.plasmid',
                'constructs.cassette',
                'constructs.cassetteuse',
                'constructs.proteinuse',
                'constructs.sequencingresult',
                'constructs.expressiontag',
            ]
            constructs_records.sort(
                key=lambda r: model_order.index(r['model']) if r['model'] in model_order else 999
            )

            constructs_file = temp_dir / 'compounds_constructs.json'
            with open(constructs_file, 'w') as f:
                json.dump(constructs_records, f, indent=2)
            transformed_files.append(constructs_file)
            self.stdout.write(f"  Transformed {len(constructs_records)} constructs records")

            # Show counts by model
            model_counts = {}
            for record in constructs_records:
                model = record['model']
                model_counts[model] = model_counts.get(model, 0) + 1
            for model, count in sorted(model_counts.items()):
                self.stdout.write(f"    {model}: {count}")

            # Load fixtures
            if not dry_run:
                self.stdout.write("\nLoading fixtures into database...")
                for fixture_file in transformed_files:
                    self.stdout.write(f"  Loading: {fixture_file.name}")
                    call_command('loaddata', str(fixture_file), verbosity=1 if verbose else 0)

                self.stdout.write(self.style.SUCCESS("\nImport completed successfully!"))
            else:
                self.stdout.write(self.style.WARNING(
                    f"\n[DRY RUN] Transformed fixtures saved to: {temp_dir}"
                ))

        finally:
            # Clean up temp files if not using output_dir
            if not output_dir and not dry_run:
                import shutil
                shutil.rmtree(temp_dir, ignore_errors=True)

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

    def _transform_records(self, data, user_lookup):
        """Transform records to new schema."""
        records = []
        for record in data:
            transformed = self._transform_record(record, user_lookup)
            if transformed:
                records.append(transformed)
        return records

    def _transform_record(self, record, user_lookup):
        """Transform a single fixture record."""
        old_model = record['model']

        if old_model not in MODEL_MAPPINGS:
            return None

        new_model = MODEL_MAPPINGS[old_model]
        pk = record['pk']
        old_fields = record['fields']
        new_fields = {}

        for old_key, value in old_fields.items():
            new_key = self._transform_field_name(old_key, new_model)

            if new_key is None:
                continue

            # Transform file paths
            if 'file' in old_key.lower() or old_key in ('snapgeneFile', 'alignmentFile'):
                value = self._transform_file_path(value, new_model)

            # Handle created_by natural key references
            if old_key == 'createdBy' and isinstance(value, list):
                # Natural key format - skip as we don't have the user
                value = None
                new_key = 'created_by'
            elif old_key == 'createdBy_id':
                # Direct FK - skip as legacy user IDs don't match
                value = None
                new_key = 'created_by'

            # Normalize field names from camelCase
            if new_key == 'cassette_use' and old_key == 'cassetteUse_id':
                new_key = 'cassette_use'
            if new_key == 'expression_tag_type' and old_key == 'expressionTagType_id':
                new_key = 'expression_tag_type'

            # Make datetime fields timezone-aware
            if new_key in DATETIME_FIELDS or old_key in DATETIME_FIELDS:
                value = self._make_timezone_aware(value)

            if new_key and value is not None:
                new_fields[new_key] = value

        return {
            'model': new_model,
            'pk': pk,
            'fields': new_fields
        }

    def _transform_field_name(self, old_name, model_type=None):
        """Transform field names to new schema.

        Args:
            old_name: The legacy field name
            model_type: The target model type (e.g., 'constructs.plasmid')
        """
        mappings = {
            # Common fields (camelCase -> snake_case)
            'createdAt': 'created_at',
            'updatedAt': 'updated_at',
            'createdBy': 'created_by',
            'createdBy_id': 'created_by',

            # Project fields
            'parent_id': 'parent',

            # Plasmid fields
            'ncnId': 'ncn_id',
            'snapgeneFile': 'genbank_file',
            'project_id': 'project',

            # Protein fields
            'uniprotId': 'uniprot_id',
            'protein_id': 'protein',

            # Cassette fields
            'start': 'start',
            'end': 'end',

            # CassetteUse fields
            'cassette_id': 'cassette',
            'plasmid_id': 'plasmid',
            'alignmentFile': 'alignment_file',

            # SequencingResult fields
            'cassetteUse_id': 'cassette_use',
            'file': 'file',

            # ExpressionTag fields
            'expressionTagType_id': 'expression_tag_type',
            'protease_id': 'protease',
            'location_id': 'location',

            # Drop these fields
            'original_id': None,
        }

        return mappings.get(old_name, old_name)

    def _transform_file_path(self, old_path, model_type):
        """Transform file paths to new structure.

        Legacy pattern: ConstructDatabase/{filename} or ConstructDatabase/snapgene/{filename}
        New pattern: ConstructDatabase/NCLCON-XXXXXXXX/{filename}

        Note: The actual file relocation should be handled by migrate-media.sh
        This only transforms the path references in the fixture.
        """
        if not old_path:
            return old_path

        # Most paths should already be in ConstructDatabase/
        # No transformation needed as the new upload paths use the same structure

        return old_path

    def _make_timezone_aware(self, value):
        """Add UTC timezone to naive datetime strings."""
        if not value or not isinstance(value, str):
            return value
        if value.endswith('Z') or '+' in value[-6:]:
            return value
        if re.match(r'^\d{4}-\d{2}-\d{2}', value):
            return value + 'Z'
        return value
