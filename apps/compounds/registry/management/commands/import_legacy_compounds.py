"""
Management command to import legacy RegisterCompounds and AssayCompounds fixtures.

This command transforms fixtures from the old RegisterCompounds/AssayCompounds
Django apps and loads them into the new compounds.registry/compounds.assays schema.

Usage:
    python manage.py import_legacy_compounds \\
        --auth-fixture auth.json \\
        --registry-fixture RegisterCompounds.json \\
        --assays-fixture AssayCompounds.json

    python manage.py import_legacy_compounds \\
        --auth-fixture auth.json \\
        --registry-fixture RegisterCompounds.json \\
        --assays-fixture AssayCompounds.json \\
        --dry-run --verbose
"""

import json
import os
import re
import tempfile
from datetime import datetime
from pathlib import Path

from django.core.management import call_command
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.utils import timezone


# Datetime fields that need timezone info
DATETIME_FIELDS = {
    'created_at', 'updated_at', 'modified_at', 'registered_at', 'uploaded_at',
    'creation_date', 'reg_date', 'modified_date',
}

# Model name mappings: old -> new
MODEL_MAPPINGS = {
    # RegisterCompounds -> compounds.registry
    'RegisterCompounds.regprojects': 'registry.target',
    'RegisterCompounds.regsuppliers': 'registry.supplier',
    'RegisterCompounds.regdata': 'registry.compound',
    'RegisterCompounds.regbatch': 'registry.batch',
    'RegisterCompounds.regbatchqcfile': 'registry.batchqcfile',
    'RegisterCompounds.regdatatemplate': 'registry.compoundtemplate',

    # AssayCompounds -> compounds.assays
    'AssayCompounds.dilutionseries': 'assays.dilutionseries',
    'AssayCompounds.protocol': 'assays.protocol',
    'AssayCompounds.experiment': 'assays.assay',
    'AssayCompounds.dataseries': 'assays.dataseries',
    'AssayCompounds.analysis': 'assays.analysisresult',
    'AssayCompounds.hypothesis': 'assays.hypothesis',
}


class Command(BaseCommand):
    """Import legacy RegisterCompounds and AssayCompounds fixtures."""

    help = "Import legacy compounds fixtures (RegisterCompounds, AssayCompounds)"

    def add_arguments(self, parser):
        parser.add_argument(
            '--auth-fixture',
            type=str,
            help='Path to auth.json fixture (optional, for user import)',
        )
        parser.add_argument(
            '--registry-fixture',
            type=str,
            help='Path to RegisterCompounds.json fixture',
        )
        parser.add_argument(
            '--assays-fixture',
            type=str,
            help='Path to AssayCompounds.json fixture',
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
        registry_fixture = options.get('registry_fixture')
        assays_fixture = options.get('assays_fixture')
        dry_run = options['dry_run']
        verbose = options['verbose']
        output_dir = options.get('output_dir')

        if not registry_fixture and not assays_fixture:
            raise CommandError("Must provide at least one of --registry-fixture or --assays-fixture")

        self.verbose = verbose
        self.stdout.write("\nImporting legacy compounds fixtures")
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
            temp_dir = Path(tempfile.mkdtemp(prefix='compounds_import_'))

        self.stdout.write(f"  Working directory: {temp_dir}")

        try:
            transformed_files = []

            # Transform registry fixture
            if registry_fixture:
                self.stdout.write(f"\nProcessing registry fixture: {registry_fixture}")
                registry_data = self._load_fixture(registry_fixture)
                registry_records = self._transform_records(registry_data, user_lookup, 'registry')

                # Sort by model for FK dependencies
                model_order = ['registry.supplier', 'registry.target', 'registry.compound',
                               'registry.batch', 'registry.batchqcfile', 'registry.compoundtemplate']
                registry_records.sort(
                    key=lambda r: model_order.index(r['model']) if r['model'] in model_order else 999
                )

                registry_file = temp_dir / 'compounds_registry.json'
                with open(registry_file, 'w') as f:
                    json.dump(registry_records, f, indent=2)
                transformed_files.append(registry_file)
                self.stdout.write(f"  Transformed {len(registry_records)} registry records")

            # Transform assays fixture
            if assays_fixture:
                self.stdout.write(f"\nProcessing assays fixture: {assays_fixture}")
                assays_data = self._load_fixture(assays_fixture)
                assays_records = self._transform_records(assays_data, user_lookup, 'assays')

                # Sort by model for FK dependencies
                model_order = ['assays.dilutionseries', 'assays.protocol', 'assays.assay',
                               'assays.analysisresult', 'assays.dataseries', 'assays.hypothesis']
                assays_records.sort(
                    key=lambda r: model_order.index(r['model']) if r['model'] in model_order else 999
                )

                assays_file = temp_dir / 'compounds_assays.json'
                with open(assays_file, 'w') as f:
                    json.dump(assays_records, f, indent=2)
                transformed_files.append(assays_file)
                self.stdout.write(f"  Transformed {len(assays_records)} assay records")

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

    def _transform_records(self, data, user_lookup, app_type):
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
            if 'file' in old_key.lower() or old_key in ('svgFile', 'dataFile'):
                value = self._transform_file_path(value, new_model)

            # Handle created_by natural key references
            if old_key == 'created_by' and isinstance(value, list):
                if value and len(value) > 0:
                    if 'legacy_registered_by' not in new_fields and new_model == 'registry.compound':
                        new_fields['legacy_registered_by'] = value[0]
                value = None

            # Normalize stereo_comment
            if new_key == 'stereo_comment':
                value = (value.lower().replace(' ', '_').replace('-', '_')) if value else 'unset'

            # Normalize analysis_method
            if new_key == 'analysis_method' and value:
                value = value.lower().replace('-', '_').replace(' ', '_')

            # Normalize status
            if new_key == 'status' and value:
                value = value.lower()

            # Unwrap legacy {"data": ...} wrapper
            if new_key in ('extracted_data', 'skip_points') and isinstance(value, dict):
                if 'data' in value:
                    value = value['data']
                if new_key == 'extracted_data' and isinstance(value, dict):
                    try:
                        max_idx = max(int(k) for k in value.keys()) if value else -1
                        value = [value.get(str(i)) for i in range(max_idx + 1)]
                    except (ValueError, TypeError):
                        pass

            # Make datetime fields timezone-aware
            if new_key in DATETIME_FIELDS or old_key in DATETIME_FIELDS:
                value = self._make_timezone_aware(value)

            if new_key:
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
            model_type: The target model type (e.g., 'assays.dataseries')
        """
        mappings = {
            'project_name': 'name',
            'creation_date': 'created_at',
            'reg_id': 'reg_number',
            'project_id': 'target',
            'reg_date': 'registered_at',
            'labbook_id': 'labbook_number',
            'page_id': 'page_number',
            'compound_id': 'compound_number',
            'user_name': 'legacy_registered_by',
            'mw': 'molecular_weight',
            'inchi1': 'inchi',
            'inchi1_qualifier': None,
            'stereo_comments': 'stereo_comment',
            'pklFile': None,  # Drop pickle files - no longer used, mol objects generated on-demand from SMILES
            'svgFile': 'svg_file',
            'rdkitSmiles': 'rdkit_smiles',
            'aliases': None,
            'original_id': None,
            'modified_date': 'modified_at',
            'regdata_id': 'compound',
            'batch_number': 'batch_number',
            'mwt': 'molecular_weight',
            'batch_id': 'batch',
            'mol2d': 'mol2d',
            'shortName': None,
            'created_by': 'created_by',
            'created_at': 'created_at',
            'analysisMethod': 'analysis_method',
            'preferredDilutions': 'preferred_dilutions',
            'pherastarTableChoice': None,  # Field removed
            'dataFile': 'data_file',
            'experiment': 'assay',
            'compoundName': 'compound_name',
            'startColumn': 'start_column',
            'endColumn': 'end_column',
            'dilutionSeries': 'dilution_series',
            'extractedData': 'extracted_data',
            'skipPoints': 'skip_points',
            'modelMinyURL': 'model_url',
            'updated_at': 'updated_at',
            'productRegId': 'product_compound',
            'completionNotes': 'completion_notes',
            'project': 'target',
        }

        # Model-specific overrides
        # DataSeries uses plot_image for image files (not svg_file)
        if model_type == 'assays.dataseries' and old_name == 'svgFile':
            return 'plot_image'

        return mappings.get(old_name, old_name)

    def _transform_file_path(self, old_path, model_type):
        """Transform file paths to new structure.

        Batch QC files are relocated from:
            RegBatchQCFile_NCL-XXXXX/{uuid}_{filename}
        to:
            RegisterCompounds/BatchQCFiles/NCL-XXXXX/{uuid}_{filename}
        """
        if not old_path:
            return old_path

        # Handle batch QC file path transformation
        # Old: RegBatchQCFile_NCL-00029551/uuid_filename.pdf
        # New: RegisterCompounds/BatchQCFiles/NCL-00029551/uuid_filename.pdf
        if old_path.startswith('RegBatchQCFile_'):
            # Extract compound ID and rest of path
            # e.g., RegBatchQCFile_NCL-00029551/uuid_file.pdf
            parts = old_path.split('/', 1)
            if len(parts) == 2:
                dir_name = parts[0]  # RegBatchQCFile_NCL-00029551
                file_part = parts[1]  # uuid_file.pdf
                compound_id = dir_name.replace('RegBatchQCFile_', '')  # NCL-00029551
                return f'RegisterCompounds/BatchQCFiles/{compound_id}/{file_part}'
            else:
                # Just the directory name without a file
                compound_id = old_path.replace('RegBatchQCFile_', '')
                return f'RegisterCompounds/BatchQCFiles/{compound_id}'

        # Other path mappings (unchanged)
        path_mappings = {
            'RegisterCompounds/svg/': 'compounds/registry/svg/',
            'AssayCompounds/Experiments/': 'compounds/assays/data/',
            'AssayCompounds/svg/': 'compounds/assays/svg/',
        }

        for old_prefix, new_prefix in path_mappings.items():
            if old_prefix in old_path:
                return old_path.replace(old_prefix, new_prefix)

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
