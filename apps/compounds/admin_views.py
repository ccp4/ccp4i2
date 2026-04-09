"""
Admin API views for compounds app.

Provides endpoints for administrative tasks like importing legacy fixtures.
These endpoints require platform admin permissions.
"""

import json
import logging
import os
import tempfile
from pathlib import Path

from django.core.management import call_command
from django.db import transaction
from rest_framework import status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.response import Response

from users.permissions import IsPlatformAdmin

# Note: IsPlatformAdmin already handles the no-auth case (returns True when
# CCP4I2_REQUIRE_AUTH is not set), so we don't need IsAuthenticated here.

logger = logging.getLogger(__name__)


# Import the transformation logic from the management command
from compounds.registry.management.commands.import_legacy_compounds import (
    MODEL_MAPPINGS,
    DATETIME_FIELDS,
)


def _make_timezone_aware(value):
    """Convert datetime string to timezone-aware format."""
    from django.utils import timezone
    if not value:
        return value
    if isinstance(value, str):
        if value.endswith('Z'):
            return value
        if '+' in value or value.endswith('Z'):
            return value
        return value + 'Z'
    return value


def _transform_field_name(old_name, model_type=None):
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
        'pklFile': None,  # Drop pickle files - no longer used
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


def _normalize_absolute_path(path):
    """
    Strip absolute path prefixes to get relative media path.

    Handles paths like:
    - /var/app/media/RegisterCompounds/svg/file.svg
    - /home/user/project/media/RegBatchQCFile_NCL-123/file.pdf
    - /data/media/AssayCompounds/Experiments/file.xlsx

    Returns the relative path from MEDIA_ROOT.
    """
    if not path or not path.startswith('/'):
        return path

    # Known directory markers that indicate start of relative media path
    media_markers = [
        'RegisterCompounds/',
        'RegBatchQCFile_',
        'AssayCompounds/',
        'compounds/',  # Already transformed paths
    ]

    for marker in media_markers:
        if marker in path:
            # Find where the marker starts and return from there
            idx = path.find(marker)
            relative_path = path[idx:]
            logger.debug(f"Normalized absolute path: {path} -> {relative_path}")
            return relative_path

    # Check for /media/ directory pattern
    if '/media/' in path:
        idx = path.find('/media/') + len('/media/')
        relative_path = path[idx:]
        logger.debug(f"Normalized path via /media/: {path} -> {relative_path}")
        return relative_path

    # Path is absolute but doesn't match known patterns - log warning
    logger.warning(f"Unable to normalize absolute path: {path}")
    return path


def _transform_file_path(old_path, model_type):
    """Transform file paths to new structure."""
    if not old_path:
        return old_path

    # First normalize any absolute paths to relative
    old_path = _normalize_absolute_path(old_path)

    path_mappings = {
        'RegisterCompounds/svg/': 'compounds/registry/svg/',
        'RegBatchQCFile_': 'compounds/registry/qc/',
        'AssayCompounds/Experiments/': 'compounds/assays/data/',
        'AssayCompounds/svg/': 'compounds/assays/svg/',
    }

    for old_prefix, new_prefix in path_mappings.items():
        if old_prefix in old_path:
            return old_path.replace(old_prefix, new_prefix)

    return old_path


def _transform_record(record, user_lookup, app_filter=None):
    """Transform a single fixture record."""
    old_model = record['model']

    if old_model not in MODEL_MAPPINGS:
        return None

    new_model = MODEL_MAPPINGS[old_model]

    # Filter by app if specified
    if app_filter and not new_model.startswith(app_filter):
        return None

    pk = record['pk']
    old_fields = record['fields']
    new_fields = {}

    for old_key, value in old_fields.items():
        new_key = _transform_field_name(old_key, new_model)

        if new_key is None:
            continue

        # Transform file paths
        if 'file' in old_key.lower() or old_key in ('svgFile', 'dataFile'):
            value = _transform_file_path(value, new_model)

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
            value = _make_timezone_aware(value)

        if new_key:
            new_fields[new_key] = value

    return {
        'model': new_model,
        'pk': pk,
        'fields': new_fields
    }


def _transform_records(records, user_lookup, app_filter=None):
    """Transform a list of fixture records."""
    transformed = []
    for record in records:
        result = _transform_record(record, user_lookup, app_filter)
        if result:
            transformed.append(result)
    return transformed


def _import_users_from_fixture(fixture_content, dry_run=False):
    """
    Import users from legacy auth.json fixture content.

    Returns dict with 'created', 'updated', 'skipped', 'errors' counts.
    """
    from django.contrib.auth import get_user_model
    from django.utils import timezone
    from users.models import UserProfile

    User = get_user_model()

    try:
        # Handle Django debug output in fixture
        start = fixture_content.find('[')
        if start > 0:
            fixture_content = fixture_content[start:]
        data = json.loads(fixture_content)
    except json.JSONDecodeError as e:
        return {'error': f'Invalid JSON: {str(e)}'}

    # Filter to auth.user records
    user_records = [r for r in data if r.get('model') == 'auth.user']

    stats = {'created': 0, 'updated': 0, 'skipped': 0, 'errors': [], 'total': len(user_records)}

    for record in user_records:
        try:
            fields = record.get('fields', {})
            pk = record.get('pk')
            username = fields.get('username', '')
            email = fields.get('email', '')

            if not username:
                stats['skipped'] += 1
                continue

            # Normalize email
            if not email and '@' in username:
                email = username

            # Check for existing user by PK first (to preserve FK references)
            existing = None
            try:
                existing = User.objects.get(pk=pk)
            except User.DoesNotExist:
                # Try by username or email
                existing = User.objects.filter(username=username).first()
                if not existing and email:
                    existing = User.objects.filter(email__iexact=email).first()

            # Parse date_joined
            date_joined = fields.get('date_joined')
            if date_joined:
                if not date_joined.endswith('Z') and '+' not in date_joined[-6:]:
                    date_joined = date_joined + 'Z'
                if date_joined.endswith('Z'):
                    date_joined = date_joined[:-1] + '+00:00'

            user_data = {
                'email': email,
                'first_name': fields.get('first_name', '')[:30],
                'last_name': fields.get('last_name', '')[:150],
                'is_active': fields.get('is_active', True),
                'is_staff': fields.get('is_staff', False),
                'is_superuser': fields.get('is_superuser', False),
            }

            if not dry_run:
                if existing:
                    for key, value in user_data.items():
                        setattr(existing, key, value)
                    existing.save()
                    user = existing
                    stats['updated'] += 1
                else:
                    # Create with specific PK to preserve FK references
                    user = User.objects.create(
                        pk=pk,
                        username=username,
                        **user_data
                    )
                    stats['created'] += 1

                # Update profile
                profile, _ = UserProfile.objects.get_or_create(user=user)
                profile.legacy_username = username
                profile.legacy_display_name = f"{fields.get('first_name', '')} {fields.get('last_name', '')}".strip()
                profile.imported_at = timezone.now()
                profile.save()
            else:
                # Dry run - just count
                if existing:
                    stats['updated'] += 1
                else:
                    stats['created'] += 1

        except Exception as e:
            stats['errors'].append(f"Error importing user {username}: {str(e)}")

    return stats


@api_view(['POST'])
@permission_classes([IsPlatformAdmin])
def import_legacy_fixtures(request):
    """
    Import legacy RegisterCompounds and/or AssayCompounds fixtures.

    Accepts multipart form data with JSON fixture files:
    - users_fixture: auth.json (optional, should be imported first)
    - registry_fixture: RegisterCompounds.json (optional)
    - assays_fixture: AssayCompounds.json (optional)
    - dry_run: If true, transform and validate without loading (optional)

    Returns summary of imported records or validation results.
    """
    users_file = request.FILES.get('users_fixture')
    registry_file = request.FILES.get('registry_fixture')
    assays_file = request.FILES.get('assays_fixture')
    dry_run = request.data.get('dry_run', 'false').lower() == 'true'

    if not users_file and not registry_file and not assays_file:
        return Response(
            {'error': 'Must provide at least one fixture file'},
            status=status.HTTP_400_BAD_REQUEST
        )

    results = {
        'dry_run': dry_run,
        'users': None,
        'registry': None,
        'assays': None,
        'errors': [],
    }

    user_lookup = {}

    try:
        # Process users fixture FIRST (to establish FK targets)
        if users_file:
            try:
                users_content = users_file.read().decode('utf-8')
                users_result = _import_users_from_fixture(users_content, dry_run)

                if 'error' in users_result:
                    results['errors'].append(f'Users fixture: {users_result["error"]}')
                else:
                    results['users'] = {
                        'total_records': users_result['total'],
                        'created': users_result['created'],
                        'updated': users_result['updated'],
                        'skipped': users_result['skipped'],
                    }
                    if users_result.get('errors'):
                        results['errors'].extend(users_result['errors'])
            except Exception as e:
                results['errors'].append(f'Error processing users_fixture: {str(e)}')

        with tempfile.TemporaryDirectory(prefix='compounds_import_') as temp_dir:
            temp_path = Path(temp_dir)
            transformed_files = []

            # Process registry fixture
            if registry_file:
                try:
                    registry_content = registry_file.read().decode('utf-8')
                    # Handle Django debug output in fixture (same as users)
                    start = registry_content.find('[')
                    if start > 0:
                        registry_content = registry_content[start:]
                    registry_data = json.loads(registry_content)

                    registry_records = _transform_records(registry_data, user_lookup, 'registry')

                    # Sort by model for FK dependencies
                    model_order = ['registry.supplier', 'registry.target', 'registry.compound',
                                   'registry.batch', 'registry.batchqcfile', 'registry.compoundtemplate']
                    registry_records.sort(
                        key=lambda r: model_order.index(r['model']) if r['model'] in model_order else 999
                    )

                    registry_output = temp_path / 'compounds_registry.json'
                    with open(registry_output, 'w') as f:
                        json.dump(registry_records, f, indent=2)
                    transformed_files.append(str(registry_output))

                    # Count by model type
                    model_counts = {}
                    for rec in registry_records:
                        model = rec['model']
                        model_counts[model] = model_counts.get(model, 0) + 1

                    results['registry'] = {
                        'total_records': len(registry_records),
                        'by_model': model_counts,
                    }
                except json.JSONDecodeError as e:
                    results['errors'].append(f'Invalid JSON in registry_fixture: {str(e)}')
                except Exception as e:
                    results['errors'].append(f'Error processing registry_fixture: {str(e)}')

            # Process assays fixture
            if assays_file:
                try:
                    assays_content = assays_file.read().decode('utf-8')
                    # Handle Django debug output in fixture (same as users)
                    start = assays_content.find('[')
                    if start > 0:
                        assays_content = assays_content[start:]
                    assays_data = json.loads(assays_content)

                    assays_records = _transform_records(assays_data, user_lookup, 'assays')

                    # Sort by model for FK dependencies
                    model_order = ['assays.dilutionseries', 'assays.protocol', 'assays.assay',
                                   'assays.analysisresult', 'assays.dataseries', 'assays.hypothesis']
                    assays_records.sort(
                        key=lambda r: model_order.index(r['model']) if r['model'] in model_order else 999
                    )

                    assays_output = temp_path / 'compounds_assays.json'
                    with open(assays_output, 'w') as f:
                        json.dump(assays_records, f, indent=2)
                    transformed_files.append(str(assays_output))

                    # Count by model type
                    model_counts = {}
                    for rec in assays_records:
                        model = rec['model']
                        model_counts[model] = model_counts.get(model, 0) + 1

                    results['assays'] = {
                        'total_records': len(assays_records),
                        'by_model': model_counts,
                    }
                except json.JSONDecodeError as e:
                    results['errors'].append(f'Invalid JSON in assays_fixture: {str(e)}')
                except Exception as e:
                    results['errors'].append(f'Error processing assays_fixture: {str(e)}')

            # Load fixtures if not dry run and no errors
            if not dry_run and not results['errors'] and transformed_files:
                try:
                    with transaction.atomic():
                        for fixture_file in transformed_files:
                            call_command('loaddata', fixture_file, verbosity=0)
                        results['loaded'] = True
                        logger.info(
                            f"Legacy fixtures imported by {getattr(request.user, 'email', 'local')}: "
                            f"users={results.get('users', {}).get('total_records', 0) if results.get('users') else 0}, "
                            f"registry={results.get('registry', {}).get('total_records', 0) if results.get('registry') else 0}, "
                            f"assays={results.get('assays', {}).get('total_records', 0) if results.get('assays') else 0}"
                        )
                except Exception as e:
                    results['errors'].append(f'Error loading fixtures: {str(e)}')
                    results['loaded'] = False

    except Exception as e:
        logger.exception("Error in import_legacy_fixtures")
        results['errors'].append(f'Unexpected error: {str(e)}')

    if results['errors']:
        return Response(results, status=status.HTTP_400_BAD_REQUEST)

    return Response(results)


@api_view(['GET'])
@permission_classes([IsPlatformAdmin])
def import_status(request):
    """
    Get current database counts for users and compounds models.

    Useful for checking import results or current state.
    """
    from django.contrib.auth import get_user_model
    from users.models import UserProfile
    from compounds.registry.models import Supplier, Target, Compound, Batch, BatchQCFile, CompoundTemplate
    from compounds.assays.models import (
        DilutionSeries, Protocol, Assay, DataSeries, AnalysisResult, Hypothesis
    )
    from compounds.constructs.models import (
        ConstructProject, Plasmid, Protein, ProteinSynonym, ProteinUse,
        Cassette, CassetteUse, SequencingResult, ExpressionTagType,
        Protease, ExpressionTagLocation, ExpressionTag
    )

    User = get_user_model()

    return Response({
        'users': {
            'total': User.objects.count(),
            'with_legacy_username': UserProfile.objects.exclude(legacy_username='').count(),
        },
        'registry': {
            'suppliers': Supplier.objects.count(),
            'targets': Target.objects.count(),
            'compounds': Compound.objects.count(),
            'batches': Batch.objects.count(),
            'batch_qc_files': BatchQCFile.objects.count(),
            'compound_templates': CompoundTemplate.objects.count(),
        },
        'assays': {
            'dilution_series': DilutionSeries.objects.count(),
            'protocols': Protocol.objects.count(),
            'assays': Assay.objects.count(),
            'data_series': DataSeries.objects.count(),
            'analysis_results': AnalysisResult.objects.count(),
            'hypotheses': Hypothesis.objects.count(),
        },
        'constructs': {
            'projects': ConstructProject.objects.count(),
            'plasmids': Plasmid.objects.count(),
            'proteins': Protein.objects.count(),
            'protein_synonyms': ProteinSynonym.objects.count(),
            'protein_uses': ProteinUse.objects.count(),
            'cassettes': Cassette.objects.count(),
            'cassette_uses': CassetteUse.objects.count(),
            'sequencing_results': SequencingResult.objects.count(),
            'expression_tag_types': ExpressionTagType.objects.count(),
            'proteases': Protease.objects.count(),
            'expression_tag_locations': ExpressionTagLocation.objects.count(),
            'expression_tags': ExpressionTag.objects.count(),
        },
    })


@api_view(['POST'])
@permission_classes([IsPlatformAdmin])
def import_constructs_fixtures(request):
    """
    Import legacy ConstructDatabase fixtures.

    Accepts multipart form data with JSON fixture files:
    - users_fixture: auth.json (optional, for user references)
    - constructs_fixture: ConstructDatabase.json (required)
    - dry_run: If true, transform and validate without loading (optional)

    Returns summary of imported records or validation results.
    """
    from compounds.constructs.management.commands.import_legacy_constructs import (
        MODEL_MAPPINGS as CONSTRUCTS_MODEL_MAPPINGS,
        DATETIME_FIELDS as CONSTRUCTS_DATETIME_FIELDS,
    )

    users_file = request.FILES.get('users_fixture')
    constructs_file = request.FILES.get('constructs_fixture')
    dry_run = request.data.get('dry_run', 'false').lower() == 'true'

    if not constructs_file:
        return Response(
            {'error': 'Must provide constructs_fixture file'},
            status=status.HTTP_400_BAD_REQUEST
        )

    results = {
        'dry_run': dry_run,
        'users': None,
        'constructs': None,
        'errors': [],
    }

    user_lookup = {}

    try:
        # Process users fixture if provided (to build lookup)
        if users_file:
            try:
                users_content = users_file.read().decode('utf-8')
                users_result = _import_users_from_fixture(users_content, dry_run)

                if 'error' in users_result:
                    results['errors'].append(f'Users fixture: {users_result["error"]}')
                else:
                    results['users'] = {
                        'total_records': users_result['total'],
                        'created': users_result['created'],
                        'updated': users_result['updated'],
                        'skipped': users_result['skipped'],
                    }
                    if users_result.get('errors'):
                        results['errors'].extend(users_result['errors'])
            except Exception as e:
                results['errors'].append(f'Error processing users_fixture: {str(e)}')

        with tempfile.TemporaryDirectory(prefix='constructs_import_') as temp_dir:
            temp_path = Path(temp_dir)

            # Process constructs fixture
            try:
                constructs_content = constructs_file.read().decode('utf-8')
                # Handle Django debug output in fixture
                start = constructs_content.find('[')
                if start > 0:
                    constructs_content = constructs_content[start:]
                constructs_data = json.loads(constructs_content)

                # Transform records
                constructs_records = []
                for record in constructs_data:
                    old_model = record.get('model', '')
                    if old_model not in CONSTRUCTS_MODEL_MAPPINGS:
                        continue

                    new_model = CONSTRUCTS_MODEL_MAPPINGS[old_model]
                    pk = record['pk']
                    old_fields = record.get('fields', {})
                    new_fields = {}

                    # Field name mappings
                    field_mappings = {
                        'createdAt': 'created_at',
                        'updatedAt': 'updated_at',
                        'createdBy': None,  # Skip - FK refs won't match
                        'createdBy_id': None,
                        'parent_id': 'parent',
                        'ncnId': 'ncn_id',
                        'snapgeneFile': 'genbank_file',
                        'project_id': 'project',
                        'uniprotId': 'uniprot_id',
                        'protein_id': 'protein',
                        'cassette_id': 'cassette',
                        'plasmid_id': 'plasmid',
                        'alignmentFile': 'alignment_file',
                        'cassetteUse_id': 'cassette_use',
                        'expressionTagType_id': 'expression_tag_type',
                        'protease_id': 'protease',
                        'location_id': 'location',
                        'original_id': None,  # Drop
                    }

                    for old_key, value in old_fields.items():
                        if old_key in field_mappings:
                            new_key = field_mappings[old_key]
                        else:
                            new_key = old_key

                        if new_key is None:
                            continue

                        # Make datetime fields timezone-aware
                        if new_key in CONSTRUCTS_DATETIME_FIELDS or old_key in CONSTRUCTS_DATETIME_FIELDS:
                            if value and isinstance(value, str):
                                if not value.endswith('Z') and '+' not in value[-6:]:
                                    value = value + 'Z'

                        if new_key and value is not None:
                            new_fields[new_key] = value

                    constructs_records.append({
                        'model': new_model,
                        'pk': pk,
                        'fields': new_fields
                    })

                # Sort by model for FK dependencies
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

                constructs_output = temp_path / 'compounds_constructs.json'
                with open(constructs_output, 'w') as f:
                    json.dump(constructs_records, f, indent=2)

                # Count by model type
                model_counts = {}
                for rec in constructs_records:
                    model = rec['model']
                    model_counts[model] = model_counts.get(model, 0) + 1

                results['constructs'] = {
                    'total_records': len(constructs_records),
                    'by_model': model_counts,
                }

                # Load fixtures if not dry run and no errors
                if not dry_run and not results['errors']:
                    try:
                        with transaction.atomic():
                            call_command('loaddata', str(constructs_output), verbosity=0)
                            results['loaded'] = True
                            logger.info(
                                f"Constructs fixtures imported by {getattr(request.user, 'email', 'local')}: "
                                f"total={len(constructs_records)}"
                            )
                    except Exception as e:
                        results['errors'].append(f'Error loading fixtures: {str(e)}')
                        results['loaded'] = False

            except json.JSONDecodeError as e:
                results['errors'].append(f'Invalid JSON in constructs_fixture: {str(e)}')
            except Exception as e:
                results['errors'].append(f'Error processing constructs_fixture: {str(e)}')

    except Exception as e:
        logger.exception("Error in import_constructs_fixtures")
        results['errors'].append(f'Unexpected error: {str(e)}')

    if results['errors']:
        return Response(results, status=status.HTTP_400_BAD_REQUEST)

    return Response(results)


# =============================================================================
# DANGER ZONE - Data Reset Endpoint
# =============================================================================
#
# This endpoint exists ONLY for the initial data migration phase.
# It should be DISABLED in production by removing the environment variable.
#
# To enable: Set CCP4I2_ALLOW_DB_RESET=true in the environment
# To disable: Remove or unset CCP4I2_ALLOW_DB_RESET
#
# =============================================================================

RESET_CONFIRMATION_PHRASE = "DELETE ALL COMPOUNDS DATA"


@api_view(['POST'])
@permission_classes([IsPlatformAdmin])
def reset_compounds_data(request):
    """
    ⚠️  DANGER: Delete all compounds registry and assays data.

    This endpoint is ONLY for use during initial data migration.
    It is gated by the CCP4I2_ALLOW_DB_RESET environment variable.

    Required POST parameters:
    - confirmation: Must be exactly "DELETE ALL COMPOUNDS DATA"

    Optional parameters:
    - include_users: If "true", also clears legacy user data (default: false)

    Returns counts of deleted records.
    """
    # Gate 1: Environment variable must be set
    if not os.environ.get('CCP4I2_ALLOW_DB_RESET', '').lower() == 'true':
        logger.warning(
            f"Reset attempt blocked - CCP4I2_ALLOW_DB_RESET not enabled. "
            f"User: {getattr(request.user, 'email', 'anonymous')}"
        )
        return Response(
            {
                'error': 'Database reset is disabled',
                'hint': 'Set CCP4I2_ALLOW_DB_RESET=true to enable (migration phase only)',
            },
            status=status.HTTP_403_FORBIDDEN
        )

    # Gate 2: Confirmation phrase must match exactly
    confirmation = request.data.get('confirmation', '')
    if confirmation != RESET_CONFIRMATION_PHRASE:
        return Response(
            {
                'error': 'Confirmation phrase does not match',
                'required': RESET_CONFIRMATION_PHRASE,
                'received': confirmation,
            },
            status=status.HTTP_400_BAD_REQUEST
        )

    include_users = request.data.get('include_users', 'false').lower() == 'true'

    # Import models
    from compounds.registry.models import (
        Supplier, Target, Compound, Batch, BatchQCFile, CompoundTemplate
    )
    from compounds.assays.models import (
        DilutionSeries, Protocol, Assay, DataSeries, AnalysisResult, Hypothesis
    )

    results = {
        'deleted': {
            'registry': {},
            'assays': {},
        },
        'users': None,
        'warning': '⚠️  This action is irreversible!',
    }

    user_email = getattr(request.user, 'email', 'anonymous')

    try:
        with transaction.atomic():
            # Delete in reverse dependency order

            # Assays (DataSeries -> AnalysisResult -> Assay -> Protocol -> DilutionSeries)
            # Hypothesis references Assay, so delete first
            results['deleted']['assays']['hypotheses'] = Hypothesis.objects.count()
            Hypothesis.objects.all().delete()

            results['deleted']['assays']['data_series'] = DataSeries.objects.count()
            DataSeries.objects.all().delete()

            results['deleted']['assays']['analysis_results'] = AnalysisResult.objects.count()
            AnalysisResult.objects.all().delete()

            results['deleted']['assays']['assays'] = Assay.objects.count()
            Assay.objects.all().delete()

            results['deleted']['assays']['protocols'] = Protocol.objects.count()
            Protocol.objects.all().delete()

            results['deleted']['assays']['dilution_series'] = DilutionSeries.objects.count()
            DilutionSeries.objects.all().delete()

            # Registry (BatchQCFile -> Batch -> Compound -> Target/Supplier)
            results['deleted']['registry']['batch_qc_files'] = BatchQCFile.objects.count()
            BatchQCFile.objects.all().delete()

            results['deleted']['registry']['batches'] = Batch.objects.count()
            Batch.objects.all().delete()

            results['deleted']['registry']['compound_templates'] = CompoundTemplate.objects.count()
            CompoundTemplate.objects.all().delete()

            results['deleted']['registry']['compounds'] = Compound.objects.count()
            Compound.objects.all().delete()

            results['deleted']['registry']['targets'] = Target.objects.count()
            Target.objects.all().delete()

            results['deleted']['registry']['suppliers'] = Supplier.objects.count()
            Supplier.objects.all().delete()

            # Optionally clear legacy user data
            if include_users:
                from users.models import UserProfile
                results['users'] = {
                    'profiles_cleared': UserProfile.objects.exclude(legacy_username='').count()
                }
                UserProfile.objects.all().update(
                    legacy_username='',
                    legacy_display_name='',
                    imported_at=None
                )

            logger.warning(
                f"🚨 COMPOUNDS DATA RESET by {user_email}: "
                f"registry={sum(results['deleted']['registry'].values())}, "
                f"assays={sum(results['deleted']['assays'].values())}, "
                f"include_users={include_users}"
            )

    except Exception as e:
        logger.exception(f"Error during compounds data reset by {user_email}")
        return Response(
            {'error': f'Reset failed: {str(e)}'},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR
        )

    return Response(results)
