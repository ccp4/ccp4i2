"""
Admin API views for CCP4i2 app.

Provides endpoints for administrative tasks like importing legacy fixtures.
These endpoints require platform admin permissions.
"""

import json
import logging
import tempfile
from io import StringIO
from pathlib import Path

from django.core.management import call_command
from django.db import transaction
from rest_framework import status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.response import Response

from users.permissions import IsPlatformAdmin

logger = logging.getLogger(__name__)


@api_view(['POST'])
@permission_classes([IsPlatformAdmin])
def import_legacy_ccp4i2(request):
    """
    Import legacy CCP4i2 dumpdata fixtures.

    Accepts multipart form data with:
    - ccp4i2_fixture: Legacy CCP4i2 dumpdata JSON file
    - dry_run: If true, validate without loading (optional)
    - remap_from: Directory path to remap from (optional)
    - remap_to: Directory path to remap to (optional)

    Returns summary of imported records or validation results.
    """
    fixture_file = request.FILES.get('ccp4i2_fixture')
    dry_run = request.data.get('dry_run', 'false').lower() == 'true'
    remap_from = request.data.get('remap_from', '').strip()
    remap_to = request.data.get('remap_to', '').strip()

    if not fixture_file:
        return Response(
            {'error': 'Must provide ccp4i2_fixture file'},
            status=status.HTTP_400_BAD_REQUEST
        )

    results = {
        'dry_run': dry_run,
        'stats': None,
        'errors': [],
    }

    try:
        # Save uploaded file temporarily
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as temp_file:
            content = fixture_file.read().decode('utf-8')
            temp_file.write(content)
            temp_path = temp_file.name

        # Capture output from management command
        out = StringIO()
        err = StringIO()

        # Build command arguments
        cmd_args = [temp_path]
        if dry_run:
            cmd_args.append('--dry-run')
        cmd_args.append('--continue-on-error')

        if remap_from and remap_to:
            cmd_args.extend(['--remap-dirs', remap_from, remap_to])

        try:
            call_command('import_legacy_ccp4i2', *cmd_args, stdout=out, stderr=err)

            # Parse output for stats
            output = out.getvalue()
            results['output'] = output

            # Extract stats from output (rough parsing)
            stats = {}
            for line in output.split('\n'):
                if ':' in line and 'imported' in line.lower():
                    parts = line.split(':')
                    if len(parts) >= 2:
                        key = parts[0].strip().strip('-').strip().lower().replace(' ', '_')
                        # Extract number from value
                        value_part = parts[1].split()[0] if parts[1].split() else '0'
                        try:
                            stats[key] = int(value_part)
                        except ValueError:
                            pass

            results['stats'] = stats

            # Check stderr for errors
            error_output = err.getvalue()
            if error_output:
                results['errors'] = [line for line in error_output.split('\n') if line.strip()]

            if not dry_run and not results['errors']:
                results['loaded'] = True
                logger.info(
                    f"Legacy CCP4i2 fixture imported by {getattr(request.user, 'email', 'local')}"
                )

        except Exception as e:
            results['errors'].append(f'Error running import command: {str(e)}')
            results['loaded'] = False

        # Clean up temp file
        Path(temp_path).unlink(missing_ok=True)

    except Exception as e:
        logger.exception("Error in import_legacy_ccp4i2")
        results['errors'].append(f'Unexpected error: {str(e)}')

    if results['errors']:
        return Response(results, status=status.HTTP_400_BAD_REQUEST)

    return Response(results)


@api_view(['GET'])
@permission_classes([IsPlatformAdmin])
def ccp4i2_import_status(request):
    """
    Get current database counts for CCP4i2 models.

    Useful for checking import results or current state.
    """
    from ccp4i2.db.models import (
        Project, ProjectGroup, ProjectGroupMembership, ProjectTag,
        Job, File, FileType, FileUse, FileImport, FileExport,
        JobValueKey, JobFloatValue, JobCharValue, XData,
    )

    return Response({
        'projects': {
            'total': Project.objects.count(),
            'groups': ProjectGroup.objects.count(),
            'memberships': ProjectGroupMembership.objects.count(),
            'tags': ProjectTag.objects.count(),
        },
        'jobs': {
            'total': Job.objects.count(),
            'value_keys': JobValueKey.objects.count(),
            'float_values': JobFloatValue.objects.count(),
            'char_values': JobCharValue.objects.count(),
            'xdata': XData.objects.count(),
        },
        'files': {
            'total': File.objects.count(),
            'types': FileType.objects.count(),
            'uses': FileUse.objects.count(),
            'imports': FileImport.objects.count(),
            'exports': FileExport.objects.count(),
        },
    })
