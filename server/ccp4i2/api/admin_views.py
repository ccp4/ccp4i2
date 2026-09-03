"""
Admin API views for CCP4i2 app.

Provides endpoints for administrative tasks like importing legacy fixtures.
These endpoints require platform admin permissions.
"""

import json
import logging
import os
import tempfile
from io import StringIO
from pathlib import Path

from django.core.management import call_command
from django.db import transaction
from django.conf import settings
from rest_framework import status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAdminUser
from rest_framework.response import Response

logger = logging.getLogger(__name__)


def _remap_dirs(request):
    """Extract (from, to) remap tuple from the request, or None."""
    remap_from = request.data.get('remap_from', '').strip()
    remap_to = request.data.get('remap_to', '').strip()
    return (remap_from, remap_to) if remap_from and remap_to else None


def _resolve_sqlite_source(request):
    """Resolve the legacy sqlite source from an upload or a server-side path.

    Returns ``(path, temp_path)`` where ``temp_path`` is the caller-owned temp
    file to unlink afterwards (None for a server-side path), or ``None`` if
    neither a file nor a path was supplied.
    """
    db_file = request.FILES.get('sqlite_db')
    db_path = request.data.get('db_path', '').strip()
    if not db_file and not db_path:
        return None
    if db_file:
        with tempfile.NamedTemporaryFile(suffix='.sqlite', delete=False) as tmp:
            for chunk in db_file.chunks():
                tmp.write(chunk)
            return tmp.name, tmp.name
    return Path(db_path).expanduser(), None


def _inplace_migration_disabled_response():
    """403 for adopting a legacy installation's projects WHERE THEY LIE.

    The distinction that matters is the mode, not the endpoint. Copying a legacy
    project to a new root reads the old installation and writes somewhere else,
    so it cannot damage what the user is still working with; adopting in place
    makes the old directories live under the new app while the old app may still
    be using them. Agreed with Paul Bond and Stuart McNicholas, 2026-09-02: the
    second is withdrawn for the alpha, the first is not.

    Callers read the flag with ``getattr(..., False)`` deliberately: a settings
    module that has never heard of it (test_settings, a bespoke deployment) must
    get the SAFE answer, not an AttributeError and a 500.
    """
    return Response(
        {
            "success": False,
            "error": "in_place_migration_disabled",
            "message": (
                "Adopting a legacy CCP4i2 installation's projects where they sit "
                "is disabled in this alpha, because it operates on projects your "
                "existing CCP4i2 may still be using. Copy them instead: pass "
                "copy_files with a dest_root, or use Projects -> Import, which "
                "accepts a project directory or a project zip and never writes "
                "to the original."
            ),
        },
        status=status.HTTP_403_FORBIDDEN,
    )


def _wrong_database_kind_response(db_path):
    """403/400 when the file is not a legacy CCP4i2 database.

    Without this, pointing the legacy importer at a database written by THIS
    application fails one query at a time on missing tables, and the user sees a
    stack of OperationalErrors rather than "you have picked the wrong file, and
    here is the right route for it".
    """
    from ccp4i2.db.import_sqlite import describe_sqlite_database

    described = describe_sqlite_database(db_path)
    if described["kind"] == "legacy":
        return None
    if not Path(db_path).expanduser().is_file():
        # "You gave me a path that isn't there" is a different complaint, and
        # the endpoints already answer it with 404. Don't shadow it with a
        # 400 about the file's contents.
        return None
    return Response(
        {
            "success": False,
            "error": f"not_a_legacy_database:{described['kind']}",
            "message": described["message"]
            or "That file is not a legacy CCP4i2 database.",
            "detected": described,
        },
        status=status.HTTP_400_BAD_REQUEST,
    )


def _in_place_migration_refused(request):
    """True when this request asks to adopt projects in place, and may not."""
    if getattr(settings, "CCP4I2_ALLOW_INPLACE_MIGRATION", False):
        return False
    return str(request.data.get("copy_files", "false")).lower() != "true"


@api_view(['POST'])
@permission_classes([IsAdminUser])
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
    # Loading a dumpdata fixture has no copy mode: it adopts whatever directories
    # the fixture names, so there is no safe variant to let through.
    if not getattr(settings, "CCP4I2_ALLOW_INPLACE_MIGRATION", False):
        return _inplace_migration_disabled_response()

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


@api_view(['POST'])
@permission_classes([IsAdminUser])
def import_sqlite(request):
    """
    Import legacy CCP4i2 data from a SQLite database file.

    Accepts multipart form data with:
    - sqlite_db: SQLite database file upload
    - dry_run: If true, validate without loading (optional)
    - remap_from: Directory path to remap from (optional)
    - remap_to: Directory path to remap to (optional)

    Or JSON body with:
    - db_path: Server-side path to SQLite database (e.g. ~/.CCP4I2/db/database.sqlite)
    - dry_run, remap_from, remap_to as above
    """
    if _in_place_migration_refused(request):
        return _inplace_migration_disabled_response()

    from ccp4i2.db.import_sqlite import SQLiteImporter, StructuralIssuesError

    src = _resolve_sqlite_source(request)
    if src is None:
        return Response(
            {'error': 'bad_input',
             'reason': 'Must provide sqlite_db file upload or db_path'},
            status=status.HTTP_400_BAD_REQUEST,
        )
    actual_path, temp_path = src

    wrong_kind = _wrong_database_kind_response(actual_path)
    if wrong_kind is not None:
        if temp_path:
            os.unlink(temp_path)
        return wrong_kind

    dry_run = str(request.data.get('dry_run', 'false')).lower() == 'true'
    copy_files = str(request.data.get('copy_files', 'false')).lower() == 'true'
    allow_structural = str(
        request.data.get('allow_structural_issues', 'false')
    ).lower() == 'true'
    dest_root = request.data.get('dest_root', '').strip() or None

    try:
        importer = SQLiteImporter(
            db_path=actual_path,
            remap_dirs=_remap_dirs(request),
            dry_run=dry_run,
            continue_on_error=True,
            copy_files=copy_files,
            dest_root=Path(dest_root).expanduser() if dest_root else None,
            allow_structural_issues=allow_structural,
        )
        result = importer.run()

        if not dry_run and not result['errors']:
            logger.info(
                f"SQLite import by {getattr(request.user, 'email', 'local')}: "
                f"{result['stats']}"
            )

        if result['errors']:
            return Response(
                {'error': 'import_failed', 'reason': 'Import completed with errors',
                 **result},
                status=status.HTTP_400_BAD_REQUEST,
            )
        return Response(result, status=status.HTTP_200_OK)

    except StructuralIssuesError as e:
        # Blocking structural issues not acknowledged — the UI maps this back to
        # the validation step. 400 + machine-readable discriminator (house style).
        return Response(
            {'error': 'structural_issues_unacknowledged',
             'reason': 'Blocking structural issues must be acknowledged before import',
             'structure': e.structure, 'plan': e.plan},
            status=status.HTTP_400_BAD_REQUEST,
        )
    except FileNotFoundError as e:
        return Response({'error': 'db_not_found', 'reason': str(e)},
                        status=status.HTTP_404_NOT_FOUND)
    except Exception as e:
        logger.exception("Error in import_sqlite")
        return Response(
            {'error': 'internal_error', 'reason': str(e)},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    finally:
        if temp_path:
            Path(temp_path).unlink(missing_ok=True)


@api_view(['POST'])
@permission_classes([IsAdminUser])
def validate_sqlite(request):
    """
    Validate a legacy CCP4i2 SQLite database against the filesystem.

    Checks that project directories, job directories, files, and imported
    file sources all exist on disk. Also checks referential integrity and
    data quality. Read-only — does not modify any database.

    Accepts multipart form data with:
    - sqlite_db: SQLite database file upload

    Or JSON body with:
    - db_path: Server-side path to SQLite database

    Optional:
    - remap_from / remap_to: Remap directory prefixes before checking
    - copy_files: intended migration mode — needed to compute the per-project
      plan and decide which structural issues are in scope (default false)
    - dest_root: override for the copy destination root
    """
    from ccp4i2.db.import_sqlite import SQLiteValidator

    src = _resolve_sqlite_source(request)
    if src is None:
        return Response(
            {'error': 'bad_input',
             'reason': 'Must provide sqlite_db file upload or db_path'},
            status=status.HTTP_400_BAD_REQUEST,
        )
    actual_path, temp_path = src

    wrong_kind = _wrong_database_kind_response(actual_path)
    if wrong_kind is not None:
        if temp_path:
            os.unlink(temp_path)
        return wrong_kind

    copy_files = str(request.data.get('copy_files', 'false')).lower() == 'true'
    dest_root = request.data.get('dest_root', '').strip() or None

    try:
        validator = SQLiteValidator(
            db_path=actual_path,
            remap_dirs=_remap_dirs(request),
            copy_files=copy_files,
            dest_root=Path(dest_root).expanduser() if dest_root else None,
        )
        report = validator.run()
        return Response(report)

    except FileNotFoundError as e:
        return Response({'error': 'db_not_found', 'reason': str(e)},
                        status=status.HTTP_404_NOT_FOUND)
    except Exception as e:
        logger.exception("Error in validate_sqlite")
        return Response(
            {'error': 'internal_error', 'reason': str(e)},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    finally:
        if temp_path:
            Path(temp_path).unlink(missing_ok=True)


@api_view(['GET'])
@permission_classes([IsAdminUser])
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
