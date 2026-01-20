# Data migration to populate original_filename for existing BatchQCFile records
# Uses reverse heuristics to extract original filename from stored file paths

import re
from pathlib import Path

from django.db import migrations


def extract_original_filename(file_path: str) -> str | None:
    """Extract original filename from file path using reverse heuristics.

    BatchQCFile pattern: RegisterCompounds/BatchQCFiles/{compound_id}/{uuid}_{filename}
    The UUID is 36 chars (8-4-4-4-12 with hyphens), followed by underscore and filename.
    """
    if not file_path:
        return None

    # Get the filename part (after last /)
    filename = Path(file_path).name

    # Pattern: {uuid}_{original_filename}
    # UUID is 36 chars, then underscore, then original filename
    if '_' in filename and len(filename) > 37:
        potential_uuid = filename[:36]
        if re.match(
            r'^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$',
            potential_uuid,
            re.IGNORECASE
        ):
            return filename[37:]  # Skip UUID and underscore

    # Fallback: just use the basename if no UUID pattern found
    return filename


def populate_original_filenames(apps, schema_editor):
    """Populate original_filename for existing BatchQCFile records."""
    BatchQCFile = apps.get_model('registry', 'BatchQCFile')

    updated_count = 0
    for qc_file in BatchQCFile.objects.filter(original_filename__isnull=True):
        if qc_file.file:
            original_filename = extract_original_filename(qc_file.file.name)
            if original_filename:
                qc_file.original_filename = original_filename
                qc_file.save(update_fields=['original_filename'])
                updated_count += 1

    if updated_count:
        print(f"\n  Populated original_filename for {updated_count} BatchQCFile records")


def reverse_populate(apps, schema_editor):
    """Reverse migration - clear original_filename values."""
    BatchQCFile = apps.get_model('registry', 'BatchQCFile')
    BatchQCFile.objects.update(original_filename=None)


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0008_batchqcfile_original_filename'),
    ]

    operations = [
        migrations.RunPython(
            populate_original_filenames,
            reverse_code=reverse_populate,
        ),
    ]
