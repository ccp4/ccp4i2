"""
Data migration to fix genbank_file paths after blob storage migration.

The legacy storage had a double-nested structure:
  ConstructDatabase/ConstructDatabase/NCLCON-*/filename.gb

Django expects:
  ConstructDatabase/NCLCON-*/filename.gb

This migration updates all Plasmid.genbank_file paths to remove the extra nesting.
"""

from django.db import migrations


def fix_genbank_file_paths(apps, schema_editor):
    """Remove double-nested ConstructDatabase from file paths."""
    Plasmid = apps.get_model('constructs', 'Plasmid')

    # Find all plasmids with the double-nested path
    plasmids_to_fix = Plasmid.objects.filter(
        genbank_file__startswith='ConstructDatabase/ConstructDatabase/'
    )

    count = 0
    for plasmid in plasmids_to_fix:
        old_path = plasmid.genbank_file
        # Replace the double-nested path with single path
        new_path = old_path.replace(
            'ConstructDatabase/ConstructDatabase/',
            'ConstructDatabase/',
            1  # Only replace the first occurrence
        )
        plasmid.genbank_file = new_path
        plasmid.save(update_fields=['genbank_file'])
        count += 1

    if count:
        print(f"\n  Fixed {count} genbank_file paths")


def reverse_fix(apps, schema_editor):
    """Reverse the path fix (add back the nesting)."""
    Plasmid = apps.get_model('constructs', 'Plasmid')

    # Find all plasmids with the single path that should have been double-nested
    # Only fix paths that start with ConstructDatabase/NCLCON (the expected pattern)
    plasmids_to_reverse = Plasmid.objects.filter(
        genbank_file__startswith='ConstructDatabase/NCLCON'
    )

    count = 0
    for plasmid in plasmids_to_reverse:
        old_path = plasmid.genbank_file
        new_path = old_path.replace(
            'ConstructDatabase/',
            'ConstructDatabase/ConstructDatabase/',
            1
        )
        plasmid.genbank_file = new_path
        plasmid.save(update_fields=['genbank_file'])
        count += 1

    if count:
        print(f"\n  Reversed {count} genbank_file paths")


class Migration(migrations.Migration):

    dependencies = [
        ('constructs', '0002_seed_reference_data'),
    ]

    operations = [
        migrations.RunPython(fix_genbank_file_paths, reverse_fix),
    ]
