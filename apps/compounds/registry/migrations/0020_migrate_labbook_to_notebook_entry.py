"""
Data migration to convert existing compound labbook references to LabNotebookEntry.

For each Compound with (supplier_id, labbook_number, page_number) set:
1. Find-or-create a LabNotebookEntry with sequence_number derived from labbook_number
2. Set compound.notebook_entry to the created entry

The sequence_number heuristic: use labbook_number * 1000 + page_number to ensure
uniqueness when multiple pages exist in the same labbook.
"""

from django.db import migrations


def migrate_labbook_to_notebook_entry(apps, schema_editor):
    """
    Convert existing compound labbook references to LabNotebookEntry records.
    """
    Compound = apps.get_model('registry', 'Compound')
    LabNotebookEntry = apps.get_model('registry', 'LabNotebookEntry')

    # Find compounds with labbook info but no notebook_entry
    compounds_with_labbook = Compound.objects.filter(
        supplier__isnull=False,
        labbook_number__isnull=False,
        notebook_entry__isnull=True
    ).select_related('supplier')

    created_entries = 0
    updated_compounds = 0
    warnings = []

    for compound in compounds_with_labbook.iterator():
        # Calculate sequence_number: labbook_number * 1000 + page_number
        # This ensures uniqueness for paper notebook records
        page = compound.page_number or 0
        sequence_number = compound.labbook_number * 1000 + page

        # Find or create the LabNotebookEntry
        entry, created = LabNotebookEntry.objects.get_or_create(
            supplier=compound.supplier,
            sequence_number=sequence_number,
            defaults={
                'labbook_number': compound.labbook_number,
                'page_number': compound.page_number,
            }
        )

        if created:
            created_entries += 1

        # Link compound to the entry
        compound.notebook_entry = entry
        compound.save(update_fields=['notebook_entry'])
        updated_compounds += 1

    # Report compounds with labbook info but no supplier (can't create entry)
    orphan_count = Compound.objects.filter(
        supplier__isnull=True,
        labbook_number__isnull=False,
        notebook_entry__isnull=True
    ).count()

    if orphan_count > 0:
        warnings.append(
            f"WARNING: {orphan_count} compounds have labbook_number set but no supplier - "
            f"cannot create LabNotebookEntry for these. Manual review required."
        )

    print(f"Created {created_entries} LabNotebookEntry records")
    print(f"Updated {updated_compounds} compounds with notebook_entry references")
    for warning in warnings:
        print(warning)


def reverse_migration(apps, schema_editor):
    """
    Reverse the migration by clearing notebook_entry references.

    Note: This does NOT delete the LabNotebookEntry records to preserve
    any data that may have been manually added.
    """
    Compound = apps.get_model('registry', 'Compound')
    updated = Compound.objects.filter(
        notebook_entry__isnull=False
    ).update(notebook_entry=None)
    print(f"Cleared notebook_entry from {updated} compounds")


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0019_eln_notebook_entry'),
    ]

    operations = [
        migrations.RunPython(
            migrate_labbook_to_notebook_entry,
            reverse_code=reverse_migration,
        ),
    ]
