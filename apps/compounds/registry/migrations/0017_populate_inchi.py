"""
Data migration to populate InChI field for existing compounds.

This migration calculates InChI for all compounds that have SMILES but
no InChI value set. Uses RDKit for calculation.
"""

from django.db import migrations


def populate_inchi(apps, schema_editor):
    """Calculate and populate InChI for compounds with SMILES but no InChI."""
    Compound = apps.get_model('registry', 'Compound')

    try:
        from rdkit import Chem
        from rdkit.Chem import inchi as rdkit_inchi
    except ImportError:
        print("RDKit not available - skipping InChI population")
        return

    # Get compounds with SMILES but no InChI
    compounds = Compound.objects.filter(
        smiles__isnull=False,
    ).exclude(smiles='').filter(inchi__isnull=True)

    total = compounds.count()
    print(f"Populating InChI for {total} compounds...")

    updated = 0
    errors = 0

    for compound in compounds.iterator():
        try:
            mol = Chem.MolFromSmiles(compound.smiles)
            if mol:
                compound.inchi = rdkit_inchi.MolToInchi(mol)
                compound.save(update_fields=['inchi'])
                updated += 1
            else:
                errors += 1
        except Exception as e:
            errors += 1

        # Progress indicator
        if (updated + errors) % 1000 == 0:
            print(f"  Processed {updated + errors}/{total}...")

    print(f"InChI population complete: {updated} updated, {errors} errors")


def reverse_inchi(apps, schema_editor):
    """Reverse migration - clear InChI values."""
    Compound = apps.get_model('registry', 'Compound')
    Compound.objects.all().update(inchi=None)


class Migration(migrations.Migration):
    dependencies = [
        ('registry', '0016_add_molecular_weight'),
    ]

    operations = [
        migrations.RunPython(populate_inchi, reverse_inchi),
    ]
