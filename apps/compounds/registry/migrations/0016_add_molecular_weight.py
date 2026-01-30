# Generated migration for molecular_weight field

from django.db import migrations, models


def backfill_molecular_weight(apps, schema_editor):
    """
    Calculate and populate molecular_weight for all existing MolecularProperties.
    """
    MolecularProperties = apps.get_model('registry', 'MolecularProperties')

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        # RDKit not available, skip backfill
        return

    for props in MolecularProperties.objects.select_related('compound').all():
        compound = props.compound
        smiles = compound.smiles or compound.rdkit_smiles
        if not smiles:
            continue

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                props.molecular_weight = Descriptors.ExactMolWt(mol)
                props.save(update_fields=['molecular_weight'])
        except Exception:
            continue


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0015_molecularproperties_thresholds'),
    ]

    operations = [
        migrations.AddField(
            model_name='molecularproperties',
            name='molecular_weight',
            field=models.FloatField(
                blank=True,
                help_text='Exact molecular weight in Daltons (Lipinski: ≤500)',
                null=True
            ),
        ),
        migrations.RunPython(backfill_molecular_weight, migrations.RunPython.noop),
    ]
