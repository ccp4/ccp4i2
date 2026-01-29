# Generated migration for MolecularProperties and MolecularPropertyThreshold models

from django.db import migrations, models
import uuid


def calculate_properties_for_existing_compounds(apps, schema_editor):
    """
    Calculate molecular properties for all existing compounds.

    This data migration runs RDKit calculations for compounds that
    don't yet have MolecularProperties records.
    """
    Compound = apps.get_model('registry', 'Compound')
    MolecularProperties = apps.get_model('registry', 'MolecularProperties')

    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import Descriptors, Lipinski
    except ImportError:
        print("RDKit not available - skipping property calculation")
        return

    compounds = Compound.objects.filter(
        molecular_properties__isnull=True
    ).exclude(smiles__isnull=True).exclude(smiles='')

    count = 0
    for compound in compounds.iterator():
        smiles = compound.smiles or compound.rdkit_smiles
        if not smiles:
            continue

        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            MolecularProperties.objects.create(
                compound=compound,
                heavy_atom_count=Lipinski.HeavyAtomCount(mol),
                hbd=Descriptors.NumHDonors(mol),
                hba=Descriptors.NumHAcceptors(mol),
                clogp=Descriptors.MolLogP(mol),
                tpsa=Descriptors.TPSA(mol),
                rotatable_bonds=Descriptors.NumRotatableBonds(mol),
                fraction_sp3=Descriptors.FractionCSP3(mol),
                rdkit_version=rdBase.rdkitVersion,
            )
            count += 1
        except Exception as e:
            print(f"Error calculating properties for {compound.id}: {e}")
            continue

    print(f"Calculated molecular properties for {count} compounds")


def seed_default_thresholds(apps, schema_editor):
    """
    Seed the database with Lipinski Rule of 5 default thresholds.
    """
    MolecularPropertyThreshold = apps.get_model('registry', 'MolecularPropertyThreshold')

    defaults = [
        {'property_name': 'molecular_weight', 'direction': 'above',
         'amber_threshold': 450, 'red_threshold': 500},
        {'property_name': 'clogp', 'direction': 'above',
         'amber_threshold': 4, 'red_threshold': 5},
        {'property_name': 'hbd', 'direction': 'above',
         'amber_threshold': 4, 'red_threshold': 5},
        {'property_name': 'hba', 'direction': 'above',
         'amber_threshold': 8, 'red_threshold': 10},
        {'property_name': 'tpsa', 'direction': 'above',
         'amber_threshold': 120, 'red_threshold': 140},
        {'property_name': 'rotatable_bonds', 'direction': 'above',
         'amber_threshold': 8, 'red_threshold': 10},
        {'property_name': 'fraction_sp3', 'direction': 'below',
         'amber_threshold': 0.3, 'red_threshold': 0.2},
    ]

    for threshold in defaults:
        MolecularPropertyThreshold.objects.get_or_create(
            property_name=threshold['property_name'],
            defaults=threshold
        )

    print(f"Seeded {len(defaults)} default property thresholds")


def reverse_thresholds(apps, schema_editor):
    """Remove seeded thresholds on reverse migration."""
    MolecularPropertyThreshold = apps.get_model('registry', 'MolecularPropertyThreshold')
    MolecularPropertyThreshold.objects.all().delete()


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0014_add_audit_fields'),
    ]

    operations = [
        # Create MolecularProperties table
        migrations.CreateModel(
            name='MolecularProperties',
            fields=[
                ('compound', models.OneToOneField(
                    on_delete=models.deletion.CASCADE,
                    primary_key=True,
                    related_name='molecular_properties',
                    serialize=False,
                    to='registry.compound'
                )),
                ('heavy_atom_count', models.IntegerField(
                    blank=True,
                    help_text='Number of non-hydrogen atoms (for ligand efficiency)',
                    null=True
                )),
                ('hbd', models.IntegerField(
                    blank=True,
                    help_text='Hydrogen bond donors (Lipinski: ≤5)',
                    null=True
                )),
                ('hba', models.IntegerField(
                    blank=True,
                    help_text='Hydrogen bond acceptors (Lipinski: ≤10)',
                    null=True
                )),
                ('clogp', models.FloatField(
                    blank=True,
                    help_text='Calculated LogP - Wildman-Crippen method (Lipinski: ≤5)',
                    null=True
                )),
                ('tpsa', models.FloatField(
                    blank=True,
                    help_text='Topological polar surface area in Å² (oral bioavailability: ≤140)',
                    null=True
                )),
                ('rotatable_bonds', models.IntegerField(
                    blank=True,
                    help_text='Number of rotatable bonds (oral bioavailability: ≤10)',
                    null=True
                )),
                ('fraction_sp3', models.FloatField(
                    blank=True,
                    help_text='Fraction of sp3 carbons (Fsp3) - 3D complexity indicator',
                    null=True
                )),
                ('calculated_at', models.DateTimeField(auto_now=True)),
                ('rdkit_version', models.CharField(
                    blank=True,
                    help_text='RDKit version used for calculation',
                    max_length=32,
                    null=True
                )),
            ],
            options={
                'verbose_name': 'Molecular Properties',
                'verbose_name_plural': 'Molecular Properties',
            },
        ),

        # Create MolecularPropertyThreshold table
        migrations.CreateModel(
            name='MolecularPropertyThreshold',
            fields=[
                ('id', models.UUIDField(
                    default=uuid.uuid4,
                    editable=False,
                    primary_key=True,
                    serialize=False
                )),
                ('property_name', models.CharField(
                    choices=[
                        ('molecular_weight', 'Molecular Weight'),
                        ('heavy_atom_count', 'Heavy Atom Count'),
                        ('hbd', 'H-Bond Donors'),
                        ('hba', 'H-Bond Acceptors'),
                        ('clogp', 'cLogP'),
                        ('tpsa', 'TPSA'),
                        ('rotatable_bonds', 'Rotatable Bonds'),
                        ('fraction_sp3', 'Fraction sp3'),
                    ],
                    help_text='Which molecular property this threshold applies to',
                    max_length=32,
                    unique=True
                )),
                ('direction', models.CharField(
                    choices=[
                        ('above', 'Above threshold is bad'),
                        ('below', 'Below threshold is bad'),
                    ],
                    default='above',
                    help_text='Whether values above or below threshold are concerning',
                    max_length=8
                )),
                ('amber_threshold', models.FloatField(
                    help_text='Value at which property turns amber (warning)'
                )),
                ('red_threshold', models.FloatField(
                    help_text='Value at which property turns red (danger)'
                )),
                ('enabled', models.BooleanField(
                    default=True,
                    help_text='Whether RAG coloring is active for this property'
                )),
            ],
            options={
                'verbose_name': 'Property Threshold',
                'verbose_name_plural': 'Property Thresholds',
                'ordering': ['property_name'],
            },
        ),

        # Data migrations
        migrations.RunPython(
            calculate_properties_for_existing_compounds,
            reverse_code=migrations.RunPython.noop,
        ),
        migrations.RunPython(
            seed_default_thresholds,
            reverse_code=reverse_thresholds,
        ),
    ]
