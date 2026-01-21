"""
Seed ADME protocol definitions for vendor data imports.

Creates Protocol records for NCU ADME assay types:
- NCU Liver Microsome Stability (LM)
- NCU Blood/Serum Stability (BS)
- NCU GSH Stability
- NCU Caco-2 Permeability
"""

from django.db import migrations


def seed_adme_protocols(apps, schema_editor):
    """Create ADME protocol records."""
    Protocol = apps.get_model('assays', 'Protocol')

    protocols = [
        {
            'name': 'NCU liver_microsome_stability',
            'analysis_method': 'pharmaron_adme',
            'comments': (
                'NCU Liver Microsome Stability (LM) assay.\n'
                'Measures metabolic stability using human and mouse liver microsomes.\n'
                'Key metrics: t1/2 (min), CLint, CLH, Extraction Ratio.\n'
                'File format: ADME-NCU-LM-YYYYMMDD.xlsx'
            ),
        },
        {
            'name': 'NCU blood_serum_stability',
            'analysis_method': 'pharmaron_adme',
            'comments': (
                'NCU Blood/Serum Stability (BS) assay.\n'
                'Measures compound stability in blood or serum over time.\n'
                'Key metrics: t1/2 (min), % remaining at time points.\n'
                'File format: ADME-NCU-BS-YYYYMMDD.xlsx'
            ),
        },
        {
            'name': 'NCU gsh_stability',
            'analysis_method': 'pharmaron_adme',
            'comments': (
                'NCU GSH (Glutathione) Stability assay.\n'
                'Screens for reactive metabolite formation via GSH trapping.\n'
                'Key metrics: t1/2 with/without GSH, GSH reactivity assessment.\n'
                'File format: ADME-NCU-GSH Stability-YYYYMMDD.xlsx'
            ),
        },
        {
            'name': 'NCU caco2_permeability',
            'analysis_method': 'pharmaron_adme',
            'comments': (
                'NCU Caco-2 Permeability assay.\n'
                'Measures intestinal permeability using Caco-2 cell monolayers.\n'
                'Key metrics: Papp A-B, Papp B-A, Efflux Ratio, Recovery.\n'
                'File format: ADME-NCU-Caco-2 Permeability-YYYYMMDD.xlsx'
            ),
        },
    ]

    for protocol_data in protocols:
        Protocol.objects.get_or_create(
            name=protocol_data['name'],
            defaults=protocol_data,
        )


def remove_adme_protocols(apps, schema_editor):
    """Remove ADME protocol records (reverse migration)."""
    Protocol = apps.get_model('assays', 'Protocol')

    Protocol.objects.filter(name__startswith='NCU ').filter(
        name__in=[
            'NCU liver_microsome_stability',
            'NCU blood_serum_stability',
            'NCU gsh_stability',
            'NCU caco2_permeability',
        ]
    ).delete()


class Migration(migrations.Migration):
    """Seed ADME protocol definitions."""

    dependencies = [
        ('assays', '0007_rename_svg_file_to_plot_image'),
    ]

    operations = [
        migrations.RunPython(
            seed_adme_protocols,
            reverse_code=remove_adme_protocols,
        ),
    ]
