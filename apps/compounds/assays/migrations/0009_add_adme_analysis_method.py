"""
Add 'pharmaron_adme' analysis method choice and update existing Pharmaron/NCU ADME protocols.
"""

from django.db import migrations, models


def update_ncu_protocols_to_pharmaron_adme(apps, schema_editor):
    """Update existing NCU ADME protocols to use 'pharmaron_adme' analysis method."""
    Protocol = apps.get_model('assays', 'Protocol')

    # Update NCU/Pharmaron ADME protocols that were previously set to 'table_of_values'
    ncu_protocol_names = [
        'NCU liver_microsome_stability',
        'NCU blood_serum_stability',
        'NCU gsh_stability',
        'NCU caco2_permeability',
    ]

    Protocol.objects.filter(
        name__in=ncu_protocol_names,
        analysis_method='table_of_values',
    ).update(analysis_method='pharmaron_adme')


def revert_ncu_protocols_to_tov(apps, schema_editor):
    """Revert NCU ADME protocols back to 'table_of_values' analysis method."""
    Protocol = apps.get_model('assays', 'Protocol')

    ncu_protocol_names = [
        'NCU liver_microsome_stability',
        'NCU blood_serum_stability',
        'NCU gsh_stability',
        'NCU caco2_permeability',
    ]

    Protocol.objects.filter(
        name__in=ncu_protocol_names,
        analysis_method='pharmaron_adme',
    ).update(analysis_method='table_of_values')


class Migration(migrations.Migration):
    """Add 'pharmaron_adme' analysis method and update NCU protocols."""

    dependencies = [
        ('assays', '0008_seed_adme_protocols'),
    ]

    operations = [
        # Add 'pharmaron_adme' to the analysis_method choices
        migrations.AlterField(
            model_name='protocol',
            name='analysis_method',
            field=models.CharField(
                choices=[
                    ('hill_langmuir', 'Hill-Langmuir'),
                    ('hill_langmuir_fix_hill', 'Hill-Langmuir (fixed Hill coefficient)'),
                    ('hill_langmuir_fix_hill_minmax', 'Hill-Langmuir (fixed Hill and min/max)'),
                    ('hill_langmuir_fix_minmax', 'Hill-Langmuir (fixed min/max)'),
                    ('ms_intact', 'MS-Intact'),
                    ('table_of_values', 'Table of values'),
                    ('pharmaron_adme', 'Pharmaron ADME'),
                ],
                default='hill_langmuir',
                help_text='Legacy analysis method (use fitting_method for new protocols)',
                max_length=50,
            ),
        ),
        # Update existing NCU protocols to use 'pharmaron_adme'
        migrations.RunPython(
            update_ncu_protocols_to_pharmaron_adme,
            reverse_code=revert_ncu_protocols_to_tov,
        ),
    ]
