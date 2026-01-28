"""
Add 'import_type' field to Protocol and migrate data from 'analysis_method'.

This migration:
1. Adds the new 'import_type' field with choices: raw_data, ms_intact, table_of_values, pharmaron_adme
2. Populates import_type from existing analysis_method values
3. For legacy hill_langmuir_fix_* variants, also populates fitting_parameters with appropriate constraints
4. Updates the 4PL fitting script to handle boolean fix_top/fix_bottom
"""

from django.db import migrations, models
from pathlib import Path


def load_script(filename):
    """Load script content from fitting_scripts directory."""
    script_dir = Path(__file__).parent.parent / 'fitting_scripts'
    script_path = script_dir / filename
    if script_path.exists():
        return script_path.read_text()
    return None


def migrate_analysis_method_to_import_type(apps, schema_editor):
    """Populate import_type from analysis_method and migrate fitting constraints."""
    Protocol = apps.get_model('assays', 'Protocol')

    for protocol in Protocol.objects.all():
        analysis_method = protocol.analysis_method

        # Map analysis_method to import_type
        if analysis_method.startswith('hill_langmuir'):
            protocol.import_type = 'raw_data'

            # Migrate fitting constraints from legacy variants
            params = protocol.fitting_parameters or {}

            if 'fix_hill' in analysis_method:
                params['fix_hill'] = 1.0

            if 'fix_minmax' in analysis_method or 'fix_hill_minmax' in analysis_method:
                params['fix_top'] = True
                params['fix_bottom'] = True

            if params != (protocol.fitting_parameters or {}):
                protocol.fitting_parameters = params

        elif analysis_method == 'ms_intact':
            protocol.import_type = 'ms_intact'
        elif analysis_method == 'table_of_values':
            protocol.import_type = 'table_of_values'
        elif analysis_method == 'pharmaron_adme':
            protocol.import_type = 'pharmaron_adme'
        else:
            # Default fallback
            protocol.import_type = 'raw_data'

        protocol.save()


def reverse_import_type_migration(apps, schema_editor):
    """Reverse migration: clear import_type field (it will be removed)."""
    # The field will be removed by the reverse AddField operation.
    # No explicit data cleanup needed since the field itself is removed.
    pass


def update_fitting_scripts(apps, schema_editor):
    """Update 4PL script to handle boolean fix_top/fix_bottom."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')

    four_pl_script = load_script('four_parameter_logistic.py')
    if four_pl_script:
        FittingMethod.objects.filter(slug='four-parameter-logistic').update(
            script=four_pl_script
        )


def noop(apps, schema_editor):
    """No-op reverse for script update."""
    pass


class Migration(migrations.Migration):
    """Add import_type field and migrate data from analysis_method."""

    dependencies = [
        ('assays', '0017_update_fitting_scripts'),
    ]

    operations = [
        # Add the new import_type field
        migrations.AddField(
            model_name='protocol',
            name='import_type',
            field=models.CharField(
                choices=[
                    ('raw_data', 'Raw Data (Dose-Response)'),
                    ('ms_intact', 'MS-Intact'),
                    ('table_of_values', 'Table of Values'),
                    ('pharmaron_adme', 'Pharmaron ADME'),
                ],
                default='raw_data',
                help_text='Type of data: raw dose-response for fitting, or pre-analyzed imports',
                max_length=50,
            ),
        ),
        # Migrate existing data
        migrations.RunPython(
            migrate_analysis_method_to_import_type,
            reverse_code=reverse_import_type_migration,
        ),
        # Update analysis_method help text to mark it deprecated
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
                help_text='DEPRECATED: Use import_type instead. Kept for backward compatibility.',
                max_length=50,
            ),
        ),
        # Update 4PL fitting script to handle boolean fix_top/fix_bottom
        migrations.RunPython(update_fitting_scripts, noop),
    ]
