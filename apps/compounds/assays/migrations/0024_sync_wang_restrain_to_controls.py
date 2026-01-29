# Data migration to sync Wang tight-binding script with restrain_to_controls feature
# Updates both script and input/output schemas

from django.db import migrations
from pathlib import Path


def load_script(filename):
    """Load script content from fitting_scripts directory."""
    script_dir = Path(__file__).parent.parent / 'fitting_scripts'
    script_path = script_dir / filename
    if script_path.exists():
        return script_path.read_text()
    return None


def sync_wang_script(apps, schema_editor):
    """Update Wang script and schemas with restrain_to_controls feature."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')

    wang_script = load_script('tight_binding_wang.py')
    if not wang_script:
        return  # Script not available, skip

    # Update the Wang fitting method with new script and updated schemas
    FittingMethod.objects.filter(slug='tight-binding-wang').update(
        script=wang_script,
        input_schema={
            'type': 'object',
            'properties': {
                'concentrations': {
                    'type': 'array',
                    'items': {'type': 'number'},
                    'description': 'Inhibitor concentrations in nM',
                },
                'responses': {
                    'type': 'array',
                    'items': {'type': 'number'},
                    'description': 'Response values (e.g., fluorescence, binding signal)',
                },
                'controls': {
                    'type': 'object',
                    'properties': {
                        'max': {'type': 'number', 'description': 'Uninhibited control'},
                        'min': {'type': 'number', 'description': 'Fully inhibited control'},
                    },
                },
                'parameters': {
                    'type': 'object',
                    'properties': {
                        'protein_conc': {
                            'type': 'number',
                            'description': 'Total protein concentration [P]t in nM',
                        },
                        'ligand_conc': {
                            'type': 'number',
                            'description': 'Total labeled ligand concentration [L]t in nM',
                        },
                        'ligand_kd': {
                            'type': 'number',
                            'description': 'Kd of labeled ligand in nM',
                        },
                        'fix_top': {
                            'type': ['boolean', 'number', 'null'],
                            'description': 'Fixed top asymptote (hard constraint). True = use max control.',
                        },
                        'fix_bottom': {
                            'type': ['boolean', 'number', 'null'],
                            'description': 'Fixed bottom asymptote (hard constraint). True = use min control.',
                        },
                        'restrain_to_controls': {
                            'type': 'boolean',
                            'description': 'Soft-constrain asymptotes via pseudo data points',
                            'default': False,
                        },
                        'pseudo_point_offset_logs': {
                            'type': 'number',
                            'description': 'Log units offset for pseudo points (default 3.0)',
                            'default': 3.0,
                        },
                    },
                    'required': ['protein_conc', 'ligand_conc', 'ligand_kd'],
                },
            },
            'required': ['concentrations', 'responses', 'parameters'],
        },
        output_schema={
            'type': 'object',
            'properties': {
                'ki': {'type': 'number', 'description': 'Fitted Ki value (primary output)'},
                'ic50_apparent': {'type': ['number', 'null'], 'description': 'Apparent IC50 for reference'},
                'top': {'type': 'number', 'description': 'Top asymptote'},
                'bottom': {'type': 'number', 'description': 'Bottom asymptote'},
                'r_squared': {'type': 'number', 'description': 'R-squared goodness of fit'},
                'curve_points': {'type': 'array', 'description': 'Points for plotting'},
                'flags': {'type': 'array', 'items': {'type': 'string'}},
                'kpi': {'type': 'string', 'const': 'ki'},
                'fit_successful': {'type': 'boolean'},
                'restraint_applied': {'type': 'boolean', 'description': 'Whether pseudo-point restraints were used'},
                'tight_binding_params': {
                    'type': 'object',
                    'properties': {
                        'protein_conc': {'type': 'number'},
                        'ligand_conc': {'type': 'number'},
                        'ligand_kd': {'type': 'number'},
                    },
                },
            },
        },
    )


def noop(apps, schema_editor):
    """No-op reverse migration - script will remain updated."""
    pass


class Migration(migrations.Migration):
    """
    Sync Wang tight-binding script with restrain_to_controls feature.

    Updates:
    - Script: Adds soft constraint support via pseudo data points
    - Input schema: Adds restrain_to_controls and pseudo_point_offset_logs parameters
    - Output schema: Adds restraint_applied field

    This feature matches the 4PL restrain_to_controls implementation, allowing
    the fit to be guided toward control values without hard-fixing asymptotes.
    """

    dependencies = [
        ('assays', '0023_ensure_wang_fitting_method'),
    ]

    operations = [
        migrations.RunPython(sync_wang_script, noop),
    ]
