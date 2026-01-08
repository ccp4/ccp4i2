# Data migration to add tight-binding Wang equation fitting method

from django.db import migrations
from pathlib import Path


def load_script(filename):
    """Load script content from fitting_scripts directory."""
    script_dir = Path(__file__).parent.parent / 'fitting_scripts'
    script_path = script_dir / filename
    if script_path.exists():
        return script_path.read_text()
    return None


def seed_wang_fitting_method(apps, schema_editor):
    """Add the Wang tight-binding fitting method."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')

    wang_script = load_script('tight_binding_wang.py')
    if wang_script:
        FittingMethod.objects.update_or_create(
            slug='tight-binding-wang',
            version='1.0.0',
            defaults={
                'name': 'Tight-Binding Competition (Wang)',
                'description': '''Competitive tight-binding analysis using the Wang equation.

For assays where inhibitor concentration approaches protein concentration,
standard IC50 analysis gives incorrect Ki values. The Wang equation solves
the mass balance equations exactly to determine true Ki.

Required protocol parameters (set in fitting_parameters):
- protein_conc: Total protein concentration [P]t (nM)
- ligand_conc: Total labeled ligand concentration [L]t (nM)
- ligand_kd: Dissociation constant of labeled ligand Kd_L (nM)

Primary output: Ki (inhibition constant)

Reference: Wang, Z-X. (1995) FEBS Letters 360, 111-114.

Options:
- fix_top: Set to fix the top asymptote
- fix_bottom: Set to fix the bottom asymptote''',
                'script': wang_script,
                'input_schema': {
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
                                    'type': ['number', 'null'],
                                    'description': 'Fixed top asymptote',
                                },
                                'fix_bottom': {
                                    'type': ['number', 'null'],
                                    'description': 'Fixed bottom asymptote',
                                },
                            },
                            'required': ['protein_conc', 'ligand_conc', 'ligand_kd'],
                        },
                    },
                    'required': ['concentrations', 'responses', 'parameters'],
                },
                'output_schema': {
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
                'is_active': True,
                'is_builtin': True,
            }
        )


def remove_wang_fitting_method(apps, schema_editor):
    """Remove the Wang fitting method (for reverse migration)."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')
    FittingMethod.objects.filter(slug='tight-binding-wang').delete()


class Migration(migrations.Migration):

    dependencies = [
        ('assays', '0005_seed_builtin_fitting_methods'),
    ]

    operations = [
        migrations.RunPython(seed_wang_fitting_method, remove_wang_fitting_method),
    ]
