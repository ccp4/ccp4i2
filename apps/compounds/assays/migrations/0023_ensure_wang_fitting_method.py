# Data migration to ensure the Wang tight-binding fitting method exists
# This fixes cases where migration 0006 ran before the script file was available

from django.db import migrations
from pathlib import Path


def load_script(filename):
    """Load script content from fitting_scripts directory."""
    script_dir = Path(__file__).parent.parent / 'fitting_scripts'
    script_path = script_dir / filename
    if script_path.exists():
        return script_path.read_text()
    return None


def ensure_wang_fitting_method(apps, schema_editor):
    """Ensure the Wang tight-binding fitting method exists."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')

    wang_script = load_script('tight_binding_wang.py')
    if not wang_script:
        return  # Script not available, skip

    # Use update_or_create to create if missing or update if exists
    FittingMethod.objects.update_or_create(
        slug='tight-binding-wang',
        defaults={
            'version': '1.0.0',
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

Constraint options:
- fix_top: Hard-fix the top asymptote to a specific value
- fix_bottom: Hard-fix the bottom asymptote to a specific value
- restrain_to_controls: Soft-constrain asymptotes using control values as
  pseudo data points. Guides fit toward controls while allowing deviation
  if the data strongly suggests different asymptotes.''',
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
                                'description': 'Fixed top asymptote (hard constraint)',
                            },
                            'fix_bottom': {
                                'type': ['number', 'null'],
                                'description': 'Fixed bottom asymptote (hard constraint)',
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
            'is_active': True,
            'is_builtin': True,
        }
    )


def noop(apps, schema_editor):
    """No-op reverse migration."""
    pass


class Migration(migrations.Migration):
    """
    Ensure the Wang tight-binding fitting method exists.

    This migration fixes cases where the original migration 0006 ran but
    the script file wasn't available (e.g., in certain Docker build scenarios).
    It uses update_or_create to ensure the method exists with the correct data.
    """

    dependencies = [
        ('assays', '0022_sync_fitting_scripts'),
    ]

    operations = [
        migrations.RunPython(ensure_wang_fitting_method, noop),
    ]
