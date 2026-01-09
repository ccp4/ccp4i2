# Data migration to seed built-in fitting methods

from django.db import migrations
from pathlib import Path


def load_script(filename):
    """Load script content from fitting_scripts directory."""
    script_dir = Path(__file__).parent.parent / 'fitting_scripts'
    script_path = script_dir / filename
    if script_path.exists():
        return script_path.read_text()
    return None


def seed_fitting_methods(apps, schema_editor):
    """Seed the database with built-in fitting methods."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')

    # Four Parameter Logistic (4PL)
    four_pl_script = load_script('four_parameter_logistic.py')
    if four_pl_script:
        FittingMethod.objects.update_or_create(
            slug='four-parameter-logistic',
            version='1.0.0',
            defaults={
                'name': 'Four Parameter Logistic (4PL)',
                'description': '''Standard sigmoidal dose-response curve fitting using the Hill-Langmuir equation.

The 4PL model fits data to: y = bottom + (top - bottom) / (1 + (x / IC50)^hill)

Parameters:
- top: Maximum response (asymptote at high concentrations)
- bottom: Minimum response (asymptote at low concentrations)
- IC50: Concentration at 50% response
- hill: Hill coefficient (curve steepness)

Options:
- fix_hill: Set to a number to fix the Hill coefficient
- fix_top: Set to a number to fix the top asymptote
- fix_bottom: Set to a number to fix the bottom asymptote''',
                'script': four_pl_script,
                'input_schema': {
                    'type': 'object',
                    'properties': {
                        'concentrations': {
                            'type': 'array',
                            'items': {'type': 'number'},
                            'description': 'Concentration values in order from high to low',
                        },
                        'responses': {
                            'type': 'array',
                            'items': {'type': 'number'},
                            'description': 'Response values corresponding to concentrations',
                        },
                        'controls': {
                            'type': 'object',
                            'properties': {
                                'max': {'type': 'number', 'description': 'Max control value'},
                                'min': {'type': 'number', 'description': 'Min control value'},
                            },
                        },
                        'parameters': {
                            'type': 'object',
                            'properties': {
                                'fix_hill': {'type': ['number', 'null'], 'description': 'Fixed Hill coefficient'},
                                'fix_top': {'type': ['number', 'null'], 'description': 'Fixed top asymptote'},
                                'fix_bottom': {'type': ['number', 'null'], 'description': 'Fixed bottom asymptote'},
                            },
                        },
                    },
                    'required': ['concentrations', 'responses'],
                },
                'output_schema': {
                    'type': 'object',
                    'properties': {
                        'ic50': {'type': 'number', 'description': 'IC50 value'},
                        'hill_slope': {'type': 'number', 'description': 'Hill coefficient'},
                        'top': {'type': 'number', 'description': 'Top asymptote'},
                        'bottom': {'type': 'number', 'description': 'Bottom asymptote'},
                        'r_squared': {'type': 'number', 'description': 'R-squared goodness of fit'},
                        'curve_points': {'type': 'array', 'description': 'Points for plotting the fitted curve'},
                        'flags': {'type': 'array', 'items': {'type': 'string'}, 'description': 'Quality flags'},
                        'kpi': {'type': 'string', 'description': 'Which field is the primary KPI'},
                        'fit_successful': {'type': 'boolean', 'description': 'Whether fitting succeeded'},
                    },
                },
                'is_active': True,
                'is_builtin': True,
            }
        )


def remove_fitting_methods(apps, schema_editor):
    """Remove built-in fitting methods (for reverse migration)."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')
    FittingMethod.objects.filter(is_builtin=True).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('assays', '0004_add_fitting_method_and_plate_layout'),
    ]

    operations = [
        migrations.RunPython(seed_fitting_methods, remove_fitting_methods),
    ]
