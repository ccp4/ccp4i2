# Data migration to sync fitting scripts to latest versions
# Updates both 4PL and Wang tight-binding scripts with unit-tested implementations

from django.db import migrations
from pathlib import Path


def load_script(filename):
    """Load script content from fitting_scripts directory."""
    script_dir = Path(__file__).parent.parent / 'fitting_scripts'
    script_path = script_dir / filename
    if script_path.exists():
        return script_path.read_text()
    return None


def update_fitting_scripts(apps, schema_editor):
    """Update fitting scripts to latest unit-tested versions."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')

    # Update 4PL script (includes restrain_to_controls feature)
    four_pl_script = load_script('four_parameter_logistic.py')
    if four_pl_script:
        FittingMethod.objects.filter(slug='four-parameter-logistic').update(
            script=four_pl_script
        )

    # Update Wang tight-binding script
    wang_script = load_script('tight_binding_wang.py')
    if wang_script:
        FittingMethod.objects.filter(slug='tight-binding-wang').update(
            script=wang_script
        )


def noop(apps, schema_editor):
    """No-op reverse migration - scripts will remain updated."""
    pass


class Migration(migrations.Migration):
    """
    Sync fitting scripts to latest unit-tested versions.

    Updates:
    - four_parameter_logistic.py: restrain_to_controls, curve_points, KPI handling
    - tight_binding_wang.py: curve_points generation for frontend chart rendering

    Both scripts now return curve_points in the result, enabling algorithm-agnostic
    chart rendering in the frontend.
    """

    dependencies = [
        ('assays', '0021_update_4pl_restrain_to_controls'),
    ]

    operations = [
        migrations.RunPython(update_fitting_scripts, noop),
    ]
