# Data migration to update 4PL fitting script with restrain_to_controls feature
# Adds soft constraint option using pseudo data points at extreme concentrations

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
    """Update 4PL fitting script with restrain_to_controls feature."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')

    # Update 4PL script
    four_pl_script = load_script('four_parameter_logistic.py')
    if four_pl_script:
        FittingMethod.objects.filter(slug='four-parameter-logistic').update(
            script=four_pl_script
        )


def noop(apps, schema_editor):
    """No-op reverse migration - scripts will remain updated."""
    pass


class Migration(migrations.Migration):
    """
    Update 4PL fitting script to add restrain_to_controls parameter.

    New parameters:
    - restrain_to_controls: boolean (default false) - enables soft constraints
    - pseudo_point_offset_logs: number (default 3.0) - log units offset

    When enabled, adds pseudo data points at extreme concentrations to guide
    the fit toward control values while allowing deviation if data suggests
    different asymptotes.
    """

    dependencies = [
        ('assays', '0020_add_dataseries_batch'),
    ]

    operations = [
        migrations.RunPython(update_fitting_scripts, noop),
    ]
