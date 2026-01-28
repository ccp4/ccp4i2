# Data migration to update built-in fitting method scripts
# This re-syncs the database scripts from the fitting_scripts directory

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
    """Update built-in fitting method scripts from files."""
    FittingMethod = apps.get_model('assays', 'FittingMethod')

    # Update 4PL script
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
    """Update fitting scripts to include end_percent and other improvements."""

    dependencies = [
        ('assays', '0016_add_algorithm_to_analysis_results'),
    ]

    operations = [
        migrations.RunPython(update_fitting_scripts, noop),
    ]
