"""
Remove the pherastar_table field from Protocol model.

This field was a legacy PHERAstar plate reader configuration that is no longer used.
"""

from django.db import migrations


class Migration(migrations.Migration):
    """Remove pherastar_table field from Protocol."""

    dependencies = [
        ('assays', '0009_add_adme_analysis_method'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='protocol',
            name='pherastar_table',
        ),
    ]
