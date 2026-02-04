# Migration to add sites JSONField to ProjectGroup for campaign binding site navigation

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ccp4i2', '0012_fix_fk_constraints_sqlite'),
    ]

    operations = [
        migrations.AddField(
            model_name='projectgroup',
            name='sites',
            field=models.JSONField(blank=True, default=list),
        ),
    ]
