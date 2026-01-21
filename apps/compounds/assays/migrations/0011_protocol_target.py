# Generated migration for adding target to Protocol

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0001_initial'),
        ('assays', '0010_remove_pherastar_table'),
    ]

    operations = [
        migrations.AddField(
            model_name='protocol',
            name='target',
            field=models.ForeignKey(
                blank=True,
                help_text='Default target for assays using this protocol',
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='protocols',
                to='registry.target',
            ),
        ),
    ]
