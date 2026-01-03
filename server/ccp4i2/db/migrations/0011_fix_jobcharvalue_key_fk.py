# Generated manually to fix JobCharValue.key FK reference
# Migration 0006 updated JobFloatValue.key but not JobCharValue.key

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('ccp4i2', '0010_alter_projectgroup_projects'),
    ]

    operations = [
        migrations.AlterField(
            model_name='jobcharvalue',
            name='key',
            field=models.ForeignKey(on_delete=django.db.models.deletion.RESTRICT, related_name='+', to='ccp4i2.jobvaluekey'),
        ),
    ]
