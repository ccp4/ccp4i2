# Generated manually for DataSeries batch tracking

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0014_add_audit_fields'),
        ('assays', '0019_add_audit_fields'),
    ]

    operations = [
        migrations.AddField(
            model_name='dataseries',
            name='batch',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='assay_results',
                to='registry.batch',
                help_text='Specific batch tested (if known). Parsed from compound_name if '
                          'format is COMPOUND_ID/BATCH_NUMBER (e.g., NCL-00026042/1)'
            ),
        ),
    ]
