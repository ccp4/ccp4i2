# Generated manually for audit trail fields

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('assays', '0018_add_import_type'),
    ]

    operations = [
        # Protocol audit fields
        migrations.AlterField(
            model_name='protocol',
            name='created_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='created_protocols',
                to=settings.AUTH_USER_MODEL
            ),
        ),
        migrations.AddField(
            model_name='protocol',
            name='modified_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='modified_protocols',
                to=settings.AUTH_USER_MODEL,
                help_text='User who last modified this protocol'
            ),
        ),
        migrations.AddField(
            model_name='protocol',
            name='modified_at',
            field=models.DateTimeField(auto_now=True, null=True),
        ),

        # Assay audit fields
        migrations.AlterField(
            model_name='assay',
            name='created_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='created_assays',
                to=settings.AUTH_USER_MODEL
            ),
        ),
        migrations.AddField(
            model_name='assay',
            name='modified_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='modified_assays',
                to=settings.AUTH_USER_MODEL,
                help_text='User who last modified this assay'
            ),
        ),
        migrations.AddField(
            model_name='assay',
            name='modified_at',
            field=models.DateTimeField(auto_now=True, null=True),
        ),
    ]
