# Generated manually for audit trail fields

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('registry', '0013_compound_aliases'),
    ]

    operations = [
        # Supplier audit fields
        migrations.AddField(
            model_name='supplier',
            name='created_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='created_suppliers',
                to=settings.AUTH_USER_MODEL
            ),
        ),
        migrations.AddField(
            model_name='supplier',
            name='modified_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='modified_suppliers',
                to=settings.AUTH_USER_MODEL
            ),
        ),
        migrations.AddField(
            model_name='supplier',
            name='created_at',
            field=models.DateTimeField(auto_now_add=True, null=True),
        ),
        migrations.AddField(
            model_name='supplier',
            name='modified_at',
            field=models.DateTimeField(auto_now=True, null=True),
        ),

        # Target audit fields
        migrations.AddField(
            model_name='target',
            name='created_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='created_targets',
                to=settings.AUTH_USER_MODEL
            ),
        ),
        migrations.AddField(
            model_name='target',
            name='modified_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='modified_targets',
                to=settings.AUTH_USER_MODEL
            ),
        ),
        migrations.AddField(
            model_name='target',
            name='modified_at',
            field=models.DateTimeField(auto_now=True, null=True),
        ),

        # Compound modified_by field
        migrations.AddField(
            model_name='compound',
            name='modified_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='modified_compounds',
                to=settings.AUTH_USER_MODEL,
                help_text='User who last modified this compound'
            ),
        ),

        # Batch audit fields
        migrations.AddField(
            model_name='batch',
            name='registered_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='registered_batches',
                to=settings.AUTH_USER_MODEL,
                help_text='User who registered this batch'
            ),
        ),
        migrations.AddField(
            model_name='batch',
            name='modified_by',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='modified_batches',
                to=settings.AUTH_USER_MODEL,
                help_text='User who last modified this batch'
            ),
        ),
        migrations.AddField(
            model_name='batch',
            name='modified_at',
            field=models.DateTimeField(auto_now=True, null=True),
        ),
    ]
