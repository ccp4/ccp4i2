# Generated manually for fitting method and plate layout models

import uuid
from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('assays', '0003_alter_protocol_pherastar_table'),
    ]

    operations = [
        # Create FittingMethod model
        migrations.CreateModel(
            name='FittingMethod',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('name', models.CharField(help_text='Display name for the fitting method', max_length=128)),
                ('slug', models.SlugField(help_text="URL-safe identifier, e.g., 'four-parameter-logistic'", max_length=64, unique=True)),
                ('version', models.CharField(default='1.0.0', help_text='Semantic version string', max_length=32)),
                ('description', models.TextField(blank=True, help_text='Detailed description of the fitting algorithm', null=True)),
                ('script', models.TextField(help_text='Python script implementing the fit() function')),
                ('input_schema', models.JSONField(blank=True, default=dict, help_text='JSON Schema describing expected input parameters')),
                ('output_schema', models.JSONField(blank=True, default=dict, help_text='JSON Schema describing output structure')),
                ('is_active', models.BooleanField(default=True, help_text='Inactive methods are hidden from selection but preserved for historical data')),
                ('is_builtin', models.BooleanField(default=False, help_text='Built-in methods shipped with the system')),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_fitting_methods', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Fitting Method',
                'verbose_name_plural': 'Fitting Methods',
                'ordering': ['name', '-version'],
                'unique_together': {('slug', 'version')},
            },
        ),
        # Add plate_layout field to Protocol
        migrations.AddField(
            model_name='protocol',
            name='plate_layout',
            field=models.JSONField(blank=True, default=dict, help_text='Plate layout template: control positions, sample region, replicate pattern'),
        ),
        # Add fitting_parameters field to Protocol
        migrations.AddField(
            model_name='protocol',
            name='fitting_parameters',
            field=models.JSONField(blank=True, default=dict, help_text='Default parameters passed to fitting script'),
        ),
        # Add fitting_method FK to Protocol
        migrations.AddField(
            model_name='protocol',
            name='fitting_method',
            field=models.ForeignKey(blank=True, help_text='Curve fitting algorithm for analyzing data from this protocol', null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='protocols', to='assays.fittingmethod'),
        ),
        # Update analysis_method help_text (keeping for backward compatibility)
        migrations.AlterField(
            model_name='protocol',
            name='analysis_method',
            field=models.CharField(choices=[('hill_langmuir', 'Hill-Langmuir'), ('hill_langmuir_fix_hill', 'Hill-Langmuir (fixed Hill coefficient)'), ('hill_langmuir_fix_hill_minmax', 'Hill-Langmuir (fixed Hill and min/max)'), ('hill_langmuir_fix_minmax', 'Hill-Langmuir (fixed min/max)'), ('ms_intact', 'MS-Intact'), ('table_of_values', 'Table of values')], default='hill_langmuir', help_text='Legacy analysis method (use fitting_method for new protocols)', max_length=50),
        ),
    ]
