# Generated manually for adding compound aliases field

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0012_target_saved_aggregation_view'),
    ]

    operations = [
        migrations.AddField(
            model_name='compound',
            name='aliases',
            field=models.JSONField(
                blank=True,
                default=list,
                help_text="Alternative names/identifiers for this compound (e.g., supplier codes, "
                          "abbreviations, internal project names). Used as fallback during import matching."
            ),
        ),
    ]
