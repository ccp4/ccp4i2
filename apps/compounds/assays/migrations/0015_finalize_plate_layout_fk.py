# Finalize PlateLayout FK: remove JSONField and rename FK field
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):
    """
    Finalize the plate_layout field migration:
    1. Remove the old JSONField
    2. Rename the new FK field from plate_layout_new to plate_layout
    """

    dependencies = [
        ('assays', '0014_migrate_plate_layout_data'),
    ]

    operations = [
        # Step 1: Remove the old JSONField
        migrations.RemoveField(
            model_name='protocol',
            name='plate_layout',
        ),
        # Step 2: Rename the FK field from plate_layout_new to plate_layout
        migrations.RenameField(
            model_name='protocol',
            old_name='plate_layout_new',
            new_name='plate_layout',
        ),
        # Step 3: Update the related_name
        migrations.AlterField(
            model_name='protocol',
            name='plate_layout',
            field=models.ForeignKey(
                blank=True,
                help_text='Plate layout template: control positions, sample region, replicate pattern',
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='protocols',
                to='assays.platelayout'
            ),
        ),
    ]
