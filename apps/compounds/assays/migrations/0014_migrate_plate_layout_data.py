# Data migration to create PlateLayout records from existing Protocol.plate_layout JSON
import json
from django.db import migrations


def migrate_plate_layouts_forward(apps, schema_editor):
    """Create PlateLayout records from existing Protocol.plate_layout JSON and link them."""
    Protocol = apps.get_model('assays', 'Protocol')
    PlateLayout = apps.get_model('assays', 'PlateLayout')

    # Track layouts by JSON content to deduplicate
    layouts_seen = {}

    for protocol in Protocol.objects.all():
        # Skip protocols without plate_layout config
        plate_layout_data = protocol.plate_layout
        if not plate_layout_data or plate_layout_data == {}:
            continue

        # Create a key from sorted JSON for deduplication
        config_key = json.dumps(plate_layout_data, sort_keys=True)

        if config_key not in layouts_seen:
            # Create new PlateLayout record
            # Generate unique name based on protocol name
            base_name = f"Layout from {protocol.name[:90]}"
            name = base_name
            counter = 1
            while PlateLayout.objects.filter(name=name).exists():
                name = f"{base_name} ({counter})"
                counter += 1

            layout = PlateLayout.objects.create(
                name=name,
                config=plate_layout_data,
                created_by=protocol.created_by,
            )
            layouts_seen[config_key] = layout

        # Link protocol to the PlateLayout via the new FK field
        protocol.plate_layout_new = layouts_seen[config_key]
        protocol.save(update_fields=['plate_layout_new'])


def migrate_plate_layouts_reverse(apps, schema_editor):
    """Restore JSONField data from PlateLayout records."""
    Protocol = apps.get_model('assays', 'Protocol')

    for protocol in Protocol.objects.select_related('plate_layout_new').all():
        if protocol.plate_layout_new:
            protocol.plate_layout = protocol.plate_layout_new.config
            protocol.save(update_fields=['plate_layout'])


class Migration(migrations.Migration):
    """Data migration to populate PlateLayout records from existing Protocol JSON."""

    dependencies = [
        ('assays', '0013_add_plate_layout_model'),
    ]

    operations = [
        migrations.RunPython(migrate_plate_layouts_forward, migrate_plate_layouts_reverse),
    ]
