"""
Backfill import_type='pharmaron_adme' for Protocols that were auto-created
by the ADME import endpoint before the views.py fix.

Such protocols have analysis_method='pharmaron_adme' but fell through to the
default import_type='raw_data', which caused the Protocol detail page to open
the wrong "Add Assay" sidebar.
"""

from django.db import migrations


def backfill_adme_import_type(apps, schema_editor):
    Protocol = apps.get_model('assays', 'Protocol')
    Protocol.objects.filter(
        analysis_method='pharmaron_adme',
    ).exclude(
        import_type='pharmaron_adme',
    ).update(import_type='pharmaron_adme')


def noop(apps, schema_editor):
    pass


class Migration(migrations.Migration):

    dependencies = [
        ('assays', '0024_sync_wang_restrain_to_controls'),
    ]

    operations = [
        migrations.RunPython(backfill_adme_import_type, noop),
    ]
