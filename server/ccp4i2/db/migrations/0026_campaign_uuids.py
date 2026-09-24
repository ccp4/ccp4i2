"""Give campaigns and their sites a uuid.

The recovery format keys on uuids, and a site is referred to from other rows
(evaluations, and soon event records): an integer primary key is meaningless
after a database rebuild and a name is the renameable thing whose instability
motivated migration 0024. Until now neither ProjectGroup nor CampaignSite had
one, which is why no campaign state could be snapshotted -- lose db.sqlite3
and every verdict in every campaign went with it.

Three steps, not one: an AddField with ``default=uuid4`` evaluates the
callable once and stamps the same value onto every existing row, which would
then fail the unique constraint. So the column is added nullable, each row is
given its own uuid, and only then is it made unique and non-null.
"""
from uuid import uuid4

from django.db import migrations, models


def assign_uuids(apps, schema_editor):
    for name in ("ProjectGroup", "CampaignSite"):
        model = apps.get_model("ccp4i2", name)
        for row in model.objects.filter(uuid__isnull=True).iterator():
            row.uuid = uuid4()
            row.save(update_fields=["uuid"])


class Migration(migrations.Migration):

    dependencies = [
        ("ccp4i2", "0025_site_tags_become_evaluations"),
    ]

    operations = [
        migrations.AddField(
            model_name="projectgroup",
            name="uuid",
            field=models.UUIDField(null=True, editable=False),
        ),
        migrations.AddField(
            model_name="campaignsite",
            name="uuid",
            field=models.UUIDField(null=True, editable=False),
        ),
        migrations.RunPython(assign_uuids, migrations.RunPython.noop),
        migrations.AlterField(
            model_name="projectgroup",
            name="uuid",
            field=models.UUIDField(default=uuid4, unique=True, editable=False),
        ),
        migrations.AlterField(
            model_name="campaignsite",
            name="uuid",
            field=models.UUIDField(default=uuid4, unique=True, editable=False),
        ),
    ]
