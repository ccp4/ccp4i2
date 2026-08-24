"""Materialise ProjectTag.path so the tag forest has real uniqueness.

``unique_together = ["parent", "text"]`` never constrained rows with
``parent = NULL`` — SQL treats NULLs as distinct — so duplicate top-level tags
were reachable through the importer and through direct ORM writes. Populating
the path has to repair that first: two nodes with the same ancestry are
indistinguishable to a user, so they are merged rather than renamed.
"""

import logging

from django.db import migrations, models

logger = logging.getLogger(f"ccp4i2:{__name__}")

PATH_SEPARATOR = "\x1f"


def populate_paths(apps, schema_editor):
    ProjectTag = apps.get_model("ccp4i2", "ProjectTag")

    tags = list(ProjectTag.objects.all())
    if not tags:
        return

    children_of = {}
    for tag in tags:
        children_of.setdefault(tag.parent_id, []).append(tag)

    by_id = {tag.id: tag for tag in tags}
    keeper_by_path = {}
    merged_into = {}
    to_update = []
    merges = 0

    level = children_of.get(None, [])
    while level:
        next_level = []
        for tag in level:
            parent_id = merged_into.get(tag.parent_id, tag.parent_id)
            if parent_id is None:
                path = tag.text
            else:
                path = by_id[parent_id].path + PATH_SEPARATOR + tag.text

            keeper_id = keeper_by_path.get(path)
            if keeper_id is not None and keeper_id != tag.id:
                # Same ancestry, same text: indistinguishable, so fold this row
                # into the one already holding the path.
                keeper = by_id[keeper_id]
                keeper.projects.add(*tag.projects.all())
                merged_into[tag.id] = keeper_id
                merges += 1
                next_level.extend(children_of.get(tag.id, []))
                tag.delete()
                continue

            tag.path = path
            keeper_by_path[path] = tag.id
            to_update.append(tag)
            next_level.extend(children_of.get(tag.id, []))
        level = next_level

    ProjectTag.objects.bulk_update(to_update, ["path"], batch_size=200)
    if merges:
        logger.warning(
            "Merged %s duplicate project tag(s) while materialising tag paths", merges
        )


def clear_paths(apps, schema_editor):
    ProjectTag = apps.get_model("ccp4i2", "ProjectTag")
    ProjectTag.objects.update(path="")


class Migration(migrations.Migration):

    dependencies = [
        ("ccp4i2", "0016_projectexport_status"),
    ]

    operations = [
        migrations.AddField(
            model_name="projecttag",
            name="path",
            field=models.CharField(default="", editable=False, max_length=1024),
        ),
        migrations.RunPython(populate_paths, clear_paths),
        migrations.AlterField(
            model_name="projecttag",
            name="path",
            field=models.CharField(
                default="", editable=False, max_length=1024, unique=True
            ),
        ),
    ]
