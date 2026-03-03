from django.db import migrations


def add_text_file_type(apps, schema_editor):
    FileType = apps.get_model("ccp4i2", "FileType")
    FileType.objects.get_or_create(
        name="text/plain",
        defaults={"description": "Plain text file"},
    )


def remove_text_file_type(apps, schema_editor):
    FileType = apps.get_model("ccp4i2", "FileType")
    FileType.objects.filter(name="text/plain").delete()


class Migration(migrations.Migration):

    dependencies = [
        ("ccp4i2", "0013_projectgroup_sites"),
    ]

    operations = [
        migrations.RunPython(add_text_file_type, remove_text_file_type),
    ]
