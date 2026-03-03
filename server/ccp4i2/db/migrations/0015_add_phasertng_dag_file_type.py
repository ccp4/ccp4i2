from django.db import migrations


def add_phasertng_dag_file_type(apps, schema_editor):
    FileType = apps.get_model("ccp4i2", "FileType")
    FileType.objects.get_or_create(
        name="application/phasertng-dag",
        defaults={"description": "PhaserTNG DAG file"},
    )


def remove_phasertng_dag_file_type(apps, schema_editor):
    FileType = apps.get_model("ccp4i2", "FileType")
    FileType.objects.filter(name="application/phasertng-dag").delete()


class Migration(migrations.Migration):

    dependencies = [
        ("ccp4i2", "0014_add_textdatafile_type"),
    ]

    operations = [
        migrations.RunPython(add_phasertng_dag_file_type, remove_phasertng_dag_file_type),
    ]
