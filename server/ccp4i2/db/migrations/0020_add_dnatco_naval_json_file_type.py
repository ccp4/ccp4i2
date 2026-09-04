from django.db import migrations


def add_dnatco_naval_json_file_type(apps, schema_editor):
    FileType = apps.get_model("ccp4i2", "FileType")
    FileType.objects.get_or_create(
        name="application/dnatco-naval-json",
        defaults={"description": "DNATCO NAVAL nucleic acid geometry validation"},
    )


def remove_dnatco_naval_json_file_type(apps, schema_editor):
    FileType = apps.get_model("ccp4i2", "FileType")
    FileType.objects.filter(name="application/dnatco-naval-json").delete()


class Migration(migrations.Migration):

    dependencies = [
        ("ccp4i2", "0019_widen_job_param_name"),
    ]

    operations = [
        migrations.RunPython(add_dnatco_naval_json_file_type, remove_dnatco_naval_json_file_type),
    ]
