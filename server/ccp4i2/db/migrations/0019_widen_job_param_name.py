# A job parameter name is a dotted path into the task container, and the paths
# the New Project flow writes -- "ASUCONTENTFILE.fileContent.seqList[0].source",
# 44 characters -- passed the Qt-era limit of 32 long ago. SQLite never enforced
# it, so the live database held them happily; the snapshot restore, which goes
# through the model serializer, rejected them, and a database rebuilt from its
# DATABASE.db.xml safety net came back without those files (2026-09-04).
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("ccp4i2", "0018_fileimport_source_checksum"),
    ]

    operations = [
        migrations.AlterField(
            model_name="file",
            name="job_param_name",
            field=models.CharField(blank=True, max_length=255),
        ),
        migrations.AlterField(
            model_name="fileuse",
            name="job_param_name",
            field=models.CharField(blank=True, max_length=255),
        ),
    ]
