from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("ccp4i2", "0015_add_phasertng_dag_file_type"),
    ]

    operations = [
        migrations.AddField(
            model_name="projectexport",
            name="status",
            field=models.CharField(
                choices=[
                    ("pending", "Pending"),
                    ("running", "Running"),
                    ("completed", "Completed"),
                    ("failed", "Failed"),
                ],
                default="pending",
                max_length=16,
            ),
        ),
    ]
