from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone
import uuid


class Migration(migrations.Migration):
    initial = True

    dependencies = [
        ("ccp4i2", "0012_fix_fk_constraints_sqlite"),
    ]

    operations = [
        migrations.CreateModel(
            name="StagedUpload",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                (
                    "uuid",
                    models.UUIDField(default=uuid.uuid4, editable=False, unique=True),
                ),
                (
                    "upload_type",
                    models.CharField(
                        choices=[
                            ("project_import", "Project Import"),
                            ("unmerged_data", "Unmerged Data"),
                        ],
                        max_length=32,
                    ),
                ),
                (
                    "status",
                    models.CharField(
                        choices=[
                            ("pending", "Pending upload"),
                            ("uploaded", "Uploaded, awaiting processing"),
                            ("processing", "Processing"),
                            ("completed", "Completed"),
                            ("failed", "Failed"),
                            ("expired", "Expired"),
                        ],
                        default="pending",
                        max_length=32,
                    ),
                ),
                ("original_filename", models.CharField(max_length=255)),
                ("blob_path", models.CharField(max_length=512)),
                ("sas_expiry", models.DateTimeField()),
                ("created_at", models.DateTimeField(default=django.utils.timezone.now)),
                ("completed_at", models.DateTimeField(blank=True, null=True)),
                ("error_message", models.TextField(blank=True)),
                ("requested_by", models.CharField(blank=True, max_length=255)),
                (
                    "target_job",
                    models.ForeignKey(
                        blank=True,
                        null=True,
                        on_delete=django.db.models.deletion.SET_NULL,
                        related_name="staged_uploads",
                        to="ccp4i2.job",
                    ),
                ),
            ],
            options={
                "ordering": ["-created_at"],
            },
        ),
    ]
