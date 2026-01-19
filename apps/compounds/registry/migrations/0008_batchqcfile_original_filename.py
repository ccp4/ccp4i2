# Generated manually for adding original_filename to BatchQCFile

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0007_add_user_to_supplier'),
    ]

    operations = [
        migrations.AddField(
            model_name='batchqcfile',
            name='original_filename',
            field=models.CharField(
                blank=True,
                help_text='Original filename as uploaded by user',
                max_length=255,
                null=True,
            ),
        ),
    ]
