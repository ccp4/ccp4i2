# Generated manually

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('azure_extensions', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='stagedupload',
            name='error_message',
            field=models.TextField(blank=True, null=True),
        ),
    ]
