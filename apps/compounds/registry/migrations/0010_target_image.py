# Generated manually for target dashboard feature

from django.db import migrations, models
import compounds.registry.models


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0009_populate_original_filename'),
    ]

    operations = [
        migrations.AddField(
            model_name='target',
            name='image',
            field=models.ImageField(
                blank=True,
                help_text='Branding image for the target dashboard',
                null=True,
                upload_to=compounds.registry.models._target_image_path,
            ),
        ),
    ]
