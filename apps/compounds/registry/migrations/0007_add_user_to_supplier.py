# Generated manually for adding user link to Supplier

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('registry', '0006_alter_batchqcfile_file_max_length'),
    ]

    operations = [
        migrations.AddField(
            model_name='supplier',
            name='user',
            field=models.OneToOneField(
                blank=True,
                help_text='User account linked to this supplier (for personal suppliers)',
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='supplier',
                to=settings.AUTH_USER_MODEL,
            ),
        ),
    ]
