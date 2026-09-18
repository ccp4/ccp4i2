from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ccp4i2', '0020_add_dnatco_naval_json_file_type'),
    ]

    operations = [
        migrations.AddField(
            model_name='fileimport',
            name='description',
            field=models.TextField(blank=True, default=''),
        ),
    ]
