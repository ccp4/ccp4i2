from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('ccp4i2', '0021_fileimport_description'),
    ]

    operations = [
        migrations.CreateModel(
            name='JobInteractiveSession',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('requested_at', models.DateTimeField(default=django.utils.timezone.now)),
                ('last_heartbeat', models.DateTimeField(blank=True, null=True)),
                ('dispatched', models.BooleanField(default=False)),
                ('finished', models.BooleanField(default=False)),
                ('finished_at', models.DateTimeField(blank=True, null=True)),
                ('job', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='interactive_session', to='ccp4i2.job')),
            ],
        ),
    ]
