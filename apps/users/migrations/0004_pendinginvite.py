# Generated migration for adding PendingInvite model

from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0003_userprofile_operating_level'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='PendingInvite',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('email', models.EmailField(max_length=320)),
                ('note', models.CharField(blank=True, max_length=255)),
                ('status', models.CharField(
                    choices=[('pending', 'Pending'), ('sent', 'Sent'), ('failed', 'Failed')],
                    default='pending',
                    max_length=10,
                )),
                ('failure_reason', models.CharField(blank=True, max_length=500)),
                ('requested_at', models.DateTimeField(auto_now_add=True)),
                ('sent_at', models.DateTimeField(blank=True, null=True)),
                ('guest_object_id', models.CharField(blank=True, max_length=64)),
                ('requested_by', models.ForeignKey(
                    blank=True, null=True,
                    on_delete=models.deletion.SET_NULL,
                    related_name='invites_requested',
                    to=settings.AUTH_USER_MODEL,
                )),
                ('sent_by', models.ForeignKey(
                    blank=True, null=True,
                    on_delete=models.deletion.SET_NULL,
                    related_name='invites_sent',
                    to=settings.AUTH_USER_MODEL,
                )),
            ],
            options={
                'ordering': ['-requested_at'],
                'indexes': [models.Index(fields=['status', 'requested_at'], name='users_pendi_status_req_idx')],
            },
        ),
    ]
