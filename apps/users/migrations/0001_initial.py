# Generated migration for users app

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='UserProfile',
            fields=[
                ('user', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, related_name='profile', serialize=False, to=settings.AUTH_USER_MODEL)),
                ('is_platform_admin', models.BooleanField(default=False, help_text='Can manage platform settings, import data, and administer users')),
                ('legacy_username', models.CharField(blank=True, help_text='Username from legacy system (for reference)', max_length=100)),
                ('legacy_display_name', models.CharField(blank=True, help_text='Display name from legacy system', max_length=255)),
                ('imported_at', models.DateTimeField(blank=True, help_text='When this user was imported from legacy fixtures', null=True)),
                ('first_login_at', models.DateTimeField(blank=True, help_text='When user first logged in via Azure AD', null=True)),
                ('last_seen_at', models.DateTimeField(blank=True, help_text='Last activity timestamp', null=True)),
            ],
            options={
                'verbose_name': 'User Profile',
                'verbose_name_plural': 'User Profiles',
            },
        ),
    ]
