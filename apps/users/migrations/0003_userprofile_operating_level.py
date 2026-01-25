# Generated migration for adding operating_level field to UserProfile

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0002_userprofile_role'),
    ]

    operations = [
        migrations.AddField(
            model_name='userprofile',
            name='operating_level',
            field=models.CharField(
                blank=True,
                choices=[
                    ('user', 'User (read-only)'),
                    ('contributor', 'Contributor (can add/edit/delete)'),
                    ('admin', 'Admin (full access)'),
                ],
                help_text='Current operating level (defaults to role if not set)',
                max_length=20,
                null=True,
            ),
        ),
    ]
