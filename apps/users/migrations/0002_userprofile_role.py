# Generated migration for adding role field to UserProfile

from django.db import migrations, models


def migrate_admin_to_role(apps, schema_editor):
    """
    Migrate existing users based on is_platform_admin flag:
    - is_platform_admin=True -> role='admin'
    - is_platform_admin=False -> role='contributor' (existing users can contribute)
    """
    UserProfile = apps.get_model('users', 'UserProfile')

    # Set admins
    UserProfile.objects.filter(is_platform_admin=True).update(role='admin')

    # Set non-admins as contributors (they were previously allowed to edit)
    UserProfile.objects.filter(is_platform_admin=False).update(role='contributor')


def reverse_migration(apps, schema_editor):
    """Reverse migration - sync is_platform_admin from role."""
    UserProfile = apps.get_model('users', 'UserProfile')
    UserProfile.objects.filter(role='admin').update(is_platform_admin=True)
    UserProfile.objects.exclude(role='admin').update(is_platform_admin=False)


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='userprofile',
            name='role',
            field=models.CharField(
                choices=[
                    ('user', 'User (read-only)'),
                    ('contributor', 'Contributor (can add/edit/delete)'),
                    ('admin', 'Admin (full access)'),
                ],
                default='contributor',
                help_text='Maximum authorization level for this user',
                max_length=20,
            ),
        ),
        migrations.RunPython(migrate_admin_to_role, reverse_migration),
    ]
