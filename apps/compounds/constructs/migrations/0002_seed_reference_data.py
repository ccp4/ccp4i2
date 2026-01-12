# Seed reference data for constructs app
from django.db import migrations


def seed_expression_tag_types(apps, schema_editor):
    """Seed ExpressionTagType with standard expression tags."""
    ExpressionTagType = apps.get_model('constructs', 'ExpressionTagType')
    tag_types = [
        'GST', 'SUMO', 'His6', 'His8', 'His10',
        'Avi', 'Strep', 'SBP', 'Flag', 'CFP', 'MBP'
    ]
    for name in tag_types:
        ExpressionTagType.objects.get_or_create(name=name)


def seed_proteases(apps, schema_editor):
    """Seed Protease with standard proteases."""
    Protease = apps.get_model('constructs', 'Protease')
    proteases = ['TEV', '3C', 'Thrombin']
    for name in proteases:
        Protease.objects.get_or_create(name=name)


def seed_expression_tag_locations(apps, schema_editor):
    """Seed ExpressionTagLocation with standard locations."""
    ExpressionTagLocation = apps.get_model('constructs', 'ExpressionTagLocation')
    locations = ['N-term', 'C-term', 'Internal']
    for name in locations:
        ExpressionTagLocation.objects.get_or_create(name=name)


def reverse_seed(apps, schema_editor):
    """Reverse migration - remove seeded data."""
    ExpressionTagType = apps.get_model('constructs', 'ExpressionTagType')
    Protease = apps.get_model('constructs', 'Protease')
    ExpressionTagLocation = apps.get_model('constructs', 'ExpressionTagLocation')

    # Only delete records that match our seeded values
    ExpressionTagType.objects.filter(name__in=[
        'GST', 'SUMO', 'His6', 'His8', 'His10',
        'Avi', 'Strep', 'SBP', 'Flag', 'CFP', 'MBP'
    ]).delete()

    Protease.objects.filter(name__in=['TEV', '3C', 'Thrombin']).delete()

    ExpressionTagLocation.objects.filter(
        name__in=['N-term', 'C-term', 'Internal']
    ).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('constructs', '0001_initial'),
    ]

    operations = [
        migrations.RunPython(seed_expression_tag_types, reverse_seed),
        migrations.RunPython(seed_proteases, reverse_seed),
        migrations.RunPython(seed_expression_tag_locations, reverse_seed),
    ]
