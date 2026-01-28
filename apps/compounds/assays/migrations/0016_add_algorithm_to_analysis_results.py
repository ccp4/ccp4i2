# Data migration to add algorithm identifier to existing AnalysisResult records
from django.db import migrations


def add_algorithm_field_forward(apps, schema_editor):
    """
    Add 'algorithm' field to existing AnalysisResult.results JSON.

    For existing results:
    - If tight_binding_params exists → 'tight-binding-wang'
    - If METHOD == 'Hill-Langmuir' (legacy format) → 'hill-langmuir-legacy'
    - Otherwise → 'four-parameter-logistic' (the current default)

    Legacy Hill-Langmuir results have a different parameter structure:
    - EC50 instead of IC50
    - endPercent field
    - METHOD field set to 'Hill-Langmuir'
    """
    AnalysisResult = apps.get_model('assays', 'AnalysisResult')

    for analysis in AnalysisResult.objects.all():
        results = analysis.results or {}

        # Skip if algorithm already set
        if 'algorithm' in results:
            continue

        # Infer algorithm from result content
        if results.get('tight_binding_params'):
            algorithm = 'tight-binding-wang'
        elif results.get('METHOD') == 'Hill-Langmuir':
            # Legacy format with EC50, endPercent, METHOD fields
            algorithm = 'hill-langmuir-legacy'
        else:
            algorithm = 'four-parameter-logistic'

        results['algorithm'] = algorithm
        analysis.results = results
        analysis.save(update_fields=['results'])


def add_algorithm_field_reverse(apps, schema_editor):
    """Remove 'algorithm' field from AnalysisResult.results JSON."""
    AnalysisResult = apps.get_model('assays', 'AnalysisResult')

    for analysis in AnalysisResult.objects.all():
        results = analysis.results or {}
        if 'algorithm' in results:
            del results['algorithm']
            analysis.results = results
            analysis.save(update_fields=['results'])


class Migration(migrations.Migration):
    """Data migration to add algorithm identifier to existing analysis results."""

    dependencies = [
        ('assays', '0015_finalize_plate_layout_fk'),
    ]

    operations = [
        migrations.RunPython(add_algorithm_field_forward, add_algorithm_field_reverse),
    ]
