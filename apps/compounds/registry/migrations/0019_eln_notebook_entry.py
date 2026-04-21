"""
Add ELN-linked compound registration support.

This migration adds:
- LabNotebookEntry model for ELN pages and paper notebook pages
- CompoundDocument model for per-compound document attachments
- New fields on Compound: notebook_entry, helm_notation, sequence_display
"""

import uuid

import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models

import compounds.registry.models


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('registry', '0018_add_isomer_stereo_choices'),
    ]

    operations = [
        # Create LabNotebookEntry model
        migrations.CreateModel(
            name='LabNotebookEntry',
            fields=[
                ('id', models.UUIDField(
                    default=uuid.uuid4,
                    editable=False,
                    primary_key=True,
                    serialize=False
                )),
                ('sequence_number', models.IntegerField(
                    help_text='Per-supplier sequential index (e.g., 1 for KF001, 22 for KF022)'
                )),
                ('title', models.CharField(
                    blank=True,
                    help_text='Page title, e.g., "SA157 Cyclic peptide synthesis and characterisation"',
                    max_length=255
                )),
                ('date', models.DateField(
                    blank=True,
                    help_text='Date of the notebook entry',
                    null=True
                )),
                ('url', models.URLField(
                    blank=True,
                    help_text='Hyperlink to ELN page (OneNote, LabArchives, etc.)',
                    max_length=2048
                )),
                ('labbook_number', models.IntegerField(
                    blank=True,
                    help_text='Lab notebook number (for paper notebooks)',
                    null=True
                )),
                ('page_number', models.IntegerField(
                    blank=True,
                    help_text='Page number in lab notebook (for paper notebooks)',
                    null=True
                )),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('supplier', models.ForeignKey(
                    help_text='Supplier/chemist who owns this notebook entry',
                    on_delete=django.db.models.deletion.PROTECT,
                    related_name='notebook_entries',
                    to='registry.supplier'
                )),
                ('created_by', models.ForeignKey(
                    blank=True,
                    null=True,
                    on_delete=django.db.models.deletion.SET_NULL,
                    related_name='created_notebook_entries',
                    to=settings.AUTH_USER_MODEL
                )),
            ],
            options={
                'verbose_name': 'Lab Notebook Entry',
                'verbose_name_plural': 'Lab Notebook Entries',
                'ordering': ['supplier', 'sequence_number'],
                'unique_together': {('supplier', 'sequence_number')},
            },
        ),

        # Add new fields to Compound
        migrations.AddField(
            model_name='compound',
            name='notebook_entry',
            field=models.ForeignKey(
                blank=True,
                help_text='Lab notebook entry (ELN page or paper notebook page)',
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name='compounds',
                to='registry.labnotebookentry'
            ),
        ),
        migrations.AddField(
            model_name='compound',
            name='helm_notation',
            field=models.CharField(
                blank=True,
                help_text='Optional HELM string. Phase 1: opaque text, no validation.',
                max_length=4096
            ),
        ),
        migrations.AddField(
            model_name='compound',
            name='sequence_display',
            field=models.CharField(
                blank=True,
                help_text="Human-readable sequence with modifications shown inline, e.g., "
                          "'FAM-linker-MCDWDIYRFPNHHC(1,4-Xylene)-NH2'. Free text.",
                max_length=512
            ),
        ),

        # Create CompoundDocument model
        migrations.CreateModel(
            name='CompoundDocument',
            fields=[
                ('id', models.UUIDField(
                    default=uuid.uuid4,
                    editable=False,
                    primary_key=True,
                    serialize=False
                )),
                ('kind', models.CharField(
                    choices=[
                        ('chemdraw', 'Chemdraw structure'),
                        ('qc', 'QC report'),
                        ('spectrum', 'Spectrum / chromatogram'),
                        ('other', 'Other'),
                    ],
                    help_text='Type of document',
                    max_length=16
                )),
                ('label', models.CharField(
                    blank=True,
                    help_text="Display text, e.g., 'CP1.cdxml'",
                    max_length=255
                )),
                ('url', models.URLField(
                    blank=True,
                    help_text='External URL to document (SharePoint, cloud storage, etc.)',
                    max_length=2048
                )),
                ('file', models.FileField(
                    blank=True,
                    help_text='Uploaded document file',
                    max_length=255,
                    upload_to=compounds.registry.models._compound_document_path
                )),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('compound', models.ForeignKey(
                    on_delete=django.db.models.deletion.CASCADE,
                    related_name='documents',
                    to='registry.compound'
                )),
                ('created_by', models.ForeignKey(
                    blank=True,
                    null=True,
                    on_delete=django.db.models.deletion.SET_NULL,
                    related_name='created_compound_documents',
                    to=settings.AUTH_USER_MODEL
                )),
            ],
            options={
                'verbose_name': 'Compound Document',
                'verbose_name_plural': 'Compound Documents',
                'ordering': ['compound', 'kind', 'created_at'],
            },
        ),
    ]
