# Generated manually - initial schema for constructs app
from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion
import uuid


def _next_ncn_id():
    """Placeholder - actual implementation is in models.py"""
    return 1


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        # Reference data models
        migrations.CreateModel(
            name='ExpressionTagType',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=256, unique=True)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_expression_tag_types', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Expression Tag Type',
                'verbose_name_plural': 'Expression Tag Types',
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='Protease',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=256, unique=True)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_proteases', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Protease',
                'verbose_name_plural': 'Proteases',
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='ExpressionTagLocation',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=256, unique=True)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_expression_tag_locations', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Expression Tag Location',
                'verbose_name_plural': 'Expression Tag Locations',
                'ordering': ['name'],
            },
        ),
        # ConstructProject
        migrations.CreateModel(
            name='ConstructProject',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=256)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('parent', models.ForeignKey(blank=True, help_text='Parent project for hierarchical organization', null=True, on_delete=django.db.models.deletion.CASCADE, related_name='children', to='constructs.constructproject')),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_construct_projects', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Construct Project',
                'verbose_name_plural': 'Construct Projects',
                'ordering': ['name'],
            },
        ),
        # Protein
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('uniprot_id', models.CharField(help_text='UniProt identifier (e.g., P04637)', max_length=256, unique=True)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_proteins', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Protein',
                'verbose_name_plural': 'Proteins',
                'ordering': ['uniprot_id'],
            },
        ),
        # Plasmid
        migrations.CreateModel(
            name='Plasmid',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('ncn_id', models.IntegerField(editable=False, help_text='Auto-assigned construct number', unique=True)),
                ('name', models.CharField(help_text='Descriptive name for the plasmid', max_length=256, unique=True)),
                ('genbank_file', models.FileField(blank=True, help_text='GenBank or SnapGene file (.gb, .gbk)', max_length=500, null=True, upload_to='constructs/%(formatted_id)s/')),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('project', models.ForeignKey(blank=True, help_text='Project this plasmid belongs to', null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='plasmids', to='constructs.constructproject')),
                ('parent', models.ForeignKey(blank=True, help_text='Parent plasmid (if this is a derivative)', null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='children', to='constructs.plasmid')),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_plasmids', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Plasmid',
                'verbose_name_plural': 'Plasmids',
                'ordering': ['-ncn_id'],
            },
        ),
        # ProteinSynonym
        migrations.CreateModel(
            name='ProteinSynonym',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=256, unique=True)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='synonyms', to='constructs.protein')),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_protein_synonyms', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Protein Synonym',
                'verbose_name_plural': 'Protein Synonyms',
                'ordering': ['name'],
            },
        ),
        # ProteinUse
        migrations.CreateModel(
            name='ProteinUse',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='project_uses', to='constructs.protein')),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='protein_uses', to='constructs.constructproject')),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_protein_uses', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Protein Use',
                'verbose_name_plural': 'Protein Uses',
                'ordering': ['protein__uniprot_id'],
                'unique_together': {('protein', 'project')},
            },
        ),
        # Cassette
        migrations.CreateModel(
            name='Cassette',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('start', models.IntegerField(help_text='Start amino acid position')),
                ('end', models.IntegerField(help_text='End amino acid position')),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='cassettes', to='constructs.protein')),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_cassettes', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Cassette',
                'verbose_name_plural': 'Cassettes',
                'ordering': ['protein__uniprot_id', 'start'],
            },
        ),
        # CassetteUse
        migrations.CreateModel(
            name='CassetteUse',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('alignment_file', models.FileField(blank=True, help_text='Sequence alignment file', max_length=500, null=True, upload_to='constructs/alignments/')),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('cassette', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='plasmid_uses', to='constructs.cassette')),
                ('plasmid', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='cassette_uses', to='constructs.plasmid')),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_cassette_uses', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Cassette Use',
                'verbose_name_plural': 'Cassette Uses',
                'ordering': ['plasmid__ncn_id', 'cassette__protein__uniprot_id'],
            },
        ),
        # SequencingResult
        migrations.CreateModel(
            name='SequencingResult',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('file', models.FileField(help_text='Sequencing result file (.ab1, .seq, etc.)', max_length=500, upload_to='constructs/sequencing/')),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('cassette_use', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='sequencing_results', to='constructs.cassetteuse')),
                ('plasmid', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='sequencing_results', to='constructs.plasmid')),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_sequencing_results', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Sequencing Result',
                'verbose_name_plural': 'Sequencing Results',
                'ordering': ['-created_at'],
            },
        ),
        # ExpressionTag
        migrations.CreateModel(
            name='ExpressionTag',
            fields=[
                ('id', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('expression_tag_type', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='expression_tags', to='constructs.expressiontagtype')),
                ('protease', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='expression_tags', to='constructs.protease')),
                ('cassette_use', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='expression_tags', to='constructs.cassetteuse')),
                ('location', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='expression_tags', to='constructs.expressiontaglocation')),
                ('created_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='created_expression_tags', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Expression Tag',
                'verbose_name_plural': 'Expression Tags',
                'ordering': ['cassette_use__plasmid__ncn_id'],
            },
        ),
    ]
