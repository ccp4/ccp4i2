"""
Construct Database Models

Models for plasmid constructs, proteins, cassettes, and sequencing results.
Migrated from legacy ConstructDatabase app with modernized naming.

Model mapping from legacy:
    Project -> ConstructProject (renamed to avoid confusion with CCP4i2 Project)
    Plasmid -> Plasmid
    Protein -> Protein
    ProteinUse -> ProteinUse
    ProteinSynonym -> ProteinSynonym
    Cassette -> Cassette
    CassetteUse -> CassetteUse
    SequencingResult -> SequencingResult
    ExpressionTagType -> ExpressionTagType
    Protease -> Protease
    ExpressionTagLocation -> ExpressionTagLocation
    ExpressionTag -> ExpressionTag
"""

import uuid
from pathlib import Path

from django.contrib.auth import get_user_model
from django.db import models
from django.db.models import Max
from django.db.models.signals import post_delete
from django.dispatch import receiver
from django.utils.functional import cached_property

from compounds.utils import delete_file_field

User = get_user_model()


# =============================================================================
# File Upload Paths
# =============================================================================

def _plasmid_file_path(instance, filename):
    """Generate upload path for plasmid files (GenBank, etc.).

    Legacy pattern: ConstructDatabase/NCLCON-XXXXXXXX/{filename}
    """
    return Path('ConstructDatabase') / instance.formatted_id / filename


def _sequencing_file_path(instance, filename):
    """Generate upload path for sequencing result files.

    Legacy pattern: ConstructDatabase/NCLCON-XXXXXXXX/{filename}
    (Same directory as plasmid files)
    """
    return Path('ConstructDatabase') / instance.plasmid.formatted_id / filename


def _alignment_file_path(instance, filename):
    """Generate upload path for alignment files.

    Legacy had no custom path (just filename in MEDIA_ROOT).
    Keeping ConstructDatabase structure for consistency.
    """
    return Path('ConstructDatabase') / instance.plasmid.formatted_id / filename


# =============================================================================
# Helper Functions
# =============================================================================

def _next_ncn_id():
    """Generate next plasmid NCN ID number."""
    max_dict = Plasmid.objects.aggregate(Max('ncn_id'))
    highest = max_dict['ncn_id__max']
    return (highest or 0) + 1


# =============================================================================
# Project Model
# =============================================================================

class ConstructProject(models.Model):
    """
    Hierarchical project for organizing plasmid constructs.

    Renamed from 'Project' to avoid confusion with CCP4i2's crystallographic
    Project concept.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=256)
    parent = models.ForeignKey(
        'self',
        on_delete=models.CASCADE,
        null=True,
        blank=True,
        related_name='children',
        help_text="Parent project for hierarchical organization"
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_construct_projects'
    )

    class Meta:
        ordering = ['name']
        verbose_name = 'Construct Project'
        verbose_name_plural = 'Construct Projects'

    def __str__(self):
        return self.name


# =============================================================================
# Plasmid Model
# =============================================================================

class Plasmid(models.Model):
    """
    Plasmid construct with unique NCLCON-XXXXXXXX identifier.

    Stores plasmid design files (GenBank format) and links to cassettes
    and sequencing results.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    ncn_id = models.IntegerField(
        unique=True,
        editable=False,
        default=_next_ncn_id,
        help_text="Auto-assigned construct number"
    )
    name = models.CharField(
        max_length=256,
        unique=True,
        help_text="Descriptive name for the plasmid"
    )
    project = models.ForeignKey(
        ConstructProject,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='plasmids',
        help_text="Project this plasmid belongs to"
    )
    parent = models.ForeignKey(
        'self',
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='children',
        help_text="Parent plasmid (if this is a derivative)"
    )
    genbank_file = models.FileField(
        upload_to=_plasmid_file_path,
        max_length=500,
        blank=True,
        null=True,
        help_text="GenBank or SnapGene file (.gb, .gbk)"
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_plasmids'
    )

    class Meta:
        ordering = ['-ncn_id']
        verbose_name = 'Plasmid'
        verbose_name_plural = 'Plasmids'

    def __str__(self):
        return self.formatted_id

    @cached_property
    def formatted_id(self):
        """NCLCON-XXXXXXXX format identifier."""
        return f'NCLCON-{self.ncn_id:08d}'

    def save(self, *args, **kwargs):
        # Ensure ncn_id is set
        if not self.ncn_id:
            self.ncn_id = _next_ncn_id()

        # Auto-generate name if not provided
        if not self.name:
            self.name = f'NCLCON-{self.ncn_id:08d}'

        super().save(*args, **kwargs)


@receiver(post_delete, sender=Plasmid)
def _plasmid_post_delete(sender, instance, **kwargs):
    """Delete associated file when plasmid is deleted."""
    delete_file_field(instance.genbank_file)


# =============================================================================
# Protein Models
# =============================================================================

class Protein(models.Model):
    """
    Protein definition identified by UniProt ID.

    Proteins can have multiple synonyms and are used to define cassettes.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    uniprot_id = models.CharField(
        max_length=256,
        unique=True,
        help_text="UniProt identifier (e.g., P04637)"
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_proteins'
    )

    class Meta:
        ordering = ['uniprot_id']
        verbose_name = 'Protein'
        verbose_name_plural = 'Proteins'

    def __str__(self):
        return self.uniprot_id


class ProteinSynonym(models.Model):
    """
    Alternative name for a protein.

    Allows proteins to be referenced by common names, gene names, etc.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=256, unique=True)
    protein = models.ForeignKey(
        Protein,
        on_delete=models.CASCADE,
        related_name='synonyms'
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_protein_synonyms'
    )

    class Meta:
        ordering = ['name']
        verbose_name = 'Protein Synonym'
        verbose_name_plural = 'Protein Synonyms'

    def __str__(self):
        return f'{self.name} ({self.protein.uniprot_id})'


class ProteinUse(models.Model):
    """
    Links a protein to a project.

    Tracks which proteins are being studied in which projects.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    protein = models.ForeignKey(
        Protein,
        on_delete=models.CASCADE,
        related_name='project_uses'
    )
    project = models.ForeignKey(
        ConstructProject,
        on_delete=models.CASCADE,
        related_name='protein_uses'
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_protein_uses'
    )

    class Meta:
        ordering = ['protein__uniprot_id']
        unique_together = ['protein', 'project']
        verbose_name = 'Protein Use'
        verbose_name_plural = 'Protein Uses'

    def __str__(self):
        return f'{self.protein.uniprot_id} in {self.project.name}'


# =============================================================================
# Cassette Models
# =============================================================================

class Cassette(models.Model):
    """
    A protein region/domain defined by start and end amino acid positions.

    Cassettes represent fragments of proteins that are cloned into plasmids.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    protein = models.ForeignKey(
        Protein,
        on_delete=models.CASCADE,
        related_name='cassettes'
    )
    start = models.IntegerField(help_text="Start amino acid position")
    end = models.IntegerField(help_text="End amino acid position")
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_cassettes'
    )

    class Meta:
        ordering = ['protein__uniprot_id', 'start']
        verbose_name = 'Cassette'
        verbose_name_plural = 'Cassettes'

    def __str__(self):
        return f'{self.protein.uniprot_id}({self.start}-{self.end})'

    @property
    def display_name(self):
        """Human-readable cassette name."""
        return f'{self.protein.uniprot_id}({self.start}-{self.end})'


class CassetteUse(models.Model):
    """
    Links a cassette to a plasmid.

    Represents the use of a specific protein cassette in a plasmid construct,
    optionally with an alignment file showing how it was cloned.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    cassette = models.ForeignKey(
        Cassette,
        on_delete=models.CASCADE,
        related_name='plasmid_uses'
    )
    plasmid = models.ForeignKey(
        Plasmid,
        on_delete=models.CASCADE,
        related_name='cassette_uses'
    )
    alignment_file = models.FileField(
        upload_to=_alignment_file_path,
        max_length=500,
        blank=True,
        null=True,
        help_text="Sequence alignment file"
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_cassette_uses'
    )

    class Meta:
        ordering = ['plasmid__ncn_id', 'cassette__protein__uniprot_id']
        verbose_name = 'Cassette Use'
        verbose_name_plural = 'Cassette Uses'

    def __str__(self):
        return f'{self.cassette} in {self.plasmid.formatted_id}'


@receiver(post_delete, sender=CassetteUse)
def _cassette_use_post_delete(sender, instance, **kwargs):
    """Delete associated alignment file when cassette use is deleted."""
    delete_file_field(instance.alignment_file)


# =============================================================================
# Sequencing Result Model
# =============================================================================

class SequencingResult(models.Model):
    """
    Sequencing validation result for a cassette in a plasmid.

    Stores sequencing files that validate the construct was made correctly.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    cassette_use = models.ForeignKey(
        CassetteUse,
        on_delete=models.CASCADE,
        related_name='sequencing_results'
    )
    plasmid = models.ForeignKey(
        Plasmid,
        on_delete=models.CASCADE,
        related_name='sequencing_results'
    )
    file = models.FileField(
        upload_to=_sequencing_file_path,
        max_length=500,
        help_text="Sequencing result file (.ab1, .seq, etc.)"
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_sequencing_results'
    )

    class Meta:
        ordering = ['-created_at']
        verbose_name = 'Sequencing Result'
        verbose_name_plural = 'Sequencing Results'

    def __str__(self):
        return f'Sequencing for {self.plasmid.formatted_id}'

    @property
    def filename(self):
        """Extract filename from path."""
        return Path(self.file.name).name if self.file else None


@receiver(post_delete, sender=SequencingResult)
def _sequencing_result_post_delete(sender, instance, **kwargs):
    """Delete associated file when sequencing result is deleted."""
    delete_file_field(instance.file)


# =============================================================================
# Reference Data Models (Expression Tags)
# =============================================================================

class ExpressionTagType(models.Model):
    """
    Type of expression tag (e.g., GST, SUMO, His6, etc.).

    Pre-populated reference data.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=256, unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_expression_tag_types'
    )

    class Meta:
        ordering = ['name']
        verbose_name = 'Expression Tag Type'
        verbose_name_plural = 'Expression Tag Types'

    def __str__(self):
        return self.name


class Protease(models.Model):
    """
    Protease used to cleave expression tags (e.g., TEV, 3C, Thrombin).

    Pre-populated reference data.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=256, unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_proteases'
    )

    class Meta:
        ordering = ['name']
        verbose_name = 'Protease'
        verbose_name_plural = 'Proteases'

    def __str__(self):
        return self.name


class ExpressionTagLocation(models.Model):
    """
    Location of expression tag relative to the cassette (N-term, C-term, Internal).

    Pre-populated reference data.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=256, unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_expression_tag_locations'
    )

    class Meta:
        ordering = ['name']
        verbose_name = 'Expression Tag Location'
        verbose_name_plural = 'Expression Tag Locations'

    def __str__(self):
        return self.name


class ExpressionTag(models.Model):
    """
    Expression tag on a cassette use.

    Combines tag type, location, and optional protease cleavage site.
    Note: Fixed FK bug from legacy - cassetteUse now correctly references CassetteUse.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    expression_tag_type = models.ForeignKey(
        ExpressionTagType,
        on_delete=models.CASCADE,
        related_name='expression_tags'
    )
    protease = models.ForeignKey(
        Protease,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='expression_tags'
    )
    cassette_use = models.ForeignKey(
        CassetteUse,  # Fixed: was incorrectly ExpressionTagType in legacy
        on_delete=models.CASCADE,
        related_name='expression_tags'
    )
    location = models.ForeignKey(
        ExpressionTagLocation,
        on_delete=models.CASCADE,
        related_name='expression_tags'
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_expression_tags'
    )

    class Meta:
        ordering = ['cassette_use__plasmid__ncn_id']
        verbose_name = 'Expression Tag'
        verbose_name_plural = 'Expression Tags'

    def __str__(self):
        parts = [self.expression_tag_type.name, self.location.name]
        if self.protease:
            parts.append(f'({self.protease.name} cleavable)')
        return ' '.join(parts)
