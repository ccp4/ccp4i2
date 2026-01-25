"""
Compound Registry Models

Models for compound registration, batch tracking, and supplier management.
Migrated from legacy RegisterCompounds app with cleaned-up naming.

Model mapping from legacy:
    RegProjects -> Target
    RegSuppliers -> Supplier
    RegData -> Compound
    RegBatch -> Batch
    RegBatchQCFile -> BatchQCFile
    RegDataTemplate -> CompoundTemplate
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


class Supplier(models.Model):
    """
    Chemical supplier or synthesis source.

    Tracks where compounds/batches originate from (commercial vendors,
    internal synthesis, collaborators, etc.)

    When a user is linked, this supplier represents that user's personal
    synthesis source (e.g., "Martin Noble" as a supplier for compounds
    they synthesized themselves).
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=64, unique=True)
    initials = models.CharField(
        max_length=16,
        unique=True,
        blank=True,
        null=True,
        help_text="Short code used in barcodes, e.g., 'NCL', 'AST'"
    )
    user = models.OneToOneField(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='supplier',
        help_text="User account linked to this supplier (for personal suppliers)"
    )

    class Meta:
        ordering = ['name']
        verbose_name = 'Supplier'
        verbose_name_plural = 'Suppliers'

    def __str__(self):
        return self.name


def _target_image_path(instance, filename):
    """Generate upload path for target branding images.

    Path: compounds/registry/targets/{uuid}_{filename}
    """
    return f'compounds/registry/targets/{instance.id}_{filename}'


class Target(models.Model):
    """
    Drug discovery target or campaign.

    This was previously called 'RegProjects'. Renamed to avoid confusion
    with CCP4i2's crystallographic 'Project' concept.

    A Target represents a biological target (protein, pathway) that compounds
    are designed to modulate. Multiple compounds are registered against a target,
    and assays test compound activity against targets.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=255, unique=True)
    parent = models.ForeignKey(
        'self',
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='children',
        help_text="Parent target for hierarchical organization"
    )
    image = models.ImageField(
        upload_to=_target_image_path,
        max_length=255,
        blank=True,
        null=True,
        help_text="Branding image for the target dashboard"
    )
    saved_aggregation_view = models.JSONField(
        blank=True,
        null=True,
        help_text="Saved aggregation query configuration for dashboard display. "
                  "Schema: {protocol_names: string[], compound_search: string, "
                  "output_format: 'compact'|'medium'|'long', "
                  "aggregations: ('geomean'|'count'|'stdev'|'list')[], "
                  "status: 'valid'|'invalid'|'unassigned'|''}"
    )
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ['name']
        verbose_name = 'Target'
        verbose_name_plural = 'Targets'

    def __str__(self):
        return self.name


def _next_reg_number():
    """Generate next compound registration number."""
    max_dict = Compound.objects.aggregate(Max('reg_number'))
    highest = max_dict['reg_number__max']
    # Start at 26000 if no compounds exist (preserving legacy numbering)
    return (highest or 25999) + 1


def _compound_svg_path(instance, filename):
    """Generate upload path for compound SVG images.

    Path: compounds/registry/svg/NCL-XXXXXXXX.svg
    """
    return Path('compounds') / 'registry' / 'svg' / f'{instance.formatted_id}.svg'


class Compound(models.Model):
    """
    Registered compound with unique NCL-XXXXXXXX identifier.

    Each compound has a unique registration number and SMILES structure.
    Compounds belong to a Target and can have multiple Batches.
    """

    STEREO_CHOICES = [
        ('unset', 'Unset'),
        ('achiral', 'Achiral'),
        ('racemic', 'Racemic mixture'),
        ('single_unknown', 'Single enantiomer, configuration unknown'),
        ('r_enantiomer', 'R enantiomer'),
        ('s_enantiomer', 'S enantiomer'),
        ('non_racemic_mixture', 'Non-racemic stereoisomer mixture'),
        ('four_diastereomers', 'Mixture of 4 diastereoisomers'),
        ('two_diastereomers', 'Mixture of 2 diastereoisomers'),
        ('single_diastereomer_unknown', 'Single diastereoisomer, configuration unknown'),
        ('rr_diastereomer', 'RR diastereoisomer'),
        ('rs_diastereomer', 'RS diastereoisomer'),
        ('sr_diastereomer', 'SR diastereoisomer'),
        ('ss_diastereomer', 'SS diastereoisomer'),
        ('epimer_mixture', 'Mixture of epimers'),
        ('ez_mixture', 'Mixture of E and Z isomers'),
        ('e_isomer', 'E isomer'),
        ('z_isomer', 'Z isomer'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    reg_number = models.IntegerField(
        unique=True,
        editable=False,
        default=_next_reg_number,
        help_text="Auto-assigned registration number"
    )
    target = models.ForeignKey(
        Target,
        on_delete=models.PROTECT,
        related_name='compounds',
        help_text="Drug discovery target this compound was designed for"
    )

    # Chemistry
    smiles = models.TextField(help_text="SMILES structure notation")
    rdkit_smiles = models.TextField(
        blank=True,
        null=True,
        help_text="RDKit canonical SMILES (for duplicate detection)"
    )
    inchi = models.TextField(blank=True, null=True, help_text="InChI structure identifier")
    molecular_weight = models.FloatField(null=True, blank=True)
    stereo_comment = models.CharField(
        max_length=70,
        choices=STEREO_CHOICES,
        default='unset'
    )

    # Provenance
    supplier = models.ForeignKey(
        Supplier,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='compounds'
    )
    supplier_ref = models.CharField(
        max_length=255,
        blank=True,
        null=True,
        help_text="Supplier's catalog or reference number"
    )
    labbook_number = models.IntegerField(null=True, blank=True, help_text="Lab notebook number")
    page_number = models.IntegerField(null=True, blank=True, help_text="Page number in lab notebook")
    compound_number = models.IntegerField(
        default=1,
        help_text="Compound number on page (if multiple)"
    )

    # Audit
    registered_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='registered_compounds'
    )
    legacy_registered_by = models.CharField(
        max_length=100,
        blank=True,
        null=True,
        help_text="Original user_name from legacy system (for historical records)"
    )
    registered_at = models.DateTimeField(auto_now_add=True)
    modified_at = models.DateTimeField(auto_now=True, null=True)
    comments = models.TextField(blank=True, null=True)

    # Generated files
    svg_file = models.ImageField(
        upload_to=_compound_svg_path,
        max_length=255,
        blank=True,
        null=True,
        help_text="SVG structure image"
    )

    class Meta:
        ordering = ['-reg_number']
        verbose_name = 'Compound'
        verbose_name_plural = 'Compounds'

    def __str__(self):
        return self.formatted_id

    @cached_property
    def formatted_id(self):
        """NCL-XXXXXXXX format identifier."""
        return f'NCL-{self.reg_number:08d}'

    @cached_property
    def barcode(self):
        """Generate barcode string from supplier initials and notebook reference."""
        if self.supplier and self.supplier.initials:
            return f'{self.supplier.initials}-{self.labbook_number}-{self.page_number}'
        return None

    def save(self, *args, **kwargs):
        # Ensure reg_number is set
        if not self.reg_number:
            self.reg_number = _next_reg_number()

        # Check if SMILES has changed (compare with database value if this is an update)
        smiles_changed = False
        if self.pk:
            try:
                old_instance = Compound.objects.get(pk=self.pk)
                smiles_changed = old_instance.smiles != self.smiles
            except Compound.DoesNotExist:
                smiles_changed = True  # New record
        else:
            smiles_changed = True  # New record

        # Generate canonical SMILES and molecular weight via RDKit
        # Recalculate if SMILES is new/changed, or if rdkit_smiles is missing
        if self.smiles and (smiles_changed or not self.rdkit_smiles):
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors
                mol = Chem.MolFromSmiles(self.smiles)
                if mol:
                    self.rdkit_smiles = Chem.MolToSmiles(mol, canonical=True)
                    # Recalculate molecular weight if SMILES changed or not set
                    if smiles_changed or not self.molecular_weight:
                        self.molecular_weight = Descriptors.MolWt(mol)
            except ImportError:
                pass  # RDKit not available
            except Exception:
                pass  # Invalid SMILES or other RDKit error

        super().save(*args, **kwargs)


def _batch_qc_path(instance, filename):
    """Generate upload path for batch QC files.

    Path: compounds/registry/qc/NCL-XXXXXXXX/{uuid}_{filename}

    Files are organized under a dedicated qc subdirectory within
    compounds/registry, with per-compound directories named by their formatted ID.
    """
    compound_id = instance.batch.compound.formatted_id
    return f'compounds/registry/qc/{compound_id}/{instance.id}_{filename}'


class Batch(models.Model):
    """
    A synthesis batch of a compound.

    Each compound can have multiple batches representing different
    synthesis runs, potentially from different suppliers.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    compound = models.ForeignKey(
        Compound,
        on_delete=models.CASCADE,
        related_name='batches'
    )
    batch_number = models.IntegerField(
        help_text="Auto-assigned within compound"
    )

    # Provenance
    supplier = models.ForeignKey(
        Supplier,
        on_delete=models.SET_NULL,
        null=True,
        blank=True
    )
    supplier_ref = models.CharField(max_length=255, blank=True, null=True)
    labbook_number = models.IntegerField(null=True, blank=True)
    page_number = models.IntegerField(null=True, blank=True)

    # Properties
    amount = models.DecimalField(
        max_digits=10,
        decimal_places=4,
        null=True,
        blank=True,
        help_text="Amount in mg"
    )
    salt_code = models.CharField(max_length=128, blank=True, null=True)
    molecular_weight = models.FloatField(
        null=True,
        blank=True,
        help_text="MW including salt (if different from parent)"
    )

    # Metadata
    registered_at = models.DateTimeField(auto_now_add=True)
    comments = models.TextField(blank=True, null=True)

    class Meta:
        ordering = ['compound', 'batch_number']
        # Note: unique_together removed because legacy data has some duplicates
        verbose_name = 'Batch'
        verbose_name_plural = 'Batches'

    def __str__(self):
        return f'{self.compound.formatted_id}/{self.batch_number}'

    def save(self, *args, **kwargs):
        # Auto-assign batch number if not set
        if self.batch_number is None or self.batch_number < 0:
            max_dict = Batch.objects.filter(
                compound=self.compound
            ).aggregate(Max('batch_number'))
            highest = max_dict['batch_number__max']
            self.batch_number = (highest or 0) + 1
        super().save(*args, **kwargs)


class BatchQCFile(models.Model):
    """
    QC documentation for a batch (NMR, LCMS, HPLC, etc.).

    Each batch can have multiple QC files as evidence of purity/identity.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    batch = models.ForeignKey(
        Batch,
        on_delete=models.CASCADE,
        related_name='qc_files'
    )
    file = models.FileField(upload_to=_batch_qc_path, max_length=255)
    original_filename = models.CharField(
        max_length=255,
        blank=True,
        null=True,
        help_text="Original filename as uploaded by user"
    )
    comments = models.TextField(blank=True, null=True)
    uploaded_at = models.DateTimeField(auto_now_add=True, null=True, blank=True)

    class Meta:
        verbose_name = 'Batch QC File'
        verbose_name_plural = 'Batch QC Files'

    def __str__(self):
        return f'QC for {self.batch}'

    def save(self, *args, **kwargs):
        # Capture original filename on first save
        if self.file and not self.original_filename:
            self.original_filename = Path(self.file.name).name
        super().save(*args, **kwargs)

    @property
    def filename(self):
        """Return original filename if available, otherwise extract from path."""
        if self.original_filename:
            return self.original_filename
        return Path(self.file.name).name if self.file else None


@receiver(post_delete, sender=BatchQCFile)
def _batch_qc_file_post_delete(sender, instance, **kwargs):
    """Delete associated file when BatchQCFile is deleted."""
    delete_file_field(instance.file)


def _template_svg_path(instance, filename):
    """Generate upload path for template SVG images.

    Path: compounds/registry/templates/{uuid}.svg
    """
    return Path('compounds') / 'registry' / 'templates' / f'{instance.id}.svg'


class CompoundTemplate(models.Model):
    """
    Scaffold/template structure for a target.

    Templates define core scaffolds that compounds in a project are based on.
    Useful for SAR analysis and compound design.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    target = models.ForeignKey(
        Target,
        on_delete=models.CASCADE,
        related_name='templates',
        null=True,
        blank=True
    )
    name = models.CharField(max_length=512)
    mol2d = models.TextField(help_text="2D MOL block for the template structure")
    svg_file = models.ImageField(
        upload_to=_template_svg_path,
        max_length=255,
        blank=True,
        null=True
    )

    class Meta:
        ordering = ['name']
        verbose_name = 'Compound Template'
        verbose_name_plural = 'Compound Templates'

    def __str__(self):
        return self.name
