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

from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import models
from django.db.models import Max
from django.db.models.signals import post_delete
from django.dispatch import receiver
from django.utils.functional import cached_property

from compounds.formatting import format_compound_id
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

    # Audit
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_suppliers'
    )
    modified_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='modified_suppliers'
    )
    created_at = models.DateTimeField(auto_now_add=True, null=True)
    modified_at = models.DateTimeField(auto_now=True, null=True)

    class Meta:
        ordering = ['name']
        verbose_name = 'Supplier'
        verbose_name_plural = 'Suppliers'

    def __str__(self):
        return self.name


class LabNotebookEntry(models.Model):
    """
    Represents one ELN page or one paper notebook page.

    Many compounds may reference the same entry. This model supports both:
    - ELN-style references (URL to OneNote, LabArchives, etc.)
    - Paper notebook references (labbook_number + page_number)

    The label property generates identifiers like 'KF001', 'KF022' from
    supplier initials + sequence_number.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    supplier = models.ForeignKey(
        Supplier,
        on_delete=models.PROTECT,
        related_name='notebook_entries',
        help_text="Supplier/chemist who owns this notebook entry"
    )
    sequence_number = models.IntegerField(
        help_text="Per-supplier sequential index (e.g., 1 for KF001, 22 for KF022)"
    )
    title = models.CharField(
        max_length=255,
        blank=True,
        help_text='Page title, e.g., "SA157 Cyclic peptide synthesis and characterisation"'
    )
    date = models.DateField(
        null=True,
        blank=True,
        help_text="Date of the notebook entry"
    )

    # ELN form
    url = models.URLField(
        max_length=2048,
        blank=True,
        help_text="Hyperlink to ELN page (OneNote, LabArchives, etc.)"
    )

    # Paper-notebook form (mirrors legacy Compound.labbook_number/page_number)
    labbook_number = models.IntegerField(
        null=True,
        blank=True,
        help_text="Lab notebook number (for paper notebooks)"
    )
    page_number = models.IntegerField(
        null=True,
        blank=True,
        help_text="Page number in lab notebook (for paper notebooks)"
    )

    # Audit
    created_at = models.DateTimeField(auto_now_add=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_notebook_entries'
    )

    class Meta:
        unique_together = [('supplier', 'sequence_number')]
        ordering = ['supplier', 'sequence_number']
        verbose_name = 'Lab Notebook Entry'
        verbose_name_plural = 'Lab Notebook Entries'

    def __str__(self):
        return self.label

    @cached_property
    def label(self):
        """
        Generate label like 'KF001', 'KF022'.

        Uses supplier initials + zero-padded sequence number (3 digits).
        Falls back to '#001' format if supplier has no initials.
        """
        if self.supplier and self.supplier.initials:
            return f'{self.supplier.initials}{self.sequence_number:03d}'
        return f'#{self.sequence_number:03d}'

    def clean(self):
        """Validate that ELN url and paper labbook fields are not both set."""
        from django.core.exceptions import ValidationError
        if self.url and (self.labbook_number or self.page_number):
            raise ValidationError(
                "Set either ELN url or paper labbook_number+page_number, not both"
            )


def _target_image_path(instance, filename):
    """Generate upload path for target branding images.

    Path: compounds/registry/targets/{uuid}_{filename}
    """
    return f'compounds/registry/targets/{instance.id}_{filename}'


class Gene(models.Model):
    """
    HGNC-hydrated metadata for a single gene.

    Canonical by symbol. HGNC is the source of truth for human gene
    nomenclature; this table caches it so downstream matching (e.g. NLP
    query target resolution) does not depend on runtime external calls.

    See apps/compounds/docs/TARGET_MODEL_PROPOSAL.md for the rationale.
    """

    HYDRATION_STATUS_CHOICES = [
        ('pending', 'Pending'),
        ('ok', 'OK'),
        ('partial', 'Partial'),
        ('failed', 'Failed'),
        ('manual', 'Manual override'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    symbol = models.CharField(
        max_length=32,
        unique=True,
        help_text="HGNC-approved gene symbol, e.g. 'EGFR', 'SKP2'."
    )
    hgnc_id = models.CharField(
        max_length=32,
        blank=True,
        db_index=True,
        help_text="HGNC stable ID, e.g. 'HGNC:3236'."
    )
    name = models.CharField(
        max_length=255,
        blank=True,
        help_text="Full gene name from HGNC, e.g. 'epidermal growth factor receptor'."
    )
    aliases = models.JSONField(
        default=list,
        blank=True,
        help_text="Union of HGNC alias_symbol and prev_symbol lists."
    )
    uniprot_ids = models.JSONField(
        default=list,
        blank=True,
        help_text="UniProt accession(s); typically one primary + isoforms."
    )
    ensembl_gene_id = models.CharField(
        max_length=32,
        blank=True,
        db_index=True,
        help_text="Ensembl gene ID, e.g. 'ENSG00000146648'."
    )

    hydration_status = models.CharField(
        max_length=16,
        choices=HYDRATION_STATUS_CHOICES,
        default='pending',
        help_text="State of HGNC hydration for this gene."
    )
    hydrated_at = models.DateTimeField(
        null=True,
        blank=True,
        help_text="Last successful hydration timestamp."
    )
    hydration_source = models.CharField(
        max_length=32,
        blank=True,
        help_text="'hgnc-api', 'hgnc-snapshot', 'manual', etc."
    )

    created_at = models.DateTimeField(auto_now_add=True)
    modified_at = models.DateTimeField(auto_now=True)

    class Meta:
        ordering = ['symbol']
        verbose_name = 'Gene'
        verbose_name_plural = 'Genes'

    def __str__(self):
        return self.symbol

    def save(self, *args, **kwargs):
        # Normalise symbol to uppercase — HGNC convention.
        if self.symbol:
            self.symbol = self.symbol.strip().upper()
        super().save(*args, **kwargs)


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
    genes = models.ManyToManyField(
        Gene,
        related_name='targets',
        blank=True,
        help_text="Biological gene(s) underlying this target. "
                  "Multi-gene for PPI / complex / pathway targets."
    )
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
                  "output_format: 'compact'|'medium'|'long'|'pivot'|'cards', "
                  "aggregations: ('geomean'|'count'|'stdev'|'list')[], "
                  "status: 'valid'|'invalid'|'unassigned'|'', "
                  "concentration_display?: 'natural'|'nM'|'uM'|'mM'|'pConc', "
                  "include_properties?: string[]}"
    )
    scorecard_config = models.JSONField(
        blank=True,
        null=True,
        help_text="Per-target scorecard definition for spider / radar visualisation. "
                  "Schema: {axes: ScorecardAxis[]} where each axis is "
                  "{label: string, kind: 'protocol'|'ratio'|'worst_of'|'lipinski', "
                  "target_value?: number, poor_value?: number, "
                  "threshold_scale?: 'log'|'linear', "
                  "protocol_id?: uuid (kind='protocol'), "
                  "numerator_id?: uuid, denominator_id?: uuid (kind='ratio'), "
                  "protocol_ids?: uuid[] (kind='worst_of')}. "
                  "kind='lipinski' takes no extra fields; score is 0-4."
    )

    # Audit
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_targets'
    )
    modified_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='modified_targets'
    )
    created_at = models.DateTimeField(auto_now_add=True)
    modified_at = models.DateTimeField(auto_now=True, null=True)

    class Meta:
        ordering = ['name']
        verbose_name = 'Target'
        verbose_name_plural = 'Targets'

    def __str__(self):
        return self.name

    def clean(self):
        super().clean()
        if self.scorecard_config is not None:
            _validate_scorecard_config(self.scorecard_config)


_VALID_AXIS_KINDS = {'protocol', 'ratio', 'worst_of', 'lipinski'}
_VALID_THRESHOLD_SCALES = {'log', 'linear'}


def _validate_scorecard_config(config):
    """Structural validator for Target.scorecard_config.

    Checks only the shape — doesn't verify protocol IDs exist or that
    ratio-axis kpi_units match. Semantic validation can live in the
    serializer or be deferred to client-side once the editor is in place.
    """
    from django.core.exceptions import ValidationError

    if not isinstance(config, dict):
        raise ValidationError({'scorecard_config': 'Must be an object.'})
    axes = config.get('axes')
    if axes is None:
        return  # Empty config is fine
    if not isinstance(axes, list):
        raise ValidationError({'scorecard_config': '"axes" must be a list.'})

    for i, axis in enumerate(axes):
        prefix = f'axes[{i}]'
        if not isinstance(axis, dict):
            raise ValidationError({'scorecard_config': f'{prefix}: must be an object.'})
        kind = axis.get('kind')
        if kind not in _VALID_AXIS_KINDS:
            raise ValidationError({'scorecard_config': f'{prefix}.kind: must be one of {sorted(_VALID_AXIS_KINDS)}.'})
        if not isinstance(axis.get('label'), str) or not axis['label'].strip():
            raise ValidationError({'scorecard_config': f'{prefix}.label: required non-empty string.'})

        if kind == 'protocol':
            if not axis.get('protocol_id'):
                raise ValidationError({'scorecard_config': f'{prefix}.protocol_id: required.'})
        elif kind == 'ratio':
            if not axis.get('numerator_id') or not axis.get('denominator_id'):
                raise ValidationError({'scorecard_config': f'{prefix}: ratio requires numerator_id and denominator_id.'})
            if axis['numerator_id'] == axis['denominator_id']:
                raise ValidationError({'scorecard_config': f'{prefix}: numerator and denominator must differ.'})
        elif kind == 'worst_of':
            ids = axis.get('protocol_ids')
            if not isinstance(ids, list) or len(ids) < 1:
                raise ValidationError({'scorecard_config': f'{prefix}: worst_of requires non-empty protocol_ids list.'})
        # kind == 'lipinski' has no kind-specific fields.

        # Threshold anchors are optional, but if both set they must differ
        # and the scale, if present, must be one of the allowed values.
        target_value = axis.get('target_value')
        poor_value = axis.get('poor_value')
        if target_value is not None and poor_value is not None and target_value == poor_value:
            raise ValidationError({'scorecard_config': f'{prefix}: target_value and poor_value must differ.'})
        scale = axis.get('threshold_scale')
        if scale is not None and scale not in _VALID_THRESHOLD_SCALES:
            raise ValidationError({'scorecard_config': f'{prefix}.threshold_scale: must be "log" or "linear".'})

        # Optional sector tag for visual grouping. Free string; we cap the
        # length and reject obvious junk but don't enforce a vocabulary so
        # projects can invent their own categories.
        sector = axis.get('sector')
        if sector is not None:
            if not isinstance(sector, str):
                raise ValidationError({'scorecard_config': f'{prefix}.sector: must be a string.'})
            if len(sector) > 64:
                raise ValidationError({'scorecard_config': f'{prefix}.sector: must be ≤ 64 characters.'})


def _next_reg_number():
    """Generate next compound registration number.

    Uses settings.COMPOUND_ID_START (default 1) as the first number when the
    database has no compounds yet. Existing deployments that want to preserve
    legacy numbering can set COMPOUND_ID_START=26000 via environment variable.
    """
    highest = Compound.objects.aggregate(Max('reg_number'))['reg_number__max']
    if highest is None:
        return getattr(settings, 'COMPOUND_ID_START', 1)
    return highest + 1


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
    ] + [
        (f'isomer_{i}', f'Isomer {i}') for i in range(1, 21)
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

    # ELN-linked provenance (Phase 1)
    notebook_entry = models.ForeignKey(
        'LabNotebookEntry',
        on_delete=models.PROTECT,
        null=True,
        blank=True,
        related_name='compounds',
        help_text="Lab notebook entry (ELN page or paper notebook page)"
    )
    helm_notation = models.CharField(
        max_length=4096,
        blank=True,
        help_text="Optional HELM string. Phase 1: opaque text, no validation."
    )
    sequence_display = models.CharField(
        max_length=512,
        blank=True,
        help_text="Human-readable sequence with modifications shown inline, e.g., "
                  "'FAM-linker-MCDWDIYRFPNHHC(1,4-Xylene)-NH2'. Free text."
    )

    # Audit
    registered_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='registered_compounds'
    )
    modified_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='modified_compounds',
        help_text="User who last modified this compound"
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

    # Alternative identifiers
    aliases = models.JSONField(
        default=list,
        blank=True,
        help_text="Alternative names/identifiers for this compound (e.g., supplier codes, "
                  "abbreviations, internal project names). Used as fallback during import matching."
    )

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
        """Formatted compound identifier (e.g., NCL-00026042)."""
        return format_compound_id(self.reg_number)

    @cached_property
    def barcode(self):
        """Generate barcode string from notebook entry or supplier initials.

        Returns:
            - notebook_entry.label (e.g., 'KF001') if notebook_entry is set
            - 'INITIALS-labbook-page' format for legacy paper notebook references
            - None if no provenance information is available
        """
        if self.notebook_entry:
            return self.notebook_entry.label
        if self.supplier and self.supplier.initials and self.labbook_number is not None:
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

        # Generate canonical SMILES, InChI, and molecular weight via RDKit
        # Recalculate if SMILES is new/changed, or if rdkit_smiles is missing
        if self.smiles and (smiles_changed or not self.rdkit_smiles):
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors
                mol = Chem.MolFromSmiles(self.smiles)
                if mol:
                    # Canonical SMILES
                    if smiles_changed or not self.rdkit_smiles:
                        self.rdkit_smiles = Chem.MolToSmiles(mol, canonical=True)
                    # Molecular weight
                    if smiles_changed or not self.molecular_weight:
                        self.molecular_weight = Descriptors.MolWt(mol)
                    # InChI - only if support is available
                    if smiles_changed or not self.inchi:
                        try:
                            from rdkit.Chem import inchi as _inchi_mod
                            if getattr(_inchi_mod, 'INCHI_AVAILABLE', False):
                                self.inchi = _inchi_mod.MolToInchi(mol)
                        except (ImportError, AttributeError):
                            pass
            except ImportError:
                pass  # RDKit not available
            except Exception:
                pass  # Invalid SMILES or other RDKit error

        super().save(*args, **kwargs)

        # Calculate molecular properties after compound is saved
        # (requires compound to have a PK for the OneToOne relationship)
        if smiles_changed and self.smiles:
            MolecularProperties.calculate_for_compound(self)


def _compound_document_path(instance, filename):
    """Generate upload path for compound document files.

    Path: compounds/registry/documents/NCL-XXXXXXXX/{uuid}_{filename}
    """
    compound_id = instance.compound.formatted_id
    return f'compounds/registry/documents/{compound_id}/{instance.id}_{filename}'


class CompoundDocument(models.Model):
    """
    Per-compound document attachment (chemdraw, QC report, spectrum, etc.).

    Documents can be either:
    - External URLs (e.g., links to SharePoint, cloud storage)
    - Uploaded files

    Exactly one of `url` or `file` should be set.

    Note: This is distinct from BatchQCFile which is per-batch. CompoundDocument
    is for compound-level documents like Chemdraw structure files that apply
    to all batches of a compound.
    """

    KIND_CHOICES = [
        ('chemdraw', 'Chemdraw structure'),
        ('qc', 'QC report'),
        ('spectrum', 'Spectrum / chromatogram'),
        ('other', 'Other'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    compound = models.ForeignKey(
        Compound,
        on_delete=models.CASCADE,
        related_name='documents'
    )
    kind = models.CharField(
        max_length=16,
        choices=KIND_CHOICES,
        help_text="Type of document"
    )
    label = models.CharField(
        max_length=255,
        blank=True,
        help_text="Display text, e.g., 'CP1.cdxml'"
    )

    # Either url or file should be set, not both
    url = models.URLField(
        max_length=2048,
        blank=True,
        help_text="External URL to document (SharePoint, cloud storage, etc.)"
    )
    file = models.FileField(
        upload_to=_compound_document_path,
        max_length=255,
        blank=True,
        help_text="Uploaded document file"
    )

    # Audit
    created_at = models.DateTimeField(auto_now_add=True)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_compound_documents'
    )

    class Meta:
        ordering = ['compound', 'kind', 'created_at']
        verbose_name = 'Compound Document'
        verbose_name_plural = 'Compound Documents'

    def __str__(self):
        return f'{self.compound.formatted_id} - {self.get_kind_display()}: {self.label or self.url or self.file.name}'

    def clean(self):
        """Validate that exactly one of url or file is set."""
        from django.core.exceptions import ValidationError
        has_url = bool(self.url)
        has_file = bool(self.file)
        if has_url == has_file:
            raise ValidationError("Set exactly one of url or file, not both or neither")


@receiver(post_delete, sender=CompoundDocument)
def _compound_document_post_delete(sender, instance, **kwargs):
    """Delete associated file when CompoundDocument is deleted."""
    if instance.file:
        delete_file_field(instance.file)


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

    # Audit
    registered_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='registered_batches',
        help_text="User who registered this batch"
    )
    modified_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='modified_batches',
        help_text="User who last modified this batch"
    )
    registered_at = models.DateTimeField(auto_now_add=True)
    modified_at = models.DateTimeField(auto_now=True, null=True)
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
    Documents can be either:
    - External URLs (e.g., links to SharePoint, ELN, cloud storage)
    - Uploaded files

    Exactly one of `url` or `file` should be set.
    """

    KIND_CHOICES = [
        ('lcms', 'LC-MS'),
        ('hrms', 'HR-MS'),
        ('nmr', 'NMR'),
        ('hplc', 'HPLC'),
        ('other', 'Other'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    batch = models.ForeignKey(
        Batch,
        on_delete=models.CASCADE,
        related_name='qc_files'
    )
    kind = models.CharField(
        max_length=16,
        choices=KIND_CHOICES,
        default='other',
        help_text="Type of QC document"
    )
    label = models.CharField(
        max_length=255,
        blank=True,
        help_text="Display text, e.g., 'LCMS-001'"
    )

    # Either url or file should be set, not both
    url = models.URLField(
        max_length=2048,
        blank=True,
        help_text="External URL to QC document (SharePoint, ELN, cloud storage)"
    )
    file = models.FileField(
        upload_to=_batch_qc_path,
        max_length=255,
        blank=True,
        help_text="Uploaded QC file"
    )
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


class MolecularProperties(models.Model):
    """
    Cached computed molecular properties for a compound.

    Stores drug-likeness descriptors calculated from SMILES via RDKit.
    Separated from Compound model to keep registration facts distinct
    from derived/computed values.

    Properties are recalculated when compound SMILES changes.
    """

    compound = models.OneToOneField(
        Compound,
        on_delete=models.CASCADE,
        primary_key=True,
        related_name='molecular_properties'
    )

    # Lipinski Rule of 5 properties
    molecular_weight = models.FloatField(
        null=True,
        blank=True,
        help_text="Exact molecular weight in Daltons (Lipinski: ≤500)"
    )
    heavy_atom_count = models.IntegerField(
        null=True,
        blank=True,
        help_text="Number of non-hydrogen atoms (for ligand efficiency)"
    )
    hbd = models.IntegerField(
        null=True,
        blank=True,
        help_text="Hydrogen bond donors (Lipinski: ≤5)"
    )
    hba = models.IntegerField(
        null=True,
        blank=True,
        help_text="Hydrogen bond acceptors (Lipinski: ≤10)"
    )
    clogp = models.FloatField(
        null=True,
        blank=True,
        help_text="Calculated LogP - Wildman-Crippen method (Lipinski: ≤5)"
    )
    tpsa = models.FloatField(
        null=True,
        blank=True,
        help_text="Topological polar surface area in Å² (oral bioavailability: ≤140)"
    )
    rotatable_bonds = models.IntegerField(
        null=True,
        blank=True,
        help_text="Number of rotatable bonds (oral bioavailability: ≤10)"
    )
    fraction_sp3 = models.FloatField(
        null=True,
        blank=True,
        help_text="Fraction of sp3 carbons (Fsp3) - 3D complexity indicator"
    )

    # Metadata
    calculated_at = models.DateTimeField(auto_now=True)
    rdkit_version = models.CharField(
        max_length=32,
        blank=True,
        null=True,
        help_text="RDKit version used for calculation"
    )

    class Meta:
        verbose_name = 'Molecular Properties'
        verbose_name_plural = 'Molecular Properties'

    def __str__(self):
        return f'Properties for {self.compound.formatted_id}'

    @classmethod
    def calculate_for_compound(cls, compound: Compound) -> 'MolecularProperties':
        """
        Calculate molecular properties from compound SMILES.

        Creates or updates MolecularProperties for the given compound.
        Returns the MolecularProperties instance (saved).
        """
        smiles = compound.smiles or compound.rdkit_smiles
        if not smiles:
            return None

        try:
            from rdkit import Chem, rdBase
            from rdkit.Chem import Descriptors, Lipinski

            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None

            props, _ = cls.objects.get_or_create(compound=compound)
            props.molecular_weight = Descriptors.ExactMolWt(mol)
            props.heavy_atom_count = Lipinski.HeavyAtomCount(mol)
            props.hbd = Descriptors.NumHDonors(mol)
            props.hba = Descriptors.NumHAcceptors(mol)
            props.clogp = Descriptors.MolLogP(mol)
            props.tpsa = Descriptors.TPSA(mol)
            props.rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            props.fraction_sp3 = Descriptors.FractionCSP3(mol)
            props.rdkit_version = rdBase.rdkitVersion
            props.save()
            return props

        except ImportError:
            return None
        except Exception:
            return None


class MolecularPropertyThreshold(models.Model):
    """
    Configurable RAG (Red/Amber/Green) thresholds for molecular properties.

    Allows users/admins to define warning and danger thresholds for
    drug-likeness properties. Used for visual highlighting in aggregation tables.
    """

    PROPERTY_CHOICES = [
        ('molecular_weight', 'Molecular Weight'),
        ('heavy_atom_count', 'Heavy Atom Count'),
        ('hbd', 'H-Bond Donors'),
        ('hba', 'H-Bond Acceptors'),
        ('clogp', 'cLogP'),
        ('tpsa', 'TPSA'),
        ('rotatable_bonds', 'Rotatable Bonds'),
        ('fraction_sp3', 'Fraction sp3'),
    ]

    DIRECTION_CHOICES = [
        ('above', 'Above threshold is bad'),  # e.g., MW > 500
        ('below', 'Below threshold is bad'),  # e.g., Fsp3 < 0.25
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    property_name = models.CharField(
        max_length=32,
        choices=PROPERTY_CHOICES,
        unique=True,
        help_text="Which molecular property this threshold applies to"
    )
    direction = models.CharField(
        max_length=8,
        choices=DIRECTION_CHOICES,
        default='above',
        help_text="Whether values above or below threshold are concerning"
    )
    amber_threshold = models.FloatField(
        help_text="Value at which property turns amber (warning)"
    )
    red_threshold = models.FloatField(
        help_text="Value at which property turns red (danger)"
    )
    enabled = models.BooleanField(
        default=True,
        help_text="Whether RAG coloring is active for this property"
    )

    class Meta:
        verbose_name = 'Property Threshold'
        verbose_name_plural = 'Property Thresholds'
        ordering = ['property_name']

    def __str__(self):
        return f'{self.get_property_name_display()} threshold'

    def get_rag_status(self, value: float) -> str:
        """
        Determine RAG status for a given value.

        Returns: 'green', 'amber', or 'red'
        """
        if value is None or not self.enabled:
            return 'green'

        if self.direction == 'above':
            if value >= self.red_threshold:
                return 'red'
            elif value >= self.amber_threshold:
                return 'amber'
            return 'green'
        else:  # below
            if value <= self.red_threshold:
                return 'red'
            elif value <= self.amber_threshold:
                return 'amber'
            return 'green'

    @classmethod
    def get_default_thresholds(cls) -> list[dict]:
        """
        Return Lipinski Rule of 5 defaults for seeding the database.
        """
        return [
            {'property_name': 'molecular_weight', 'direction': 'above',
             'amber_threshold': 450, 'red_threshold': 500},
            {'property_name': 'clogp', 'direction': 'above',
             'amber_threshold': 4, 'red_threshold': 5},
            {'property_name': 'hbd', 'direction': 'above',
             'amber_threshold': 4, 'red_threshold': 5},
            {'property_name': 'hba', 'direction': 'above',
             'amber_threshold': 8, 'red_threshold': 10},
            {'property_name': 'tpsa', 'direction': 'above',
             'amber_threshold': 120, 'red_threshold': 140},
            {'property_name': 'rotatable_bonds', 'direction': 'above',
             'amber_threshold': 8, 'red_threshold': 10},
            {'property_name': 'fraction_sp3', 'direction': 'below',
             'amber_threshold': 0.3, 'red_threshold': 0.2},
        ]


class ScaffoldExtension(models.Model):
    """
    Runtime extension to the curated NLP scaffold catalog
    (``apps/compounds/nlp/substructures.py``).

    Lets chemists add fragments the seed module doesn't cover, scoped
    either to a specific Target (project-local nomenclature like
    "Series 1" / "the deletion compounds") or shared across the whole
    instance (target=null). Consulted by ``resolve_scaffold`` alongside
    the Python seed; project-scoped entries override shared, shared
    overrides seed.

    Storage is SMARTS, not SMILES — the matching engine
    (``RDKit.Chem.MolFromSmarts`` + ``HasSubstructMatch``) is the same
    one the seed module uses, so positional degeneracy / atom-list
    tolerance / ring-size constraints are all expressible in the
    ``smarts`` field. The ``notes`` column captures the human gloss.
    """

    name = models.CharField(
        max_length=64,
        help_text='Canonical name chemists will type ("indazole", "Series 1")',
    )
    smarts = models.CharField(
        max_length=512,
        help_text='SMARTS pattern; matched via RDKit HasSubstructMatch',
    )
    aliases = models.JSONField(
        default=list, blank=True,
        help_text='Additional typed forms that map to this entry (e.g. plurals, abbreviations)',
    )
    target = models.ForeignKey(
        Target,
        on_delete=models.CASCADE,
        null=True, blank=True,
        related_name='scaffold_extensions',
        help_text='If set, this entry only resolves under that target/project. Null = shared across all projects.',
    )
    notes = models.TextField(
        blank=True,
        help_text='Human gloss explaining what the fragment captures (e.g. "Series 1: 5-6 fused heterocycles, atoms 2/4 may be C or N")',
    )

    # Audit
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True, blank=True,
        related_name='created_scaffold_extensions',
    )
    created_at = models.DateTimeField(auto_now_add=True)
    source_prompt = models.TextField(
        blank=True,
        help_text="Original NLP prompt that triggered the addition, if any",
    )
    llm_generated = models.BooleanField(
        default=False,
        help_text='True when the SMARTS was proposed by the LLM (slice 18+); False when typed by the chemist directly.',
    )

    class Meta:
        verbose_name = 'Scaffold extension'
        verbose_name_plural = 'Scaffold extensions'
        # One name per scope: ('benzofuran', null) is shared, ('Series 1', ARd_pk)
        # is project-local. The same name can exist twice across scopes
        # — that's the per-project-overrides-shared semantics.
        unique_together = [('name', 'target')]
        ordering = ['target__name', 'name']

    def __str__(self) -> str:
        scope = self.target.name if self.target_id else 'shared'
        return f'{self.name} ({scope})'
