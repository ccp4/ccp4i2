"""
Compound Assays Models

Models for assay protocols, experiments, data series, and analysis results.
Migrated from legacy AssayCompounds app with cleaned-up naming.

Model mapping from legacy:
    Protocol -> Protocol (unchanged)
    ProtocolDocument -> ProtocolDocument (unchanged)
    DilutionSeries -> DilutionSeries (unchanged)
    Experiment -> Assay
    DataSeries -> DataSeries (unchanged)
    Analysis -> AnalysisResult
    Hypothesis -> Hypothesis (unchanged)
"""

import uuid
from pathlib import Path

from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import models
from django.utils.functional import cached_property

# Import from registry for cross-app relationships
from compounds.registry.models import Target, Compound

User = get_user_model()


class DilutionSeries(models.Model):
    """
    Standard concentration series for dose-response experiments.

    Defines the concentrations used in assay plates, typically a
    serial dilution series like [100, 30, 10, 3, 1, 0.3, 0.1] µM.
    """

    UNIT_CHOICES = [
        ('nM', 'nanomolar'),
        ('uM', 'micromolar'),
        ('mM', 'millimolar'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    concentrations = models.JSONField(
        default=list,
        help_text="List of concentrations, e.g., [100, 30, 10, 3, 1, 0.3]"
    )
    unit = models.CharField(max_length=10, choices=UNIT_CHOICES, default='nM')

    class Meta:
        verbose_name = 'Dilution Series'
        verbose_name_plural = 'Dilution Series'

    def __str__(self):
        concs = ','.join(str(c) for c in self.concentrations)
        return f"{concs} {self.unit}"


class Protocol(models.Model):
    """
    Assay protocol definition.

    Defines the experimental method, analysis approach, and default
    parameters for running an assay.
    """

    ANALYSIS_METHOD_CHOICES = [
        ('hill_langmuir', 'Hill-Langmuir'),
        ('hill_langmuir_fix_hill', 'Hill-Langmuir (fixed Hill coefficient)'),
        ('hill_langmuir_fix_hill_minmax', 'Hill-Langmuir (fixed Hill and min/max)'),
        ('hill_langmuir_fix_minmax', 'Hill-Langmuir (fixed min/max)'),
        ('ms_intact', 'MS-Intact'),
        ('table_of_values', 'Table of values'),
    ]

    PHERASTAR_TABLE_CHOICES = [
        ('none', 'None'),
        ('table_1', 'Table 1'),
        ('table_2', 'Table 2'),
        ('table_3', 'Table 3'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=256)
    analysis_method = models.CharField(
        max_length=50,
        choices=ANALYSIS_METHOD_CHOICES,
        default='hill_langmuir'
    )
    preferred_dilutions = models.ForeignKey(
        DilutionSeries,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        help_text="Default dilution series for experiments using this protocol"
    )
    pherastar_table = models.CharField(
        max_length=32,
        choices=PHERASTAR_TABLE_CHOICES,
        default='table_3',
        blank=True,
        null=True,
        help_text="PHERAstar plate reader table selection"
    )

    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True
    )
    created_at = models.DateTimeField(auto_now_add=True)
    comments = models.TextField(blank=True, null=True)

    class Meta:
        ordering = ['name']
        verbose_name = 'Protocol'
        verbose_name_plural = 'Protocols'

    def __str__(self):
        return self.name


def _protocol_doc_path(instance, filename):
    """Generate upload path for protocol documents."""
    return Path('protocols') / str(instance.protocol.id) / filename


class ProtocolDocument(models.Model):
    """
    Supporting documentation for a protocol.

    Protocols can have multiple attached documents (SOPs, references, etc.)
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    protocol = models.ForeignKey(
        Protocol,
        on_delete=models.CASCADE,
        related_name='documents'
    )
    file = models.FileField(upload_to=_protocol_doc_path)
    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True
    )
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        verbose_name = 'Protocol Document'
        verbose_name_plural = 'Protocol Documents'

    def __str__(self):
        filename = Path(self.file.name).name if self.file else 'unnamed'
        return f'{self.protocol.name}: {filename}'


def _assay_data_path(instance, filename):
    """Generate upload path for assay data files."""
    return Path('assays') / str(instance.id) / filename


class Assay(models.Model):
    """
    An assay experiment (plate read, etc.).

    Renamed from 'Experiment' for clarity. An Assay runs a Protocol
    on one or more compounds, producing DataSeries results.

    Note: The target field allows testing compounds against a target
    different from their registration target (e.g., selectivity testing).
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    protocol = models.ForeignKey(
        Protocol,
        on_delete=models.PROTECT,
        related_name='assays'
    )
    target = models.ForeignKey(
        Target,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='assays',
        help_text="Target tested (may differ from compound's registered target)"
    )

    data_file = models.FileField(
        upload_to=_assay_data_path,
        help_text="Raw data file (Excel, CSV, etc.)"
    )
    labbook_number = models.IntegerField(null=True, blank=True)
    page_number = models.IntegerField(null=True, blank=True)

    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True
    )
    created_at = models.DateTimeField(auto_now_add=True)
    comments = models.TextField(blank=True, null=True)

    class Meta:
        ordering = ['-created_at']
        verbose_name = 'Assay'
        verbose_name_plural = 'Assays'

    def __str__(self):
        date_str = self.created_at.strftime('%Y-%m-%d') if self.created_at else 'unknown'
        return f'{self.protocol.name} - {date_str}'

    @property
    def data_filename(self):
        """Extract filename from data_file path."""
        return Path(self.data_file.name).name if self.data_file else None


def _series_svg_path(instance, filename):
    """Generate upload path for data series plot SVGs."""
    assay_id = instance.assay.id
    return Path('assays') / str(assay_id) / 'plots' / f'{instance.id}.svg'


class DataSeries(models.Model):
    """
    A single compound's dose-response data within an assay.

    Each DataSeries represents one compound's results from an Assay,
    including the raw data points and fitted analysis.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    assay = models.ForeignKey(
        Assay,
        on_delete=models.CASCADE,
        related_name='data_series'
    )
    compound = models.ForeignKey(
        Compound,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='assay_results',
        help_text="Linked compound (matched by name/barcode)"
    )
    compound_name = models.CharField(
        max_length=256,
        blank=True,
        null=True,
        help_text="Original name from data file (for matching)"
    )

    # Position in source file
    row = models.IntegerField()
    start_column = models.IntegerField()
    end_column = models.IntegerField()

    # Data
    dilution_series = models.ForeignKey(
        DilutionSeries,
        on_delete=models.SET_NULL,
        null=True,
        blank=True
    )
    extracted_data = models.JSONField(
        default=dict,
        help_text="Parsed data points from source file"
    )
    skip_points = models.JSONField(
        default=list,
        help_text="Indices of points to exclude from fitting"
    )

    # Results
    analysis = models.OneToOneField(
        'AnalysisResult',
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='data_series'
    )
    svg_file = models.ImageField(
        upload_to=_series_svg_path,
        blank=True,
        null=True,
        help_text="Fitted curve plot"
    )

    class Meta:
        ordering = ['compound_name']
        verbose_name = 'Data Series'
        verbose_name_plural = 'Data Series'

    def __str__(self):
        return f'{self.compound_name} in {self.assay}'


class AnalysisResult(models.Model):
    """
    Fitted results for a data series.

    Stores the outcome of curve fitting (EC50, Hill coefficient, etc.)
    and validation status.
    """

    STATUS_CHOICES = [
        ('valid', 'Valid'),
        ('invalid', 'Invalid'),
        ('unassigned', 'Unassigned'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    status = models.CharField(
        max_length=20,
        choices=STATUS_CHOICES,
        default='unassigned'
    )
    results = models.JSONField(
        default=dict,
        help_text="Fitted parameters: {EC50, Hill, minVal, maxVal, KPI, ...}"
    )

    class Meta:
        verbose_name = 'Analysis Result'
        verbose_name_plural = 'Analysis Results'

    def __str__(self):
        if self.status != 'valid':
            return self.status.upper()

        # Try to format the KPI value
        kpi_key = self.results.get('KPI')
        if kpi_key and kpi_key in self.results:
            value = self.results[kpi_key]
            if isinstance(value, (int, float)):
                return f"{value:.2f}"
            return str(value) if value else 'N/A'

        return 'Valid'

    @property
    def kpi_value(self):
        """Extract the key performance indicator value."""
        kpi_key = self.results.get('KPI')
        if kpi_key and kpi_key in self.results:
            return self.results[kpi_key]
        return None


def _hypothesis_svg_path(instance, filename):
    """Generate upload path for hypothesis SVG images."""
    return Path('hypotheses') / 'svg' / f'{instance.id}.svg'


class Hypothesis(models.Model):
    """
    Compound design hypothesis for a target.

    Tracks proposed compound ideas through their lifecycle from
    initial design to synthesis and testing.
    """

    STATUS_CHOICES = [
        ('pending', 'Pending'),
        ('rejected', 'Rejected'),
        ('chemistry', 'In Chemistry'),
        ('shelved', 'Shelved'),
        ('made', 'Made'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    target = models.ForeignKey(
        Target,
        on_delete=models.CASCADE,
        related_name='hypotheses'
    )
    parent = models.ForeignKey(
        'self',
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='children',
        help_text="Parent hypothesis this was derived from"
    )

    smiles = models.CharField(max_length=1024, help_text="Proposed structure")
    rationale = models.TextField(blank=True, null=True, help_text="Design rationale")
    model_url = models.CharField(
        max_length=2048,
        blank=True,
        null=True,
        help_text="Link to 3D model or docking pose"
    )

    status = models.CharField(
        max_length=20,
        choices=STATUS_CHOICES,
        default='pending'
    )
    product_compound = models.ForeignKey(
        Compound,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='from_hypotheses',
        help_text="Resulting compound if hypothesis was synthesized"
    )
    completion_notes = models.TextField(blank=True, null=True)

    svg_file = models.ImageField(
        upload_to=_hypothesis_svg_path,
        blank=True,
        null=True
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        ordering = ['-created_at']
        verbose_name = 'Hypothesis'
        verbose_name_plural = 'Hypotheses'

    def __str__(self):
        return f'{self.target.name}: {self.smiles[:50]}...' if len(self.smiles) > 50 else f'{self.target.name}: {self.smiles}'
