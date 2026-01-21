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

from django.contrib.auth import get_user_model
from django.db import models
from django.db.models.signals import post_delete
from django.dispatch import receiver
from django.utils.functional import cached_property

# Import from registry for cross-app relationships
from compounds.registry.models import Target, Compound
from compounds.utils import delete_file_field

User = get_user_model()


def _fitting_script_path(instance, filename):
    """Generate upload path for fitting script files."""
    return Path('fitting_scripts') / f'{instance.id}_{filename}'


class FittingMethod(models.Model):
    """
    Versioned curve fitting script for analyzing dose-response data.

    Each FittingMethod contains a Python script that accepts extracted
    assay data and returns fitted parameters (IC50, Hill coefficient, etc.).

    Scripts are admin-managed for security. The script interface:

        def fit(input_data: dict) -> dict:
            '''
            Args:
                input_data: {
                    "concentrations": [10000, 3333, 1111, ...],
                    "responses": [95.2, 87.1, 62.3, ...],
                    "controls": {"max": 100.0, "min": 2.3},
                    "parameters": {}  # method-specific options
                }

            Returns:
                {
                    "ic50": 234.5,
                    "hill_slope": 1.2,
                    "top": 98.7,
                    "bottom": 3.1,
                    "r_squared": 0.994,
                    "curve_points": [[x1, y1], ...],
                    "flags": [],
                    "kpi": "ic50",  # which result is the primary KPI
                }
            '''
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=128, help_text="Display name for the fitting method")
    slug = models.SlugField(
        max_length=64,
        unique=True,
        help_text="URL-safe identifier, e.g., 'four-parameter-logistic'"
    )
    version = models.CharField(
        max_length=32,
        default='1.0.0',
        help_text="Semantic version string"
    )
    description = models.TextField(
        blank=True,
        null=True,
        help_text="Detailed description of the fitting algorithm"
    )

    # Script content - stored as text for easy versioning/viewing
    script = models.TextField(
        help_text="Python script implementing the fit() function"
    )

    # Schema definitions for validation and UI generation
    input_schema = models.JSONField(
        default=dict,
        blank=True,
        help_text="JSON Schema describing expected input parameters"
    )
    output_schema = models.JSONField(
        default=dict,
        blank=True,
        help_text="JSON Schema describing output structure"
    )

    # Metadata
    is_active = models.BooleanField(
        default=True,
        help_text="Inactive methods are hidden from selection but preserved for historical data"
    )
    is_builtin = models.BooleanField(
        default=False,
        help_text="Built-in methods shipped with the system"
    )

    created_by = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='created_fitting_methods'
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        ordering = ['name', '-version']
        verbose_name = 'Fitting Method'
        verbose_name_plural = 'Fitting Methods'
        unique_together = [['slug', 'version']]

    def __str__(self):
        return f"{self.name} v{self.version}"


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

    Defines the experimental method, analysis approach, plate layout,
    and default parameters for running an assay.

    The plate_layout field stores a JSON structure defining:
    - Plate format (96, 384, 1536)
    - Control well positions (max/min)
    - Sample region boundaries
    - Dilution direction and orientation
    - Replicate pattern
    - Compound source/mapping rules

    Example plate_layout:
    {
        "plate_format": 384,
        "controls": {
            "max": {"columns": [1, 2], "rows": ["A", "B"]},
            "min": {"columns": [23, 24], "rows": ["A", "B"]}
        },
        "sample_region": {
            "start_column": 3,
            "end_column": 22,
            "start_row": "A",
            "end_row": "P"
        },
        "dilution": {
            "direction": "horizontal",
            "num_concentrations": 10
        },
        "replicate": {
            "count": 2,
            "pattern": "adjacent_rows"
        },
        "compound_source": {
            "type": "row_order",
            "id_column": "Compound_ID"
        }
    }
    """

    ANALYSIS_METHOD_CHOICES = [
        ('hill_langmuir', 'Hill-Langmuir'),
        ('hill_langmuir_fix_hill', 'Hill-Langmuir (fixed Hill coefficient)'),
        ('hill_langmuir_fix_hill_minmax', 'Hill-Langmuir (fixed Hill and min/max)'),
        ('hill_langmuir_fix_minmax', 'Hill-Langmuir (fixed min/max)'),
        ('ms_intact', 'MS-Intact'),
        ('table_of_values', 'Table of values'),
        ('pharmaron_adme', 'Pharmaron ADME'),
    ]

    PLATE_FORMAT_CHOICES = [
        (24, '24-well'),
        (96, '96-well'),
        (384, '384-well'),
        (1536, '1536-well'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=256)

    # Legacy analysis method - kept for backward compatibility
    analysis_method = models.CharField(
        max_length=50,
        choices=ANALYSIS_METHOD_CHOICES,
        default='hill_langmuir',
        help_text="Legacy analysis method (use fitting_method for new protocols)"
    )

    # New: Link to versioned fitting script
    fitting_method = models.ForeignKey(
        FittingMethod,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='protocols',
        help_text="Curve fitting algorithm for analyzing data from this protocol"
    )

    preferred_dilutions = models.ForeignKey(
        DilutionSeries,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        help_text="Default dilution series for experiments using this protocol"
    )

    # New: Plate layout configuration
    plate_layout = models.JSONField(
        default=dict,
        blank=True,
        help_text="Plate layout template: control positions, sample region, replicate pattern"
    )

    # Default fitting parameters (can be overridden per-assay)
    fitting_parameters = models.JSONField(
        default=dict,
        blank=True,
        help_text="Default parameters passed to fitting script"
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

    def get_effective_fitting_method(self):
        """Return the fitting method to use, falling back to legacy analysis_method."""
        if self.fitting_method:
            return self.fitting_method
        # Could return a default FittingMethod based on analysis_method
        return None


def _protocol_doc_path(instance, filename):
    """Generate upload path for protocol documents.

    Legacy pattern: AssayCompounds/Protocols/Protocol_{uuid}/{filename}
    """
    return Path('AssayCompounds') / 'Protocols' / f'Protocol_{instance.protocol.id}' / filename


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


@receiver(post_delete, sender=ProtocolDocument)
def _protocol_document_post_delete(sender, instance, **kwargs):
    """Delete associated file when protocol document is deleted."""
    delete_file_field(instance.file)


def _assay_data_path(instance, filename):
    """Generate upload path for assay data files.

    Legacy pattern: AssayCompounds/Experiments/Experiment_{uuid}/{filename}
    Note: 'Assay' was called 'Experiment' in legacy.
    """
    return Path('AssayCompounds') / 'Experiments' / f'Experiment_{instance.id}' / filename


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


@receiver(post_delete, sender=Assay)
def _assay_post_delete(sender, instance, **kwargs):
    """Delete associated data file when assay is deleted."""
    delete_file_field(instance.data_file)


def _series_plot_path(instance, filename):
    """Generate upload path for data series plot images.

    Legacy pattern: Same directory as parent Experiment's dataFile.
    AssayCompounds/Experiments/Experiment_{uuid}/{filename}

    Supports various image formats (SVG, PNG, JPG, etc.).
    For Hill-Langmuir: auto-generated plots pass '{id}.svg' as filename
    For Table-Of-Values: preserves original filename for matching via
                         analysis.results['Image File']
    """
    assay_id = instance.assay.id
    return Path('AssayCompounds') / 'Experiments' / f'Experiment_{assay_id}' / filename


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
    plot_image = models.ImageField(
        upload_to=_series_plot_path,
        blank=True,
        null=True,
        help_text="Fitted curve plot image"
    )

    class Meta:
        ordering = ['compound_name']
        verbose_name = 'Data Series'
        verbose_name_plural = 'Data Series'

    def __str__(self):
        return f'{self.compound_name} in {self.assay}'


@receiver(post_delete, sender=DataSeries)
def _data_series_post_delete(sender, instance, **kwargs):
    """Delete associated plot image when data series is deleted."""
    delete_file_field(instance.plot_image)


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
    """Generate upload path for hypothesis SVG images.

    Legacy pattern: AssayCompounds/svg/{uuid}.svg
    """
    return Path('AssayCompounds') / 'svg' / f'{instance.id}.svg'


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
