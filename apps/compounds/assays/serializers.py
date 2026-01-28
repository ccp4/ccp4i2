"""
Assays Serializers

DRF serializers for assay models.
"""

from django.urls import reverse
from rest_framework import serializers

from compounds.validators import validate_protocol_document, validate_assay_data_file
from .models import (
    DilutionSeries,
    PlateLayout,
    Protocol,
    ProtocolDocument,
    Assay,
    DataSeries,
    AnalysisResult,
    Hypothesis,
    FittingMethod,
)
from compounds.registry.models import Target
from compounds.registry.serializers import CompoundListSerializer


class ProtectedFileField(serializers.Field):
    """
    Custom field that returns a protected API URL for file downloads
    instead of the direct storage URL.

    This ensures files are served through authenticated endpoints.
    """

    def __init__(self, url_name, model_field='id', url_kwarg=None, **kwargs):
        """
        Args:
            url_name: Name of the URL pattern for the protected endpoint
            model_field: Field on the model to get the ID from (default: 'id')
            url_kwarg: Name of the URL kwarg (default: same as model_field)
        """
        self.url_name = url_name
        self.model_field = model_field
        self.url_kwarg = url_kwarg or model_field
        kwargs['read_only'] = True
        super().__init__(**kwargs)

    def to_representation(self, value):
        """Return the protected URL if file exists, None otherwise."""
        if not value:
            return None

        # Django's FieldFile has an 'instance' attribute pointing to the model
        # This works correctly during both single and list serialization
        obj = getattr(value, 'instance', None)

        # Fallback to parent's instance for single object serialization
        if obj is None and self.parent:
            instance = getattr(self.parent, 'instance', None)
            if instance and not (hasattr(instance, '__iter__') and not isinstance(instance, dict)):
                obj = instance

        if obj is None:
            return None

        obj_id = getattr(obj, self.model_field, None)
        if not obj_id:
            return None

        # Build the protected URL - return path through the Next.js proxy
        # The Django URL is /api/compounds/media/... but frontend needs /api/proxy/compounds/media/...
        # This ensures URLs work both in Docker (where Django is at server:8000) and locally
        url = reverse(f'compounds:{self.url_name}', kwargs={self.url_kwarg: obj_id})
        # Convert /api/compounds/... to /api/proxy/compounds/...
        if url.startswith('/api/compounds/'):
            url = '/api/proxy/compounds/' + url[len('/api/compounds/'):]
        return url


class FittingMethodSerializer(serializers.ModelSerializer):
    """Serializer for FittingMethod - exposes metadata without the script code."""

    class Meta:
        model = FittingMethod
        fields = ['id', 'name', 'slug', 'version', 'description', 'is_active', 'is_builtin']
        read_only_fields = ['id', 'slug', 'version', 'is_builtin']


class FittingMethodDetailSerializer(serializers.ModelSerializer):
    """Full serializer for FittingMethod including script code for editing."""

    class Meta:
        model = FittingMethod
        fields = [
            'id', 'name', 'slug', 'version', 'description',
            'script', 'input_schema', 'output_schema',
            'is_active', 'is_builtin',
            'created_at', 'updated_at',
        ]
        read_only_fields = ['id', 'slug', 'is_builtin', 'created_at', 'updated_at']


class DilutionSeriesSerializer(serializers.ModelSerializer):
    display_name = serializers.CharField(source='__str__', read_only=True)

    class Meta:
        model = DilutionSeries
        fields = ['id', 'concentrations', 'unit', 'display_name']


class PlateLayoutSerializer(serializers.ModelSerializer):
    """Serializer for PlateLayout - reusable plate configurations."""
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)
    protocols_count = serializers.SerializerMethodField()
    # Convenience field to expose plate_format from config
    plate_format = serializers.IntegerField(source='config.plate_format', read_only=True)

    class Meta:
        model = PlateLayout
        fields = [
            'id', 'name', 'description', 'config',
            'plate_format', 'is_builtin',
            'created_by', 'created_by_email',
            'created_at', 'updated_at',
            'protocols_count',
        ]
        read_only_fields = ['id', 'created_by', 'created_at', 'updated_at', 'is_builtin']

    def get_protocols_count(self, obj):
        return obj.protocols.count()

    def validate_config(self, value):
        """Convert null to empty dict."""
        return value if value is not None else {}


class ProtocolDocumentSerializer(serializers.ModelSerializer):
    """Serializer for reading protocol documents - returns protected download URLs."""
    filename = serializers.SerializerMethodField()
    # Use protected URL for file instead of direct storage URL
    file = ProtectedFileField(url_name='protocol-document-file', model_field='id', url_kwarg='document_id')

    class Meta:
        model = ProtocolDocument
        fields = ['id', 'protocol', 'file', 'filename', 'created_by', 'created_at']
        read_only_fields = ['created_by', 'created_at']

    def get_filename(self, obj):
        if obj.file:
            from pathlib import Path
            return Path(obj.file.name).name
        return None


class ProtocolDocumentCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating/uploading protocol documents."""

    class Meta:
        model = ProtocolDocument
        fields = ['id', 'protocol', 'file']
        read_only_fields = ['id']

    def validate_file(self, value):
        """Validate file extension and size."""
        validate_protocol_document(value)
        return value


class ProtocolSerializer(serializers.ModelSerializer):
    preferred_dilutions_display = serializers.SerializerMethodField()
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)
    modified_by_email = serializers.CharField(source='modified_by.email', read_only=True)
    fitting_method_name = serializers.CharField(
        source='fitting_method.name', read_only=True
    )
    target_name = serializers.CharField(source='target.name', read_only=True)
    # Plate layout FK with denormalized fields for display
    plate_layout = serializers.PrimaryKeyRelatedField(
        queryset=PlateLayout.objects.all(),
        allow_null=True,
        required=False
    )
    plate_layout_name = serializers.CharField(source='plate_layout.name', read_only=True)
    plate_layout_config = serializers.JSONField(source='plate_layout.config', read_only=True)
    # Explicitly allow null for ForeignKey fields
    fitting_method = serializers.PrimaryKeyRelatedField(
        queryset=FittingMethod.objects.all(),
        allow_null=True,
        required=False
    )
    preferred_dilutions = serializers.PrimaryKeyRelatedField(
        queryset=DilutionSeries.objects.all(),
        allow_null=True,
        required=False
    )
    target = serializers.PrimaryKeyRelatedField(
        queryset=Target.objects.all(),
        allow_null=True,
        required=False
    )
    # fitting_parameters doesn't allow null in DB, so convert null to empty dict
    fitting_parameters = serializers.JSONField(
        required=False,
        allow_null=True
    )
    # List view fields
    assays_count = serializers.IntegerField(source='assays.count', read_only=True)
    has_recent_assays = serializers.SerializerMethodField()

    class Meta:
        model = Protocol
        fields = [
            'id', 'name', 'import_type', 'analysis_method',
            'fitting_method', 'fitting_method_name',
            'target', 'target_name',
            'plate_layout', 'plate_layout_name', 'plate_layout_config',
            'fitting_parameters',
            'preferred_dilutions', 'preferred_dilutions_display',
            'created_by', 'created_by_email',
            'modified_by', 'modified_by_email',
            'created_at', 'modified_at',
            'comments',
            'assays_count', 'has_recent_assays',
        ]
        read_only_fields = ['created_by', 'created_at', 'modified_at']

    def validate_fitting_parameters(self, value):
        """Convert null to empty dict since DB column doesn't allow NULL."""
        return value if value is not None else {}

    def get_preferred_dilutions_display(self, obj):
        """Return string representation of dilution series, or None if not set."""
        if obj.preferred_dilutions:
            return str(obj.preferred_dilutions)
        return None

    def get_has_recent_assays(self, obj):
        """Return True if any assays were created in the last 7 days."""
        from django.utils import timezone
        from datetime import timedelta
        seven_days_ago = timezone.now() - timedelta(days=7)
        return obj.assays.filter(created_at__gte=seven_days_ago).exists()


class ProtocolDetailSerializer(ProtocolSerializer):
    """Extended serializer with documents."""
    documents = ProtocolDocumentSerializer(many=True, read_only=True)

    class Meta(ProtocolSerializer.Meta):
        fields = ProtocolSerializer.Meta.fields + ['documents']


class AnalysisResultSerializer(serializers.ModelSerializer):
    kpi_value = serializers.SerializerMethodField()

    class Meta:
        model = AnalysisResult
        fields = ['id', 'status', 'results', 'kpi_value']

    def get_kpi_value(self, obj):
        return obj.kpi_value


class DataSeriesListSerializer(serializers.ModelSerializer):
    """Compact serializer for list views with chart data."""
    compound_formatted_id = serializers.CharField(
        source='compound.formatted_id', read_only=True
    )
    batch_number = serializers.IntegerField(
        source='batch.batch_number', read_only=True
    )
    analysis_status = serializers.CharField(source='analysis.status', read_only=True)
    analysis_kpi = serializers.SerializerMethodField()
    # Include dilution_series and analysis for chart rendering in tables
    # Use SerializerMethodField to implement fallback to protocol's preferred_dilutions
    dilution_series = serializers.SerializerMethodField()
    analysis = AnalysisResultSerializer(read_only=True)

    class Meta:
        model = DataSeries
        fields = [
            'id', 'assay', 'compound', 'compound_formatted_id',
            'batch', 'batch_number',
            'compound_name',
            'row', 'start_column', 'end_column',
            'dilution_series', 'extracted_data',
            'analysis', 'analysis_status', 'analysis_kpi',
            'plot_image',
        ]

    def get_dilution_series(self, obj):
        """Return dilution_series, falling back to protocol's preferred_dilutions."""
        dilution = obj.dilution_series
        if not dilution and obj.assay and obj.assay.protocol:
            dilution = obj.assay.protocol.preferred_dilutions
        if dilution:
            return DilutionSeriesSerializer(dilution).data
        return None

    def get_analysis_kpi(self, obj):
        if obj.analysis:
            return obj.analysis.kpi_value
        return None


class DataSeriesDetailSerializer(serializers.ModelSerializer):
    """Full serializer with all data."""
    compound_formatted_id = serializers.CharField(
        source='compound.formatted_id', read_only=True
    )
    batch_number = serializers.IntegerField(
        source='batch.batch_number', read_only=True
    )
    analysis = AnalysisResultSerializer(read_only=True)
    # Use SerializerMethodField to implement fallback to protocol's preferred_dilutions
    dilution_series = serializers.SerializerMethodField()

    class Meta:
        model = DataSeries
        fields = [
            'id', 'assay', 'compound', 'compound_formatted_id',
            'batch', 'batch_number',
            'compound_name',
            'row', 'start_column', 'end_column',
            'dilution_series', 'extracted_data', 'skip_points',
            'analysis', 'plot_image',
        ]

    def get_dilution_series(self, obj):
        """Return dilution_series, falling back to protocol's preferred_dilutions."""
        dilution = obj.dilution_series
        if not dilution and obj.assay and obj.assay.protocol:
            dilution = obj.assay.protocol.preferred_dilutions
        if dilution:
            return DilutionSeriesSerializer(dilution).data
        return None


class DataSeriesCreateSerializer(serializers.ModelSerializer):
    """
    Serializer for creating data series from frontend extraction.

    Supports batch specification in two ways:
    1. Embedded in compound_name with "/" separator: "NCL-00026042/1"
    2. Explicit batch field (takes precedence)

    When compound_name includes a batch suffix (e.g., "NCL-00026042/1" or
    "NCL-00026042/001"), both compound and batch are resolved automatically.
    """

    class Meta:
        model = DataSeries
        fields = [
            'id', 'assay', 'compound', 'batch', 'compound_name',
            'row', 'start_column', 'end_column',
            'dilution_series', 'extracted_data',
        ]
        read_only_fields = ['id']

    def create(self, validated_data):
        # Try to match compound (and optionally batch) by name if not provided
        if not validated_data.get('compound') and validated_data.get('compound_name'):
            from compounds.utils import resolve_compound_with_batch
            compound, batch = resolve_compound_with_batch(validated_data['compound_name'])
            if compound:
                validated_data['compound'] = compound
            # Only set batch if not explicitly provided and we found one
            if batch and not validated_data.get('batch'):
                validated_data['batch'] = batch

        # Inherit dilution_series from protocol if not explicitly provided
        if not validated_data.get('dilution_series'):
            assay = validated_data.get('assay')
            if assay and assay.protocol and assay.protocol.preferred_dilutions:
                validated_data['dilution_series'] = assay.protocol.preferred_dilutions

        return super().create(validated_data)


class DataSeriesUpdateSerializer(serializers.ModelSerializer):
    """Serializer for updating data series fields including dilution_series and batch."""

    class Meta:
        model = DataSeries
        fields = [
            'id', 'batch', 'dilution_series', 'skip_points',
        ]
        read_only_fields = ['id']


class AssayListSerializer(serializers.ModelSerializer):
    """Compact serializer for list views."""
    protocol_name = serializers.CharField(source='protocol.name', read_only=True)
    target_name = serializers.CharField(source='target.name', read_only=True)
    data_filename = serializers.CharField(read_only=True)
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)
    modified_by_email = serializers.CharField(source='modified_by.email', read_only=True)
    data_series_count = serializers.IntegerField(source='data_series.count', read_only=True)

    class Meta:
        model = Assay
        fields = [
            'id', 'protocol', 'protocol_name',
            'target', 'target_name',
            'data_filename',
            'created_by_email', 'modified_by_email',
            'created_at', 'modified_at',
            'data_series_count',
        ]


class AssayDetailSerializer(serializers.ModelSerializer):
    """Full serializer with nested data."""
    protocol_name = serializers.CharField(source='protocol.name', read_only=True)
    target_name = serializers.CharField(source='target.name', read_only=True)
    data_filename = serializers.CharField(read_only=True)
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)
    modified_by_email = serializers.CharField(source='modified_by.email', read_only=True)
    data_series = DataSeriesListSerializer(many=True, read_only=True)
    # Use protected URL for data file instead of direct storage URL
    data_file = ProtectedFileField(url_name='assay-data-file', model_field='id', url_kwarg='assay_id')

    class Meta:
        model = Assay
        fields = [
            'id', 'protocol', 'protocol_name',
            'target', 'target_name',
            'data_file', 'data_filename',
            'labbook_number', 'page_number',
            'created_by', 'created_by_email',
            'modified_by', 'modified_by_email',
            'created_at', 'modified_at',
            'comments',
            'data_series',
        ]
        read_only_fields = ['created_by', 'created_at', 'modified_at']


class AssayCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating assays."""

    class Meta:
        model = Assay
        fields = [
            'id', 'protocol', 'target', 'data_file',
            'labbook_number', 'page_number', 'comments',
        ]
        read_only_fields = ['id']

    def validate_data_file(self, value):
        """Validate file extension and size."""
        validate_assay_data_file(value)
        return value

    def create(self, validated_data):
        request = self.context.get('request')
        if request and request.user.is_authenticated:
            validated_data['created_by'] = request.user
        # Default target from protocol if not provided explicitly
        if 'target' not in validated_data or validated_data['target'] is None:
            protocol = validated_data.get('protocol')
            if protocol and protocol.target:
                validated_data['target'] = protocol.target
        return super().create(validated_data)


class HypothesisListSerializer(serializers.ModelSerializer):
    """Compact serializer for list views."""
    target_name = serializers.CharField(source='target.name', read_only=True)
    product_formatted_id = serializers.CharField(
        source='product_compound.formatted_id', read_only=True
    )

    class Meta:
        model = Hypothesis
        fields = [
            'id', 'target', 'target_name',
            'smiles', 'status',
            'product_compound', 'product_formatted_id',
            'created_at', 'updated_at',
        ]


class HypothesisDetailSerializer(serializers.ModelSerializer):
    """Full serializer with all fields."""
    target_name = serializers.CharField(source='target.name', read_only=True)
    parent_smiles = serializers.CharField(source='parent.smiles', read_only=True)
    product_compound_detail = CompoundListSerializer(
        source='product_compound', read_only=True
    )

    class Meta:
        model = Hypothesis
        fields = [
            'id', 'target', 'target_name',
            'parent', 'parent_smiles',
            'smiles', 'rationale', 'model_url',
            'status', 'completion_notes',
            'product_compound', 'product_compound_detail',
            'svg_file',
            'created_at', 'updated_at',
        ]


# Aggregation serializers

class PredicatesSerializer(serializers.Serializer):
    """Serializer for aggregation filter predicates."""
    targets = serializers.ListField(
        child=serializers.UUIDField(),
        required=False,
        default=list,
        help_text="List of target UUIDs to filter by"
    )
    compounds = serializers.ListField(
        child=serializers.UUIDField(),
        required=False,
        default=list,
        help_text="List of compound UUIDs to filter by"
    )
    compound_search = serializers.CharField(
        required=False,
        default='',
        help_text="Text search for compound formatted_id (e.g., NCL-00026)"
    )
    protocols = serializers.ListField(
        child=serializers.UUIDField(),
        required=False,
        default=list,
        help_text="List of protocol UUIDs to filter by"
    )
    status = serializers.ChoiceField(
        choices=['valid', 'invalid', 'unassigned', ''],
        required=False,
        default='valid',
        help_text="Analysis status filter (default: valid)"
    )


class AggregationRequestSerializer(serializers.Serializer):
    """Serializer for aggregation request validation."""
    predicates = PredicatesSerializer(required=False, default=dict)
    output_format = serializers.ChoiceField(
        choices=['compact', 'long'],
        default='compact',
        help_text="Output format: 'compact' (one row per compound) or 'long' (one row per measurement)"
    )
    aggregations = serializers.ListField(
        child=serializers.ChoiceField(choices=['geomean', 'count', 'stdev', 'list']),
        default=['geomean', 'count'],
        help_text="Aggregation functions to compute"
    )


class TableOfValuesRowSerializer(serializers.Serializer):
    """Serializer for a single row in Table of Values import."""
    # All fields are dynamic - we just validate it's a dict
    pass


class TableOfValuesImportSerializer(serializers.Serializer):
    """
    Serializer for Table of Values bulk import.

    Validates the import request structure for pre-analyzed data.
    """
    compound_column = serializers.CharField(
        help_text="Name of the column containing compound identifiers"
    )
    kpi_column = serializers.CharField(
        help_text="Name of the column containing KPI names (e.g., 'KPI')"
    )
    image_column = serializers.CharField(
        required=False,
        allow_blank=True,
        default='',
        help_text="Optional: Name of the column containing image filenames"
    )
    kpi_unit = serializers.CharField(
        required=False,
        allow_blank=True,
        allow_null=True,
        default=None,
        help_text="Optional: Override the inferred KPI unit (e.g., 'nM', 'uM', '%')"
    )
    data = serializers.ListField(
        child=serializers.DictField(),
        help_text="List of row objects from the spreadsheet"
    )

    def validate(self, attrs):
        """Validate that required columns exist in the data."""
        data = attrs.get('data', [])
        compound_column = attrs.get('compound_column')
        kpi_column = attrs.get('kpi_column')

        if not data:
            raise serializers.ValidationError("No data rows provided")

        # Check first row has required columns
        first_row = data[0]
        if compound_column not in first_row:
            raise serializers.ValidationError(
                f"Compound column '{compound_column}' not found in data"
            )
        if kpi_column not in first_row:
            raise serializers.ValidationError(
                f"KPI column '{kpi_column}' not found in data"
            )

        # Validate all rows have the same KPI value
        kpi_values = set(row.get(kpi_column) for row in data if row.get(kpi_column))
        if len(kpi_values) > 1:
            raise serializers.ValidationError(
                f"All rows must have the same KPI value. Found: {', '.join(str(v) for v in kpi_values)}"
            )

        # Validate the KPI value is a column name that exists
        if kpi_values:
            kpi_value = next(iter(kpi_values))
            if kpi_value not in first_row:
                raise serializers.ValidationError(
                    f"KPI value '{kpi_value}' is not a column name in the data"
                )

        return attrs
