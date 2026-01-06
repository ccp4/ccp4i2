"""
Assays Serializers

DRF serializers for assay models.
"""

from rest_framework import serializers

from .models import (
    DilutionSeries,
    Protocol,
    ProtocolDocument,
    Assay,
    DataSeries,
    AnalysisResult,
    Hypothesis,
)
from compounds.registry.serializers import CompoundListSerializer


class DilutionSeriesSerializer(serializers.ModelSerializer):
    display_name = serializers.CharField(source='__str__', read_only=True)

    class Meta:
        model = DilutionSeries
        fields = ['id', 'concentrations', 'unit', 'display_name']


class ProtocolDocumentSerializer(serializers.ModelSerializer):
    filename = serializers.SerializerMethodField()

    class Meta:
        model = ProtocolDocument
        fields = ['id', 'protocol', 'file', 'filename', 'created_by', 'created_at']
        read_only_fields = ['created_by', 'created_at']

    def get_filename(self, obj):
        if obj.file:
            from pathlib import Path
            return Path(obj.file.name).name
        return None


class ProtocolSerializer(serializers.ModelSerializer):
    preferred_dilutions_display = serializers.CharField(
        source='preferred_dilutions.__str__', read_only=True
    )
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)

    class Meta:
        model = Protocol
        fields = [
            'id', 'name', 'analysis_method', 'pherastar_table',
            'preferred_dilutions', 'preferred_dilutions_display',
            'created_by', 'created_by_email', 'created_at',
            'comments',
        ]
        read_only_fields = ['created_by', 'created_at']


class ProtocolDetailSerializer(ProtocolSerializer):
    """Extended serializer with documents."""
    documents = ProtocolDocumentSerializer(many=True, read_only=True)
    assays_count = serializers.IntegerField(source='assays.count', read_only=True)

    class Meta(ProtocolSerializer.Meta):
        fields = ProtocolSerializer.Meta.fields + ['documents', 'assays_count']


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
    analysis_status = serializers.CharField(source='analysis.status', read_only=True)
    analysis_kpi = serializers.SerializerMethodField()
    # Include dilution_series and analysis for chart rendering in tables
    dilution_series = DilutionSeriesSerializer(read_only=True)
    analysis = AnalysisResultSerializer(read_only=True)

    class Meta:
        model = DataSeries
        fields = [
            'id', 'assay', 'compound', 'compound_formatted_id', 'compound_name',
            'row', 'start_column', 'end_column',
            'dilution_series', 'extracted_data',
            'analysis', 'analysis_status', 'analysis_kpi',
        ]

    def get_analysis_kpi(self, obj):
        if obj.analysis:
            return obj.analysis.kpi_value
        return None


class DataSeriesDetailSerializer(serializers.ModelSerializer):
    """Full serializer with all data."""
    compound_formatted_id = serializers.CharField(
        source='compound.formatted_id', read_only=True
    )
    analysis = AnalysisResultSerializer(read_only=True)
    dilution_series = DilutionSeriesSerializer(read_only=True)

    class Meta:
        model = DataSeries
        fields = [
            'id', 'assay', 'compound', 'compound_formatted_id', 'compound_name',
            'row', 'start_column', 'end_column',
            'dilution_series', 'extracted_data', 'skip_points',
            'analysis', 'svg_file',
        ]


class AssayListSerializer(serializers.ModelSerializer):
    """Compact serializer for list views."""
    protocol_name = serializers.CharField(source='protocol.name', read_only=True)
    target_name = serializers.CharField(source='target.name', read_only=True)
    data_filename = serializers.CharField(read_only=True)
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)
    data_series_count = serializers.IntegerField(source='data_series.count', read_only=True)

    class Meta:
        model = Assay
        fields = [
            'id', 'protocol', 'protocol_name',
            'target', 'target_name',
            'data_filename', 'created_by_email', 'created_at',
            'data_series_count',
        ]


class AssayDetailSerializer(serializers.ModelSerializer):
    """Full serializer with nested data."""
    protocol_name = serializers.CharField(source='protocol.name', read_only=True)
    target_name = serializers.CharField(source='target.name', read_only=True)
    data_filename = serializers.CharField(read_only=True)
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)
    data_series = DataSeriesListSerializer(many=True, read_only=True)

    class Meta:
        model = Assay
        fields = [
            'id', 'protocol', 'protocol_name',
            'target', 'target_name',
            'data_file', 'data_filename',
            'labbook_number', 'page_number',
            'created_by', 'created_by_email', 'created_at',
            'comments',
            'data_series',
        ]
        read_only_fields = ['created_by', 'created_at']


class AssayCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating assays."""

    class Meta:
        model = Assay
        fields = [
            'protocol', 'target', 'data_file',
            'labbook_number', 'page_number', 'comments',
        ]

    def create(self, validated_data):
        request = self.context.get('request')
        if request and request.user.is_authenticated:
            validated_data['created_by'] = request.user
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
