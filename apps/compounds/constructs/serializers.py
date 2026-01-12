"""
Constructs Serializers

DRF serializers for construct database models.
Follows tiered serializer pattern (List, Detail, Create variants).
"""

from rest_framework import serializers

from .models import (
    ConstructProject,
    Plasmid,
    Protein,
    ProteinSynonym,
    ProteinUse,
    Cassette,
    CassetteUse,
    SequencingResult,
    ExpressionTagType,
    Protease,
    ExpressionTagLocation,
    ExpressionTag,
)


# =============================================================================
# Reference Data Serializers
# =============================================================================

class ExpressionTagTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = ExpressionTagType
        fields = ['id', 'name', 'created_at']
        read_only_fields = ['created_at']


class ProteaseSerializer(serializers.ModelSerializer):
    class Meta:
        model = Protease
        fields = ['id', 'name', 'created_at']
        read_only_fields = ['created_at']


class ExpressionTagLocationSerializer(serializers.ModelSerializer):
    class Meta:
        model = ExpressionTagLocation
        fields = ['id', 'name', 'created_at']
        read_only_fields = ['created_at']


# =============================================================================
# Construct Project Serializers
# =============================================================================

class ConstructProjectSerializer(serializers.ModelSerializer):
    parent_name = serializers.CharField(source='parent.name', read_only=True)
    plasmid_count = serializers.IntegerField(source='plasmids.count', read_only=True)
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)

    class Meta:
        model = ConstructProject
        fields = [
            'id', 'name', 'parent', 'parent_name',
            'plasmid_count', 'created_at', 'updated_at',
            'created_by', 'created_by_email',
        ]
        read_only_fields = ['created_at', 'updated_at']


class ConstructProjectDetailSerializer(ConstructProjectSerializer):
    """Extended serializer with nested children."""
    children = ConstructProjectSerializer(many=True, read_only=True)

    class Meta(ConstructProjectSerializer.Meta):
        fields = ConstructProjectSerializer.Meta.fields + ['children']


class ConstructProjectCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating construct projects."""

    class Meta:
        model = ConstructProject
        fields = ['name', 'parent']

    def create(self, validated_data):
        request = self.context.get('request')
        if request and request.user.is_authenticated:
            validated_data['created_by'] = request.user
        return super().create(validated_data)


# =============================================================================
# Protein Serializers
# =============================================================================

class ProteinSynonymSerializer(serializers.ModelSerializer):
    protein_uniprot_id = serializers.CharField(source='protein.uniprot_id', read_only=True)

    class Meta:
        model = ProteinSynonym
        fields = ['id', 'name', 'protein', 'protein_uniprot_id', 'created_at']
        read_only_fields = ['created_at']


class ProteinListSerializer(serializers.ModelSerializer):
    """Compact serializer for protein list views."""
    synonym_count = serializers.IntegerField(source='synonyms.count', read_only=True)
    cassette_count = serializers.IntegerField(source='cassettes.count', read_only=True)

    class Meta:
        model = Protein
        fields = ['id', 'uniprot_id', 'synonym_count', 'cassette_count', 'created_at']


class ProteinDetailSerializer(serializers.ModelSerializer):
    """Full serializer with synonyms."""
    synonyms = ProteinSynonymSerializer(many=True, read_only=True)
    cassette_count = serializers.IntegerField(source='cassettes.count', read_only=True)
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)

    class Meta:
        model = Protein
        fields = [
            'id', 'uniprot_id', 'synonyms', 'cassette_count',
            'created_at', 'updated_at', 'created_by', 'created_by_email',
        ]
        read_only_fields = ['created_at', 'updated_at']


class ProteinCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating proteins."""

    class Meta:
        model = Protein
        fields = ['uniprot_id']

    def create(self, validated_data):
        request = self.context.get('request')
        if request and request.user.is_authenticated:
            validated_data['created_by'] = request.user
        return super().create(validated_data)


class ProteinUseSerializer(serializers.ModelSerializer):
    protein_uniprot_id = serializers.CharField(source='protein.uniprot_id', read_only=True)
    project_name = serializers.CharField(source='project.name', read_only=True)

    class Meta:
        model = ProteinUse
        fields = [
            'id', 'protein', 'protein_uniprot_id',
            'project', 'project_name', 'created_at',
        ]
        read_only_fields = ['created_at']


# =============================================================================
# Cassette Serializers
# =============================================================================

class CassetteListSerializer(serializers.ModelSerializer):
    """Compact serializer for cassette list views."""
    protein_uniprot_id = serializers.CharField(source='protein.uniprot_id', read_only=True)
    display_name = serializers.CharField(read_only=True)
    plasmid_count = serializers.IntegerField(source='plasmid_uses.count', read_only=True)

    class Meta:
        model = Cassette
        fields = [
            'id', 'protein', 'protein_uniprot_id',
            'start', 'end', 'display_name', 'plasmid_count', 'created_at',
        ]


class CassetteDetailSerializer(serializers.ModelSerializer):
    """Full serializer with relationships."""
    protein_uniprot_id = serializers.CharField(source='protein.uniprot_id', read_only=True)
    display_name = serializers.CharField(read_only=True)
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)

    class Meta:
        model = Cassette
        fields = [
            'id', 'protein', 'protein_uniprot_id',
            'start', 'end', 'display_name',
            'created_at', 'updated_at', 'created_by', 'created_by_email',
        ]
        read_only_fields = ['created_at', 'updated_at']


class CassetteCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating cassettes."""

    class Meta:
        model = Cassette
        fields = ['protein', 'start', 'end']

    def create(self, validated_data):
        request = self.context.get('request')
        if request and request.user.is_authenticated:
            validated_data['created_by'] = request.user
        return super().create(validated_data)


# =============================================================================
# Expression Tag Serializer
# =============================================================================

class ExpressionTagSerializer(serializers.ModelSerializer):
    """Serializer for expression tags on cassette uses."""
    expression_tag_type_name = serializers.CharField(
        source='expression_tag_type.name', read_only=True
    )
    protease_name = serializers.CharField(source='protease.name', read_only=True)
    location_name = serializers.CharField(source='location.name', read_only=True)

    class Meta:
        model = ExpressionTag
        fields = [
            'id', 'expression_tag_type', 'expression_tag_type_name',
            'protease', 'protease_name',
            'cassette_use', 'location', 'location_name',
            'created_at',
        ]
        read_only_fields = ['created_at']


# =============================================================================
# Cassette Use Serializers
# =============================================================================

class CassetteUseListSerializer(serializers.ModelSerializer):
    """Compact serializer for cassette use list views."""
    cassette_display = serializers.CharField(source='cassette.display_name', read_only=True)
    plasmid_formatted_id = serializers.CharField(
        source='plasmid.formatted_id', read_only=True
    )
    expression_tag_count = serializers.IntegerField(
        source='expression_tags.count', read_only=True
    )

    class Meta:
        model = CassetteUse
        fields = [
            'id', 'cassette', 'cassette_display',
            'plasmid', 'plasmid_formatted_id',
            'expression_tag_count', 'created_at',
        ]


class CassetteUseDetailSerializer(serializers.ModelSerializer):
    """Full serializer with expression tags."""
    cassette_display = serializers.CharField(source='cassette.display_name', read_only=True)
    plasmid_formatted_id = serializers.CharField(
        source='plasmid.formatted_id', read_only=True
    )
    expression_tags = ExpressionTagSerializer(many=True, read_only=True)
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)

    class Meta:
        model = CassetteUse
        fields = [
            'id', 'cassette', 'cassette_display',
            'plasmid', 'plasmid_formatted_id',
            'alignment_file', 'expression_tags',
            'created_at', 'updated_at', 'created_by', 'created_by_email',
        ]
        read_only_fields = ['created_at', 'updated_at']


class CassetteUseCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating cassette uses."""

    class Meta:
        model = CassetteUse
        fields = ['cassette', 'plasmid', 'alignment_file']

    def create(self, validated_data):
        request = self.context.get('request')
        if request and request.user.is_authenticated:
            validated_data['created_by'] = request.user
        return super().create(validated_data)


# =============================================================================
# Sequencing Result Serializers
# =============================================================================

class SequencingResultSerializer(serializers.ModelSerializer):
    """Serializer for sequencing results."""
    plasmid_formatted_id = serializers.CharField(
        source='plasmid.formatted_id', read_only=True
    )
    cassette_display = serializers.CharField(
        source='cassette_use.cassette.display_name', read_only=True
    )
    filename = serializers.CharField(read_only=True)

    class Meta:
        model = SequencingResult
        fields = [
            'id', 'cassette_use', 'cassette_display',
            'plasmid', 'plasmid_formatted_id',
            'file', 'filename', 'created_at',
        ]
        read_only_fields = ['created_at']


class SequencingResultCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating sequencing results."""

    class Meta:
        model = SequencingResult
        fields = ['cassette_use', 'plasmid', 'file']

    def create(self, validated_data):
        request = self.context.get('request')
        if request and request.user.is_authenticated:
            validated_data['created_by'] = request.user
        return super().create(validated_data)


# =============================================================================
# Plasmid Serializers
# =============================================================================

class PlasmidListSerializer(serializers.ModelSerializer):
    """Compact serializer for plasmid list views."""
    formatted_id = serializers.CharField(read_only=True)
    project_name = serializers.CharField(source='project.name', read_only=True)
    parent_formatted_id = serializers.CharField(
        source='parent.formatted_id', read_only=True
    )
    cassette_count = serializers.IntegerField(
        source='cassette_uses.count', read_only=True
    )
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)

    class Meta:
        model = Plasmid
        fields = [
            'id', 'ncn_id', 'formatted_id', 'name',
            'project', 'project_name',
            'parent', 'parent_formatted_id',
            'cassette_count', 'created_at', 'created_by_email',
        ]


class PlasmidDetailSerializer(serializers.ModelSerializer):
    """Full serializer with all fields and relationships."""
    formatted_id = serializers.CharField(read_only=True)
    project_name = serializers.CharField(source='project.name', read_only=True)
    parent_formatted_id = serializers.CharField(
        source='parent.formatted_id', read_only=True
    )
    cassette_uses = CassetteUseDetailSerializer(many=True, read_only=True)
    sequencing_results = SequencingResultSerializer(many=True, read_only=True)
    created_by_email = serializers.CharField(source='created_by.email', read_only=True)
    genbank_file_url = serializers.SerializerMethodField()

    class Meta:
        model = Plasmid
        fields = [
            'id', 'ncn_id', 'formatted_id', 'name',
            'project', 'project_name',
            'parent', 'parent_formatted_id',
            'genbank_file', 'genbank_file_url',
            'cassette_uses', 'sequencing_results',
            'created_at', 'updated_at', 'created_by', 'created_by_email',
        ]
        read_only_fields = ['ncn_id', 'formatted_id', 'created_at', 'updated_at']

    def get_genbank_file_url(self, obj):
        """Return URL for downloading the GenBank file."""
        if obj.genbank_file:
            request = self.context.get('request')
            if request:
                return request.build_absolute_uri(obj.genbank_file.url)
            return obj.genbank_file.url
        return None


class PlasmidCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating plasmids."""
    formatted_id = serializers.CharField(read_only=True)

    class Meta:
        model = Plasmid
        fields = ['id', 'formatted_id', 'name', 'project', 'parent', 'genbank_file']
        read_only_fields = ['id', 'formatted_id']

    def create(self, validated_data):
        request = self.context.get('request')
        if request and request.user.is_authenticated:
            validated_data['created_by'] = request.user
        return super().create(validated_data)
