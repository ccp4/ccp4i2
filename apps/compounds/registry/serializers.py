"""
Registry Serializers

DRF serializers for compound registration models.
"""

from rest_framework import serializers

from .models import Supplier, Target, Compound, Batch, BatchQCFile, CompoundTemplate


class SupplierSerializer(serializers.ModelSerializer):
    class Meta:
        model = Supplier
        fields = ['id', 'name', 'initials']


class TargetSerializer(serializers.ModelSerializer):
    parent_name = serializers.CharField(source='parent.name', read_only=True)
    compound_count = serializers.IntegerField(source='compounds.count', read_only=True)

    class Meta:
        model = Target
        fields = ['id', 'name', 'parent', 'parent_name', 'created_at', 'compound_count']


class TargetDetailSerializer(TargetSerializer):
    """Extended serializer with nested compounds."""
    compounds_count = serializers.IntegerField(source='compounds.count', read_only=True)
    children = TargetSerializer(many=True, read_only=True)

    class Meta(TargetSerializer.Meta):
        fields = TargetSerializer.Meta.fields + ['compounds_count', 'children']


class CompoundListSerializer(serializers.ModelSerializer):
    """Compact serializer for list views."""
    formatted_id = serializers.CharField(read_only=True)
    target_name = serializers.CharField(source='target.name', read_only=True)
    supplier_name = serializers.CharField(source='supplier.name', read_only=True)
    batch_count = serializers.IntegerField(source='batches.count', read_only=True)

    class Meta:
        model = Compound
        fields = [
            'id', 'reg_number', 'formatted_id',
            'target', 'target_name',
            'smiles', 'molecular_weight', 'stereo_comment',
            'supplier', 'supplier_name', 'supplier_ref',
            'registered_at',
            'batch_count',
        ]


class CompoundDetailSerializer(serializers.ModelSerializer):
    """Full serializer with all fields."""
    formatted_id = serializers.CharField(read_only=True)
    barcode = serializers.CharField(read_only=True)
    target_name = serializers.CharField(source='target.name', read_only=True)
    supplier_name = serializers.CharField(source='supplier.name', read_only=True)
    registered_by_email = serializers.CharField(
        source='registered_by.email', read_only=True
    )
    batch_count = serializers.IntegerField(source='batches.count', read_only=True)

    class Meta:
        model = Compound
        fields = [
            'id', 'reg_number', 'formatted_id', 'barcode',
            'target', 'target_name',
            'smiles', 'rdkit_smiles', 'inchi', 'molecular_weight', 'stereo_comment',
            'supplier', 'supplier_name', 'supplier_ref',
            'labbook_number', 'page_number', 'compound_number',
            'registered_by', 'registered_by_email', 'legacy_registered_by',
            'registered_at', 'modified_at',
            'comments', 'svg_file',
            'batch_count',
        ]
        read_only_fields = ['reg_number', 'formatted_id', 'barcode', 'registered_at', 'modified_at']


class CompoundCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating compounds."""

    class Meta:
        model = Compound
        fields = [
            'target', 'smiles', 'stereo_comment',
            'supplier', 'supplier_ref',
            'labbook_number', 'page_number', 'compound_number',
            'comments',
        ]

    def create(self, validated_data):
        # Set registered_by from request user
        request = self.context.get('request')
        if request and request.user.is_authenticated:
            validated_data['registered_by'] = request.user
        return super().create(validated_data)


class BatchSerializer(serializers.ModelSerializer):
    compound_formatted_id = serializers.CharField(
        source='compound.formatted_id', read_only=True
    )
    supplier_name = serializers.CharField(source='supplier.name', read_only=True)
    qc_file_count = serializers.IntegerField(source='qc_files.count', read_only=True)

    class Meta:
        model = Batch
        fields = [
            'id', 'compound', 'compound_formatted_id', 'batch_number',
            'supplier', 'supplier_name', 'supplier_ref',
            'labbook_number', 'page_number',
            'amount', 'salt_code', 'molecular_weight',
            'registered_at', 'comments',
            'qc_file_count',
        ]
        read_only_fields = ['batch_number', 'registered_at']


class BatchQCFileSerializer(serializers.ModelSerializer):
    filename = serializers.CharField(read_only=True)

    class Meta:
        model = BatchQCFile
        fields = ['id', 'batch', 'file', 'filename', 'comments', 'uploaded_at']
        read_only_fields = ['uploaded_at']


class CompoundTemplateSerializer(serializers.ModelSerializer):
    target_name = serializers.CharField(source='target.name', read_only=True)

    class Meta:
        model = CompoundTemplate
        fields = ['id', 'target', 'target_name', 'name', 'mol2d', 'svg_file']
