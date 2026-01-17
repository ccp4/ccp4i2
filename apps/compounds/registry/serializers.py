"""
Registry Serializers

DRF serializers for compound registration models.
"""

from django.urls import reverse
from rest_framework import serializers

from .models import Supplier, Target, Compound, Batch, BatchQCFile, CompoundTemplate


class ProtectedFileField(serializers.Field):
    """
    Custom field that returns a protected API URL for file downloads
    instead of the direct storage URL.

    This ensures files are served through authenticated endpoints.
    """

    def __init__(self, url_name, id_field='id', **kwargs):
        """
        Args:
            url_name: Name of the URL pattern for the protected endpoint
            id_field: Field on the model to use as the URL argument
        """
        self.url_name = url_name
        self.id_field = id_field
        kwargs['read_only'] = True
        super().__init__(**kwargs)

    def to_representation(self, value):
        """Return the protected URL if file exists, None otherwise."""
        if not value:
            return None

        # Get the parent object's ID
        obj = self.parent.instance
        if hasattr(obj, '__iter__') and not isinstance(obj, dict):
            return None

        obj_id = getattr(obj, self.id_field, None)
        if not obj_id:
            return None

        # Build the protected URL
        request = self.context.get('request')
        url = reverse(f'compounds:{self.url_name}', kwargs={f'{self.id_field}': obj_id})

        if request:
            return request.build_absolute_uri(url)
        return url


class SupplierSerializer(serializers.ModelSerializer):
    is_current_user = serializers.SerializerMethodField()
    compound_count = serializers.IntegerField(source='compounds.count', read_only=True)
    batch_count = serializers.SerializerMethodField()

    class Meta:
        model = Supplier
        fields = ['id', 'name', 'initials', 'user', 'is_current_user', 'compound_count', 'batch_count']
        read_only_fields = ['is_current_user', 'compound_count', 'batch_count']

    def get_is_current_user(self, obj):
        """Check if this supplier is linked to the current user."""
        request = self.context.get('request')
        if request and request.user.is_authenticated and obj.user:
            return obj.user.id == request.user.id
        return False

    def get_batch_count(self, obj):
        """Count batches from this supplier."""
        return Batch.objects.filter(supplier=obj).count()


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
    # Use protected URL for file instead of direct storage URL
    file = ProtectedFileField(url_name='batch-qc-file', id_field='id')

    class Meta:
        model = BatchQCFile
        fields = ['id', 'batch', 'file', 'filename', 'comments', 'uploaded_at']
        read_only_fields = ['uploaded_at']


class CompoundTemplateSerializer(serializers.ModelSerializer):
    target_name = serializers.CharField(source='target.name', read_only=True)

    class Meta:
        model = CompoundTemplate
        fields = ['id', 'target', 'target_name', 'name', 'mol2d', 'svg_file']
