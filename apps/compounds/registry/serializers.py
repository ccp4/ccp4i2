"""
Registry Serializers

DRF serializers for compound registration models.
"""

from datetime import timedelta

from django.urls import reverse
from django.utils import timezone
from rest_framework import serializers

from compounds.validators import validate_qc_file
from .models import Supplier, Target, Compound, Batch, BatchQCFile, CompoundTemplate


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


class SavedAggregationViewSerializer(serializers.Serializer):
    """Serializer for validating saved aggregation view configuration."""
    protocol_names = serializers.ListField(
        child=serializers.CharField(max_length=256),
        required=False,
        default=list
    )
    compound_search = serializers.CharField(
        required=False,
        default='',
        allow_blank=True
    )
    output_format = serializers.ChoiceField(
        choices=['compact', 'medium', 'long'],
        default='compact'
    )
    aggregations = serializers.ListField(
        child=serializers.ChoiceField(choices=['geomean', 'count', 'stdev', 'list']),
        required=False,
        default=['geomean', 'count']
    )
    status = serializers.ChoiceField(
        choices=['valid', 'invalid', 'unassigned', ''],
        required=False,
        default='valid',
        allow_blank=True
    )


class TargetSerializer(serializers.ModelSerializer):
    parent_name = serializers.CharField(source='parent.name', read_only=True)
    compound_count = serializers.IntegerField(source='compounds.count', read_only=True)
    assay_count = serializers.SerializerMethodField()
    has_recent_compounds = serializers.SerializerMethodField()
    has_recent_assays = serializers.SerializerMethodField()
    latest_activity = serializers.DateTimeField(read_only=True, required=False)
    image = ProtectedFileField(url_name='target-image', model_field='id', url_kwarg='target_id')

    class Meta:
        model = Target
        fields = [
            'id', 'name', 'parent', 'parent_name', 'created_at',
            'compound_count', 'assay_count',
            'has_recent_compounds', 'has_recent_assays',
            'latest_activity',
            'image',
            'saved_aggregation_view',
        ]

    def get_assay_count(self, obj):
        """Count assays for this target."""
        return obj.assays.count()

    def get_has_recent_compounds(self, obj):
        """Check if any compounds were registered in the last 7 days."""
        cutoff = timezone.now() - timedelta(days=7)
        return obj.compounds.filter(registered_at__gte=cutoff).exists()

    def get_has_recent_assays(self, obj):
        """Check if any assays were created in the last 7 days."""
        cutoff = timezone.now() - timedelta(days=7)
        return obj.assays.filter(created_at__gte=cutoff).exists()


class TargetDetailSerializer(TargetSerializer):
    """Extended serializer with nested compounds."""
    compounds_count = serializers.IntegerField(source='compounds.count', read_only=True)
    children = TargetSerializer(many=True, read_only=True)

    class Meta(TargetSerializer.Meta):
        fields = TargetSerializer.Meta.fields + ['compounds_count', 'children']


class DashboardCompoundSerializer(serializers.ModelSerializer):
    """Compact compound serializer for dashboard carousel."""
    formatted_id = serializers.CharField(read_only=True)

    class Meta:
        model = Compound
        fields = ['id', 'formatted_id', 'smiles', 'registered_at', 'molecular_weight']


class DashboardAssaySerializer(serializers.Serializer):
    """Assay serializer for dashboard carousel."""
    id = serializers.UUIDField()
    protocol_name = serializers.CharField()
    created_at = serializers.DateTimeField()
    data_series_count = serializers.IntegerField()


class DashboardProjectSerializer(serializers.Serializer):
    """CCP4i2 project serializer for dashboard carousel."""
    id = serializers.IntegerField()
    name = serializers.CharField()
    last_access = serializers.DateTimeField()
    job_count = serializers.IntegerField()
    matching_compound_ids = serializers.ListField(child=serializers.CharField())


class TargetDashboardSerializer(TargetSerializer):
    """Target serializer with dashboard data (recent compounds, assays)."""
    recent_compounds = DashboardCompoundSerializer(many=True, read_only=True)
    recent_assays = DashboardAssaySerializer(many=True, read_only=True)

    class Meta(TargetSerializer.Meta):
        fields = TargetSerializer.Meta.fields + ['recent_compounds', 'recent_assays']


class CompoundListSerializer(serializers.ModelSerializer):
    """Compact serializer for list views."""
    formatted_id = serializers.CharField(read_only=True)
    target_name = serializers.CharField(source='target.name', read_only=True)
    supplier_name = serializers.CharField(source='supplier.name', read_only=True)
    # Use annotated field from queryset (avoids N+1), falls back to count for non-list views
    batch_count = serializers.SerializerMethodField()

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

    def get_batch_count(self, obj):
        """Use annotated batch_count if available, otherwise compute."""
        if hasattr(obj, 'batch_count'):
            return obj.batch_count
        return obj.batches.count()


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
    formatted_id = serializers.CharField(read_only=True)

    class Meta:
        model = Compound
        fields = [
            'id', 'formatted_id',
            'target', 'smiles', 'stereo_comment',
            'supplier', 'supplier_ref',
            'labbook_number', 'page_number', 'compound_number',
            'comments',
        ]
        read_only_fields = ['id', 'formatted_id']

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
    """Serializer for reading batch QC files - returns protected download URLs."""
    filename = serializers.CharField(read_only=True)
    # Use protected URL for file instead of direct storage URL
    file = ProtectedFileField(url_name='batch-qc-file', model_field='id', url_kwarg='qc_file_id')

    class Meta:
        model = BatchQCFile
        fields = ['id', 'batch', 'file', 'filename', 'comments', 'uploaded_at']
        read_only_fields = ['uploaded_at']


class BatchQCFileCreateSerializer(serializers.ModelSerializer):
    """Serializer for creating/uploading batch QC files."""

    class Meta:
        model = BatchQCFile
        fields = ['id', 'batch', 'file', 'comments']
        read_only_fields = ['id']

    def validate_file(self, value):
        """Validate file extension and size."""
        validate_qc_file(value)
        return value


class CompoundTemplateSerializer(serializers.ModelSerializer):
    target_name = serializers.CharField(source='target.name', read_only=True)

    class Meta:
        model = CompoundTemplate
        fields = ['id', 'target', 'target_name', 'name', 'mol2d', 'svg_file']
