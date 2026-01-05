"""
Registry ViewSets

DRF ViewSets for compound registration models with reversion support.
"""

import reversion
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import viewsets, filters, status
from rest_framework.decorators import action
from rest_framework.response import Response
from reversion.models import Version

from .models import Supplier, Target, Compound, Batch, BatchQCFile, CompoundTemplate
from .serializers import (
    SupplierSerializer,
    TargetSerializer,
    TargetDetailSerializer,
    CompoundListSerializer,
    CompoundDetailSerializer,
    CompoundCreateSerializer,
    BatchSerializer,
    BatchQCFileSerializer,
    CompoundTemplateSerializer,
)


class ReversionMixin:
    """Mixin to add reversion support to ViewSets."""

    def perform_create(self, serializer):
        with reversion.create_revision():
            instance = serializer.save()
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment(f"Created via API")
            return instance

    def perform_update(self, serializer):
        with reversion.create_revision():
            instance = serializer.save()
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment(f"Updated via API")
            return instance

    def perform_destroy(self, instance):
        with reversion.create_revision():
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment(f"Deleted via API")
            instance.delete()

    @action(detail=True, methods=['get'])
    def history(self, request, pk=None):
        """Get version history for this object."""
        obj = self.get_object()
        versions = Version.objects.get_for_object(obj)
        return Response([{
            'id': v.id,
            'date': v.revision.date_created,
            'user': v.revision.user.email if v.revision.user else None,
            'comment': v.revision.comment,
        } for v in versions[:50]])  # Limit to 50 most recent


class SupplierViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Suppliers."""
    queryset = Supplier.objects.all()
    serializer_class = SupplierSerializer
    filter_backends = [filters.SearchFilter, filters.OrderingFilter]
    search_fields = ['name', 'initials']
    ordering_fields = ['name']
    ordering = ['name']


class TargetViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Targets (drug discovery campaigns)."""
    queryset = Target.objects.all()
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['parent']
    search_fields = ['name']
    ordering_fields = ['name', 'created_at']
    ordering = ['name']

    def get_serializer_class(self):
        if self.action == 'retrieve':
            return TargetDetailSerializer
        return TargetSerializer

    @action(detail=True, methods=['get'])
    def compounds(self, request, pk=None):
        """List compounds for this target."""
        target = self.get_object()
        compounds = target.compounds.all()
        serializer = CompoundListSerializer(compounds, many=True)
        return Response(serializer.data)


class CompoundViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Compounds."""
    queryset = Compound.objects.select_related('target', 'supplier', 'registered_by')
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['target', 'supplier', 'stereo_comment']
    search_fields = ['smiles', 'supplier_ref', 'comments']
    ordering_fields = ['reg_number', 'registered_at', 'molecular_weight']
    ordering = ['-reg_number']

    def get_serializer_class(self):
        if self.action == 'create':
            return CompoundCreateSerializer
        if self.action == 'retrieve':
            return CompoundDetailSerializer
        return CompoundListSerializer

    @action(detail=False, methods=['get'])
    def by_reg_number(self, request):
        """Lookup compound by registration number."""
        reg_number = request.query_params.get('reg_number')
        if not reg_number:
            return Response(
                {'error': 'reg_number parameter required'},
                status=status.HTTP_400_BAD_REQUEST
            )
        try:
            compound = Compound.objects.get(reg_number=int(reg_number))
            serializer = CompoundDetailSerializer(compound)
            return Response(serializer.data)
        except (ValueError, Compound.DoesNotExist):
            return Response(
                {'error': 'Compound not found'},
                status=status.HTTP_404_NOT_FOUND
            )

    @action(detail=True, methods=['get'])
    def batches(self, request, pk=None):
        """List batches for this compound."""
        compound = self.get_object()
        batches = compound.batches.all()
        serializer = BatchSerializer(batches, many=True)
        return Response(serializer.data)


class BatchViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Batches."""
    queryset = Batch.objects.select_related('compound', 'supplier')
    serializer_class = BatchSerializer
    filter_backends = [DjangoFilterBackend, filters.OrderingFilter]
    filterset_fields = ['compound', 'supplier']
    ordering_fields = ['batch_number', 'registered_at']
    ordering = ['compound', 'batch_number']

    @action(detail=True, methods=['get'])
    def qc_files(self, request, pk=None):
        """List QC files for this batch."""
        batch = self.get_object()
        files = batch.qc_files.all()
        serializer = BatchQCFileSerializer(files, many=True)
        return Response(serializer.data)


class BatchQCFileViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Batch QC Files."""
    queryset = BatchQCFile.objects.select_related('batch', 'batch__compound')
    serializer_class = BatchQCFileSerializer
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ['batch']


class CompoundTemplateViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Compound Templates."""
    queryset = CompoundTemplate.objects.select_related('target')
    serializer_class = CompoundTemplateSerializer
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['target']
    search_fields = ['name']
