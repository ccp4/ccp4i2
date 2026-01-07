"""
Assays ViewSets

DRF ViewSets for assay models with reversion support.
"""

import reversion
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import viewsets, filters, status
from rest_framework.decorators import action
from rest_framework.response import Response
from reversion.models import Version

from .analysis import analyse_assay, analyse_data_series

from .models import (
    DilutionSeries,
    Protocol,
    ProtocolDocument,
    Assay,
    DataSeries,
    AnalysisResult,
    Hypothesis,
)
from .serializers import (
    DilutionSeriesSerializer,
    ProtocolSerializer,
    ProtocolDetailSerializer,
    ProtocolDocumentSerializer,
    AssayListSerializer,
    AssayDetailSerializer,
    AssayCreateSerializer,
    DataSeriesListSerializer,
    DataSeriesDetailSerializer,
    DataSeriesCreateSerializer,
    AnalysisResultSerializer,
    HypothesisListSerializer,
    HypothesisDetailSerializer,
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
        } for v in versions[:50]])


class DilutionSeriesViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Dilution Series."""
    queryset = DilutionSeries.objects.all()
    serializer_class = DilutionSeriesSerializer
    filter_backends = [filters.OrderingFilter]
    ordering = ['unit']


class ProtocolViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Protocols."""
    queryset = Protocol.objects.select_related('preferred_dilutions', 'created_by')
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['analysis_method']
    search_fields = ['name', 'comments']
    ordering_fields = ['name', 'created_at']
    ordering = ['name']

    def get_serializer_class(self):
        if self.action == 'retrieve':
            return ProtocolDetailSerializer
        return ProtocolSerializer

    def perform_create(self, serializer):
        with reversion.create_revision():
            instance = serializer.save(created_by=self.request.user if self.request.user.is_authenticated else None)
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment("Created via API")
            return instance

    @action(detail=True, methods=['get'])
    def assays(self, request, pk=None):
        """List assays using this protocol."""
        protocol = self.get_object()
        assays = protocol.assays.all()
        serializer = AssayListSerializer(assays, many=True)
        return Response(serializer.data)


class AssayViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Assays."""
    queryset = Assay.objects.select_related('protocol', 'target', 'created_by')
    filter_backends = [DjangoFilterBackend, filters.OrderingFilter]
    filterset_fields = ['protocol', 'target']
    ordering_fields = ['created_at']
    ordering = ['-created_at']

    def get_serializer_class(self):
        if self.action == 'create':
            return AssayCreateSerializer
        if self.action == 'retrieve':
            return AssayDetailSerializer
        return AssayListSerializer

    @action(detail=True, methods=['get'])
    def data_series(self, request, pk=None):
        """List all data series for this assay."""
        assay = self.get_object()
        series = assay.data_series.select_related('compound', 'analysis')
        serializer = DataSeriesListSerializer(series, many=True)
        return Response(serializer.data)

    @action(detail=True, methods=['post'])
    def discover_series(self, request, pk=None):
        """
        Parse the data file and discover data series.

        This is a placeholder - actual implementation would parse
        the Excel/CSV file and create DataSeries objects.
        """
        assay = self.get_object()
        # TODO: Implement actual parsing logic from legacy code
        return Response({
            'status': 'not_implemented',
            'message': 'Data series discovery not yet implemented'
        })

    @action(detail=True, methods=['post'])
    def analyse_all(self, request, pk=None):
        """
        Run curve fitting analysis on all data series in this assay.

        Returns summary of analysis results.
        """
        import logging
        logger = logging.getLogger(__name__)

        try:
            assay = self.get_object()
            logger.info(f"Starting analysis for assay {assay.id} with {assay.data_series.count()} data series")
            results = analyse_assay(assay)
            logger.info(f"Analysis completed: {results['successful']} successful, {results['failed']} failed")
            return Response({
                'status': 'completed',
                'total': results['total'],
                'successful': results['successful'],
                'failed': results['failed'],
                'valid': results['valid'],
                'invalid': results['invalid'],
            })
        except Exception as e:
            logger.exception(f"Analysis failed for assay {pk}: {e}")
            return Response({
                'status': 'error',
                'error': str(e),
            }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


class DataSeriesViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Data Series."""
    queryset = DataSeries.objects.select_related('assay', 'compound', 'analysis', 'dilution_series')
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['assay', 'compound']
    search_fields = ['compound_name']
    ordering_fields = ['compound_name']
    ordering = ['compound_name']

    def get_serializer_class(self):
        if self.action == 'create':
            return DataSeriesCreateSerializer
        if self.action == 'retrieve':
            return DataSeriesDetailSerializer
        return DataSeriesListSerializer

    @action(detail=True, methods=['post'])
    def analyse(self, request, pk=None):
        """
        Run curve fitting analysis on this data series.

        Returns the analysis result.
        """
        series = self.get_object()
        fitting_method = series.assay.protocol.get_effective_fitting_method()

        try:
            analysis = analyse_data_series(series, fitting_method)
            return Response({
                'status': analysis.status,
                'results': analysis.results,
                'kpi_value': analysis.kpi_value,
            })
        except Exception as e:
            return Response({
                'status': 'error',
                'error': str(e),
            }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


class AnalysisResultViewSet(viewsets.ReadOnlyModelViewSet):
    """Read-only operations for Analysis Results."""
    queryset = AnalysisResult.objects.all()
    serializer_class = AnalysisResultSerializer
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ['status']


class HypothesisViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Hypotheses."""
    queryset = Hypothesis.objects.select_related('target', 'parent', 'product_compound')
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['target', 'status']
    search_fields = ['smiles', 'rationale']
    ordering_fields = ['created_at', 'updated_at', 'status']
    ordering = ['-created_at']

    def get_serializer_class(self):
        if self.action == 'retrieve':
            return HypothesisDetailSerializer
        return HypothesisListSerializer

    @action(detail=False, methods=['get'])
    def by_status(self, request):
        """Group hypotheses by status."""
        from django.db.models import Count
        status_counts = (
            Hypothesis.objects
            .values('status')
            .annotate(count=Count('id'))
            .order_by('status')
        )
        return Response(list(status_counts))
