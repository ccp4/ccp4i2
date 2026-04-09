"""
Constructs ViewSets

DRF ViewSets for construct database models with reversion support.
"""

import logging

import reversion
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import viewsets, filters, status
from rest_framework.decorators import action
from rest_framework.response import Response
from reversion.models import Version

from users.permissions import IsContributorOrReadOnly

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
from .serializers import (
    ConstructProjectSerializer,
    ConstructProjectDetailSerializer,
    ConstructProjectCreateSerializer,
    PlasmidListSerializer,
    PlasmidDetailSerializer,
    PlasmidCreateSerializer,
    ProteinListSerializer,
    ProteinDetailSerializer,
    ProteinCreateSerializer,
    ProteinSynonymSerializer,
    ProteinUseSerializer,
    CassetteListSerializer,
    CassetteDetailSerializer,
    CassetteCreateSerializer,
    CassetteUseListSerializer,
    CassetteUseDetailSerializer,
    CassetteUseCreateSerializer,
    SequencingResultSerializer,
    SequencingResultCreateSerializer,
    ExpressionTagTypeSerializer,
    ProteaseSerializer,
    ExpressionTagLocationSerializer,
    ExpressionTagSerializer,
)

logger = logging.getLogger(__name__)


class ReversionMixin:
    """
    Mixin to add reversion support and audit field tracking to ViewSets.

    Automatically:
    - Creates revision entries for create/update/delete operations
    - Sets created_by on new records
    - Sets modified_by on updates
    - Generates descriptive revision comments including model name and identifier
    """

    def _get_instance_identifier(self, instance):
        """Get a human-readable identifier for the instance."""
        # Try common identifier patterns
        if hasattr(instance, 'name'):
            return instance.name
        if hasattr(instance, 'plasmid_name'):
            return instance.plasmid_name
        return str(instance.pk)[:8]

    def _get_model_name(self, instance):
        """Get the model name for the instance."""
        return instance._meta.verbose_name

    def _set_audit_user(self, instance, field_name):
        """Set an audit user field if it exists and user is authenticated."""
        if hasattr(instance, field_name) and self.request.user.is_authenticated:
            setattr(instance, field_name, self.request.user)
            return True
        return False

    def perform_create(self, serializer):
        with reversion.create_revision():
            instance = serializer.save()

            # Set audit fields if not already set by serializer
            user_set = False
            if hasattr(instance, 'created_by') and getattr(instance, 'created_by') is None:
                if self._set_audit_user(instance, 'created_by'):
                    user_set = True

            if user_set:
                instance.save(update_fields=['created_by'])

            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)

            model_name = self._get_model_name(instance)
            identifier = self._get_instance_identifier(instance)
            reversion.set_comment(f"Created {model_name}: {identifier}")
            return instance

    def perform_update(self, serializer):
        with reversion.create_revision():
            instance = serializer.instance

            # Set modified_by before save
            if hasattr(instance, 'modified_by') and self.request.user.is_authenticated:
                serializer.validated_data['modified_by'] = self.request.user

            instance = serializer.save()

            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)

            model_name = self._get_model_name(instance)
            identifier = self._get_instance_identifier(instance)
            reversion.set_comment(f"Updated {model_name}: {identifier}")
            return instance

    def perform_destroy(self, instance):
        model_name = self._get_model_name(instance)
        identifier = self._get_instance_identifier(instance)

        with reversion.create_revision():
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment(f"Deleted {model_name}: {identifier}")
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


# =============================================================================
# Reference Data ViewSets
# =============================================================================

class ExpressionTagTypeViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Expression Tag Types."""
    queryset = ExpressionTagType.objects.all()
    serializer_class = ExpressionTagTypeSerializer
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [filters.SearchFilter, filters.OrderingFilter]
    search_fields = ['name']
    ordering_fields = ['name']
    ordering = ['name']


class ProteaseViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Proteases."""
    queryset = Protease.objects.all()
    serializer_class = ProteaseSerializer
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [filters.SearchFilter, filters.OrderingFilter]
    search_fields = ['name']
    ordering_fields = ['name']
    ordering = ['name']


class ExpressionTagLocationViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Expression Tag Locations."""
    queryset = ExpressionTagLocation.objects.all()
    serializer_class = ExpressionTagLocationSerializer
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [filters.SearchFilter, filters.OrderingFilter]
    search_fields = ['name']
    ordering_fields = ['name']
    ordering = ['name']


class ExpressionTagViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Expression Tags."""
    queryset = ExpressionTag.objects.select_related(
        'expression_tag_type', 'protease', 'cassette_use', 'location'
    )
    serializer_class = ExpressionTagSerializer
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.OrderingFilter]
    filterset_fields = ['cassette_use', 'expression_tag_type', 'location']
    ordering = ['cassette_use__plasmid__ncn_id']


# =============================================================================
# Construct Project ViewSet
# =============================================================================

class ConstructProjectViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Construct Projects."""
    queryset = ConstructProject.objects.select_related('parent', 'created_by')
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['parent']
    search_fields = ['name']
    ordering_fields = ['name', 'created_at']
    ordering = ['name']

    def get_serializer_class(self):
        if self.action == 'create':
            return ConstructProjectCreateSerializer
        if self.action == 'retrieve':
            return ConstructProjectDetailSerializer
        return ConstructProjectSerializer

    @action(detail=True, methods=['get'])
    def plasmids(self, request, pk=None):
        """List plasmids for this project."""
        project = self.get_object()
        plasmids = project.plasmids.all()
        serializer = PlasmidListSerializer(plasmids, many=True)
        return Response(serializer.data)


# =============================================================================
# Protein ViewSets
# =============================================================================

class ProteinViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Proteins."""
    queryset = Protein.objects.prefetch_related('synonyms', 'cassettes')
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [filters.SearchFilter, filters.OrderingFilter]
    search_fields = ['uniprot_id', 'synonyms__name']
    ordering_fields = ['uniprot_id', 'created_at']
    ordering = ['uniprot_id']

    def get_serializer_class(self):
        if self.action == 'create':
            return ProteinCreateSerializer
        if self.action == 'retrieve':
            return ProteinDetailSerializer
        return ProteinListSerializer

    @action(detail=True, methods=['get'])
    def cassettes(self, request, pk=None):
        """List cassettes for this protein."""
        protein = self.get_object()
        cassettes = protein.cassettes.all()
        serializer = CassetteListSerializer(cassettes, many=True)
        return Response(serializer.data)

    @action(detail=False, methods=['get'])
    def by_uniprot_id(self, request):
        """Lookup protein by UniProt ID."""
        uniprot_id = request.query_params.get('uniprot_id')
        if not uniprot_id:
            return Response(
                {'error': 'uniprot_id parameter required'},
                status=status.HTTP_400_BAD_REQUEST
            )
        try:
            protein = Protein.objects.get(uniprot_id=uniprot_id)
            serializer = ProteinDetailSerializer(protein)
            return Response(serializer.data)
        except Protein.DoesNotExist:
            return Response(
                {'error': 'Protein not found'},
                status=status.HTTP_404_NOT_FOUND
            )


class ProteinSynonymViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Protein Synonyms."""
    queryset = ProteinSynonym.objects.select_related('protein')
    serializer_class = ProteinSynonymSerializer
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['protein']
    search_fields = ['name']
    ordering_fields = ['name']
    ordering = ['name']


class ProteinUseViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Protein Uses (protein-project links)."""
    queryset = ProteinUse.objects.select_related('protein', 'project')
    serializer_class = ProteinUseSerializer
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.OrderingFilter]
    filterset_fields = ['protein', 'project']
    ordering = ['protein__uniprot_id']


# =============================================================================
# Cassette ViewSets
# =============================================================================

class CassetteViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Cassettes."""
    queryset = Cassette.objects.select_related('protein', 'created_by')
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['protein']
    search_fields = ['protein__uniprot_id']
    ordering_fields = ['protein__uniprot_id', 'start', 'created_at']
    ordering = ['protein__uniprot_id', 'start']

    def get_serializer_class(self):
        if self.action == 'create':
            return CassetteCreateSerializer
        if self.action == 'retrieve':
            return CassetteDetailSerializer
        return CassetteListSerializer

    @action(detail=True, methods=['get'])
    def plasmids(self, request, pk=None):
        """List plasmids using this cassette."""
        cassette = self.get_object()
        cassette_uses = cassette.plasmid_uses.select_related('plasmid')
        plasmids = [cu.plasmid for cu in cassette_uses]
        serializer = PlasmidListSerializer(plasmids, many=True)
        return Response(serializer.data)


class CassetteUseViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Cassette Uses."""
    queryset = CassetteUse.objects.select_related(
        'cassette', 'cassette__protein', 'plasmid', 'created_by'
    ).prefetch_related('expression_tags')
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.OrderingFilter]
    filterset_fields = ['cassette', 'plasmid']
    ordering = ['plasmid__ncn_id']

    def get_serializer_class(self):
        if self.action == 'create':
            return CassetteUseCreateSerializer
        if self.action == 'retrieve':
            return CassetteUseDetailSerializer
        return CassetteUseListSerializer

    @action(detail=True, methods=['get'])
    def expression_tags(self, request, pk=None):
        """List expression tags for this cassette use."""
        cassette_use = self.get_object()
        tags = cassette_use.expression_tags.select_related(
            'expression_tag_type', 'protease', 'location'
        )
        serializer = ExpressionTagSerializer(tags, many=True)
        return Response(serializer.data)


# =============================================================================
# Sequencing Result ViewSet
# =============================================================================

class SequencingResultViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Sequencing Results."""
    queryset = SequencingResult.objects.select_related(
        'cassette_use', 'cassette_use__cassette', 'plasmid'
    )
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.OrderingFilter]
    filterset_fields = ['cassette_use', 'plasmid']
    ordering = ['-created_at']

    def get_serializer_class(self):
        if self.action == 'create':
            return SequencingResultCreateSerializer
        return SequencingResultSerializer


# =============================================================================
# Plasmid ViewSet
# =============================================================================

class PlasmidViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Plasmids."""
    queryset = Plasmid.objects.select_related('project', 'parent', 'created_by')
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['project', 'parent']
    search_fields = ['ncn_id', 'name', 'formatted_id']
    ordering_fields = ['ncn_id', 'created_at', 'name']
    ordering = ['-ncn_id']

    def get_serializer_class(self):
        if self.action == 'create':
            return PlasmidCreateSerializer
        if self.action == 'retrieve':
            return PlasmidDetailSerializer
        return PlasmidListSerializer

    @action(detail=False, methods=['get'])
    def by_ncn_id(self, request):
        """Lookup plasmid by NCN ID number."""
        ncn_id = request.query_params.get('ncn_id')
        if not ncn_id:
            return Response(
                {'error': 'ncn_id parameter required'},
                status=status.HTTP_400_BAD_REQUEST
            )
        try:
            plasmid = Plasmid.objects.get(ncn_id=int(ncn_id))
            serializer = PlasmidDetailSerializer(
                plasmid, context={'request': request}
            )
            return Response(serializer.data)
        except (ValueError, Plasmid.DoesNotExist):
            return Response(
                {'error': 'Plasmid not found'},
                status=status.HTTP_404_NOT_FOUND
            )

    @action(detail=True, methods=['get'])
    def genbank_content(self, request, pk=None):
        """
        Return GenBank file content for SeqViz visualization.

        Returns the raw content of the plasmid's GenBank file as text.
        """
        plasmid = self.get_object()
        if not plasmid.genbank_file:
            return Response(
                {'error': 'No GenBank file attached'},
                status=status.HTTP_404_NOT_FOUND
            )
        try:
            content = plasmid.genbank_file.read().decode('utf-8')
            plasmid.genbank_file.seek(0)  # Reset file pointer
            return Response({'content': content})
        except Exception as e:
            logger.error(f"Error reading GenBank file for plasmid {plasmid.id}: {e}")
            return Response(
                {'error': 'Failed to read GenBank file'},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR
            )

    @action(detail=True, methods=['get'])
    def cassette_uses(self, request, pk=None):
        """List cassette uses for this plasmid."""
        plasmid = self.get_object()
        cassette_uses = plasmid.cassette_uses.select_related(
            'cassette', 'cassette__protein'
        ).prefetch_related('expression_tags')
        serializer = CassetteUseDetailSerializer(cassette_uses, many=True)
        return Response(serializer.data)

    @action(detail=True, methods=['get'])
    def sequencing_results(self, request, pk=None):
        """List sequencing results for this plasmid."""
        plasmid = self.get_object()
        results = plasmid.sequencing_results.select_related(
            'cassette_use', 'cassette_use__cassette'
        )
        serializer = SequencingResultSerializer(results, many=True)
        return Response(serializer.data)

    @action(detail=True, methods=['get'])
    def children(self, request, pk=None):
        """List plasmids derived from this plasmid."""
        plasmid = self.get_object()
        children = plasmid.children.all()
        serializer = PlasmidListSerializer(children, many=True)
        return Response(serializer.data)
