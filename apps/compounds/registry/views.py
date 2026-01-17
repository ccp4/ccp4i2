"""
Registry ViewSets

DRF ViewSets for compound registration models with reversion support.
"""

import logging

import reversion
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import viewsets, filters, status
from rest_framework.decorators import action
from rest_framework.response import Response
from reversion.models import Version

from .models import Supplier, Target, Compound, Batch, BatchQCFile, CompoundTemplate

logger = logging.getLogger(__name__)

# Try to import RDKit for structure searches
try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not available - structure search will be disabled")
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

    @action(detail=False, methods=['get', 'post'])
    def my_supplier(self, request):
        """
        Get or create a supplier linked to the current user.

        GET: Returns the current user's linked supplier, or null if none exists.
             Also returns suggested_name for creating a new one.

        POST: Creates a supplier linked to the current user.
              Optional body: {"name": "Custom Name", "initials": "CN"}
              If not provided, uses user's full name or email.

        Only works when authentication is enabled.
        """
        if not request.user.is_authenticated:
            return Response(
                {'error': 'Authentication required'},
                status=status.HTTP_401_UNAUTHORIZED
            )

        # Build suggested name from user info
        if request.user.first_name and request.user.last_name:
            suggested_name = f"{request.user.first_name} {request.user.last_name}"
        elif request.user.first_name:
            suggested_name = request.user.first_name
        elif request.user.email:
            # Use email prefix as name
            suggested_name = request.user.email.split('@')[0].replace('.', ' ').replace('_', ' ').title()
        else:
            suggested_name = f"User {request.user.id}"

        # Generate initials from suggested name
        name_parts = suggested_name.split()
        if len(name_parts) >= 2:
            suggested_initials = (name_parts[0][0] + name_parts[-1][0]).upper()
        else:
            suggested_initials = suggested_name[:2].upper()

        if request.method == 'GET':
            # Try to find existing supplier linked to this user
            try:
                supplier = Supplier.objects.get(user=request.user)
                serializer = SupplierSerializer(supplier, context={'request': request})
                return Response({
                    'supplier': serializer.data,
                    'suggested_name': suggested_name,
                    'suggested_initials': suggested_initials,
                })
            except Supplier.DoesNotExist:
                return Response({
                    'supplier': None,
                    'suggested_name': suggested_name,
                    'suggested_initials': suggested_initials,
                })

        else:  # POST - create supplier
            # Check if user already has a supplier
            if Supplier.objects.filter(user=request.user).exists():
                return Response(
                    {'error': 'You already have a linked supplier'},
                    status=status.HTTP_400_BAD_REQUEST
                )

            # Get name from request or use suggested
            name = request.data.get('name', suggested_name)
            initials = request.data.get('initials', suggested_initials)

            # Check for name uniqueness
            if Supplier.objects.filter(name=name).exists():
                return Response(
                    {'error': f'A supplier named "{name}" already exists'},
                    status=status.HTTP_400_BAD_REQUEST
                )

            # Check for initials uniqueness (if provided)
            if initials and Supplier.objects.filter(initials=initials).exists():
                return Response(
                    {'error': f'A supplier with initials "{initials}" already exists'},
                    status=status.HTTP_400_BAD_REQUEST
                )

            # Create the supplier
            with reversion.create_revision():
                supplier = Supplier.objects.create(
                    name=name,
                    initials=initials if initials else None,
                    user=request.user
                )
                reversion.set_user(request.user)
                reversion.set_comment("Created personal supplier via API")

            serializer = SupplierSerializer(supplier, context={'request': request})
            return Response({
                'supplier': serializer.data,
                'suggested_name': suggested_name,
                'suggested_initials': suggested_initials,
            }, status=status.HTTP_201_CREATED)


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
    # Use dict format to enable __in lookup for reg_number (batch SMILES lookup)
    filterset_fields = {
        'target': ['exact'],
        'supplier': ['exact'],
        'stereo_comment': ['exact'],
        'reg_number': ['exact', 'in'],
    }
    search_fields = ['reg_number', 'smiles', 'supplier_ref', 'comments']
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

    @action(detail=False, methods=['post'])
    def bulk_create(self, request):
        """
        Create multiple compounds in a single request.

        Request body:
        {
            "compounds": [
                {
                    "target": "<uuid>",
                    "smiles": "...",
                    "stereo_comment": "...",
                    "supplier": "<uuid>",
                    "supplier_ref": "...",
                    "comments": "...",
                    "labbook_number": 123,
                    "page_number": 45
                },
                ...
            ]
        }

        Returns:
        {
            "created": [{"id": "...", "formatted_id": "NCL-..."}],
            "errors": [{"index": 0, "error": "..."}]
        }
        """
        compounds_data = request.data.get('compounds', [])

        if not compounds_data:
            return Response(
                {'error': 'No compounds provided'},
                status=status.HTTP_400_BAD_REQUEST
            )

        if len(compounds_data) > 1000:
            return Response(
                {'error': 'Maximum 1000 compounds per request'},
                status=status.HTTP_400_BAD_REQUEST
            )

        created = []
        errors = []

        for idx, compound_data in enumerate(compounds_data):
            try:
                serializer = CompoundCreateSerializer(
                    data=compound_data,
                    context={'request': request}
                )
                if serializer.is_valid():
                    with reversion.create_revision():
                        instance = serializer.save()
                        if request.user.is_authenticated:
                            reversion.set_user(request.user)
                        reversion.set_comment("Created via bulk import")

                    created.append({
                        'id': str(instance.id),
                        'formatted_id': instance.formatted_id,
                        'index': idx,
                    })
                else:
                    errors.append({
                        'index': idx,
                        'error': serializer.errors,
                    })
            except Exception as e:
                errors.append({
                    'index': idx,
                    'error': str(e),
                })

        return Response({
            'created': created,
            'errors': errors,
            'total_created': len(created),
            'total_errors': len(errors),
        })

    @action(detail=False, methods=['get'])
    def structure_search(self, request):
        """
        Search compounds by structure (substructure or superstructure).

        Query parameters:
            smiles: SMILES or SMARTS pattern to search for
            mode: 'substructure' or 'superstructure' (default: 'substructure')
            target: Optional target UUID to filter by

        Substructure search: Find compounds containing the query as a substructure
                            (query is a fragment of the compound)
        Superstructure search: Find compounds that are substructures of the query
                              (compound is a fragment of the query)
        """
        if not RDKIT_AVAILABLE:
            return Response(
                {'error': 'Structure search not available - RDKit not installed'},
                status=status.HTTP_503_SERVICE_UNAVAILABLE
            )

        smiles = request.query_params.get('smiles')
        mode = request.query_params.get('mode', 'substructure')
        target_id = request.query_params.get('target')

        if not smiles:
            return Response(
                {'error': 'smiles parameter required'},
                status=status.HTTP_400_BAD_REQUEST
            )

        if mode not in ('substructure', 'superstructure'):
            return Response(
                {'error': 'mode must be "substructure" or "superstructure"'},
                status=status.HTTP_400_BAD_REQUEST
            )

        # Parse the query pattern
        try:
            # Try as SMARTS first (for query patterns with wildcards)
            query_mol = Chem.MolFromSmarts(smiles)
            if query_mol is None:
                # Try as SMILES
                query_mol = Chem.MolFromSmiles(smiles)
            if query_mol is None:
                return Response(
                    {'error': 'Invalid SMILES/SMARTS pattern'},
                    status=status.HTTP_400_BAD_REQUEST
                )
        except Exception as e:
            return Response(
                {'error': f'Failed to parse structure: {str(e)}'},
                status=status.HTTP_400_BAD_REQUEST
            )

        # Get compounds to search through
        queryset = self.get_queryset()
        if target_id:
            queryset = queryset.filter(target_id=target_id)

        # Perform structure search
        matches = []
        for compound in queryset:
            if not compound.smiles:
                continue

            try:
                compound_mol = Chem.MolFromSmiles(compound.smiles)
                if compound_mol is None:
                    continue

                if mode == 'substructure':
                    # Find compounds containing the query as a substructure
                    # i.e., query is a fragment of the compound
                    if compound_mol.HasSubstructMatch(query_mol):
                        matches.append(compound)
                else:  # superstructure
                    # Find compounds that are substructures of the query
                    # i.e., compound is a fragment of the query
                    if query_mol.HasSubstructMatch(compound_mol):
                        matches.append(compound)

            except Exception as e:
                logger.warning(f"Error processing compound {compound.id}: {e}")
                continue

        serializer = CompoundListSerializer(matches, many=True)
        return Response({
            'query': smiles,
            'mode': mode,
            'count': len(matches),
            'matches': serializer.data
        })


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
