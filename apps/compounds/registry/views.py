"""
Registry ViewSets

DRF ViewSets for compound registration models with reversion support.

Permission model:
- Read access: Anyone (including 'user' operating level)
- Create access: Requires 'contributor' or 'admin' operating level
- Update/Delete access: Compounds require 'admin' only; other models allow 'contributor'
"""

import logging
import re

import reversion
from django.db.models import Count, Max
from django.db.models.functions import Coalesce, Greatest
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import viewsets, filters, status
from rest_framework.decorators import action
from rest_framework.parsers import MultiPartParser
from rest_framework.response import Response
from reversion.models import Version

from compounds.formatting import format_compound_id, get_compound_pattern
from users.permissions import IsContributorOrReadOnly, IsContributorCreateAdminUpdate, can_administer
from .filters import CompoundSearchFilter
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
    TargetDashboardSerializer,
    DashboardProjectSerializer,
    SavedAggregationViewSerializer,
    CompoundListSerializer,
    CompoundDetailSerializer,
    CompoundCreateSerializer,
    BatchSerializer,
    BatchQCFileSerializer,
    BatchQCFileCreateSerializer,
    CompoundTemplateSerializer,
)


class ReversionMixin:
    """
    Mixin to add reversion support and audit field tracking to ViewSets.

    Automatically:
    - Creates revision entries for create/update/delete operations
    - Sets created_by/registered_by on new records
    - Sets modified_by on updates
    - Generates descriptive revision comments including model name and identifier
    """

    def _get_instance_identifier(self, instance):
        """Get a human-readable identifier for the instance."""
        # Try common identifier patterns
        if hasattr(instance, 'formatted_id'):
            return instance.formatted_id
        if hasattr(instance, 'name'):
            return instance.name
        if hasattr(instance, 'batch_number') and hasattr(instance, 'compound'):
            return f"{instance.compound.formatted_id}/{instance.batch_number}"
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
            for field in ['created_by', 'registered_by']:
                if hasattr(instance, field) and getattr(instance, field) is None:
                    if self._set_audit_user(instance, field):
                        user_set = True

            if user_set:
                instance.save(update_fields=[f for f in ['created_by', 'registered_by']
                                             if hasattr(instance, f)])

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
        } for v in versions[:50]])  # Limit to 50 most recent


class SupplierViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Suppliers."""
    queryset = Supplier.objects.all()
    serializer_class = SupplierSerializer
    permission_classes = [IsContributorOrReadOnly]
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
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['parent']
    search_fields = ['name']
    ordering_fields = ['name', 'created_at', 'latest_activity']
    ordering = ['-latest_activity', 'name']  # Most recent activity first, then alphabetically

    def get_queryset(self):
        """
        Annotate targets with latest_activity date for sorting.

        This computes the most recent date from either:
        - The most recently registered compound
        - The most recently created assay
        """
        return Target.objects.annotate(
            latest_compound=Max('compounds__registered_at'),
            latest_assay=Max('assays__created_at'),
        ).annotate(
            # Use Greatest to get the most recent of compound or assay dates
            # Coalesce handles nulls - falls back to created_at if no activity
            latest_activity=Coalesce(
                Greatest('latest_compound', 'latest_assay'),
                'latest_compound',
                'latest_assay',
                'created_at',
            )
        )

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

    @action(detail=True, methods=['get'])
    def dashboard(self, request, pk=None):
        """
        Get dashboard data for a target.

        Returns target info with:
        - recent_compounds: 50 most recently registered compounds
        - recent_assays: 50 most recent assays for this target
        """
        target = self.get_object()

        # Get 50 most recent compounds
        recent_compounds = target.compounds.order_by('-registered_at')[:50]

        # Get 50 most recent assays for this target
        from compounds.assays.models import Assay
        recent_assays_qs = (
            Assay.objects.filter(target=target)
            .select_related('protocol')
            .annotate(data_series_count=Count('data_series'))
            .order_by('-created_at')[:50]
        )

        # Transform assays to dashboard format
        recent_assays = [
            {
                'id': assay.id,
                'protocol_name': assay.protocol.name,
                'created_at': assay.created_at,
                'data_series_count': assay.data_series_count,
            }
            for assay in recent_assays_qs
        ]

        # Build response with dashboard data
        serializer = TargetDashboardSerializer(target)
        data = serializer.data
        data['recent_compounds'] = CompoundListSerializer(recent_compounds, many=True).data
        data['recent_assays'] = recent_assays

        return Response(data)

    @action(detail=True, methods=['post'], parser_classes=[MultiPartParser])
    def upload_image(self, request, pk=None):
        """
        Upload a branding image for the target dashboard.

        Expects multipart form data with 'image' field.
        """
        target = self.get_object()

        if 'image' not in request.FILES:
            return Response(
                {'error': 'No image file provided'},
                status=status.HTTP_400_BAD_REQUEST
            )

        image_file = request.FILES['image']

        # Validate file type
        allowed_types = ['image/jpeg', 'image/png', 'image/gif', 'image/webp']
        if image_file.content_type not in allowed_types:
            return Response(
                {'error': f'Invalid file type. Allowed: {", ".join(allowed_types)}'},
                status=status.HTTP_400_BAD_REQUEST
            )

        # Delete old image if exists
        if target.image:
            target.image.delete(save=False)

        # Save new image
        with reversion.create_revision():
            target.image = image_file
            target.save()
            if request.user.is_authenticated:
                reversion.set_user(request.user)
            reversion.set_comment("Uploaded dashboard image")

        serializer = TargetSerializer(target, context={'request': request})
        return Response(serializer.data)

    @action(detail=True, methods=['delete'])
    def delete_image(self, request, pk=None):
        """Delete the target's branding image."""
        target = self.get_object()

        if not target.image:
            return Response(
                {'error': 'Target has no image'},
                status=status.HTTP_404_NOT_FOUND
            )

        with reversion.create_revision():
            target.image.delete(save=True)
            if request.user.is_authenticated:
                reversion.set_user(request.user)
            reversion.set_comment("Deleted dashboard image")

        return Response({'success': True})

    @action(detail=True, methods=['get'])
    def recent_projects(self, request, pk=None):
        """
        Get CCP4i2 projects that match compounds in this target.

        Projects are matched by name containing compound IDs in patterns like:
        - ncl-XXXXX or NCL-XXXXXXXX (case-insensitive)
        - The integer portion must match a compound's reg_number

        Returns projects ordered by most recent activity.
        """
        target = self.get_object()

        # Get all compound reg_numbers for this target
        reg_numbers = list(
            target.compounds.values_list('reg_number', flat=True)
        )

        if not reg_numbers:
            return Response([])

        # Build regex pattern to match compound ID formats
        # Pattern matches: PREFIX-12345, prefix-00012345, etc.
        ncl_pattern = get_compound_pattern(capturing=True)

        try:
            from ccp4i2.db.models import Project
        except ImportError:
            logger.warning("CCP4i2 models not available - cannot query projects")
            return Response([])

        # Query all projects and filter by compound ID match
        reg_number_set = set(reg_numbers)
        matching_projects = []

        for project in Project.objects.annotate(job_count=Count('jobs')).order_by('-last_access'):
            # Find all NCL-XXXXX patterns in the project name
            matches = ncl_pattern.findall(project.name)
            matched_ids = []

            for match in matches:
                try:
                    reg_num = int(match)
                    if reg_num in reg_number_set:
                        matched_ids.append(format_compound_id(reg_num))
                except ValueError:
                    continue

            if matched_ids:
                matching_projects.append({
                    'id': project.id,
                    'name': project.name,
                    'last_access': project.last_access,
                    'job_count': project.job_count,
                    'matching_compound_ids': matched_ids,
                })

            # Limit to 50 projects
            if len(matching_projects) >= 50:
                break

        serializer = DashboardProjectSerializer(matching_projects, many=True)
        return Response(serializer.data)

    @action(detail=True, methods=['get', 'post', 'delete'])
    def saved_view(self, request, pk=None):
        """
        Get, set, or delete the saved aggregation view for this target.

        GET: Returns the current saved view configuration (or null)
        POST: Saves a new view configuration (admin only)
        DELETE: Removes the saved view (admin only)

        Request body for POST:
        {
            "protocol_names": ["Protocol A"],
            "compound_search": "",
            "output_format": "compact",
            "aggregations": ["geomean", "count"],
            "status": "valid"
        }
        """
        target = self.get_object()

        if request.method == 'GET':
            return Response({
                'saved_view': target.saved_aggregation_view
            })

        # Write operations require admin
        if not can_administer(request):
            return Response(
                {'error': 'Admin access required'},
                status=status.HTTP_403_FORBIDDEN
            )

        if request.method == 'DELETE':
            with reversion.create_revision():
                target.saved_aggregation_view = None
                target.save(update_fields=['saved_aggregation_view'])
                if request.user.is_authenticated:
                    reversion.set_user(request.user)
                reversion.set_comment("Removed saved aggregation view")
            return Response({'success': True})

        # POST - validate and save
        serializer = SavedAggregationViewSerializer(data=request.data)
        if not serializer.is_valid():
            return Response(
                {'error': serializer.errors},
                status=status.HTTP_400_BAD_REQUEST
            )

        with reversion.create_revision():
            target.saved_aggregation_view = serializer.validated_data
            target.save(update_fields=['saved_aggregation_view'])
            if request.user.is_authenticated:
                reversion.set_user(request.user)
            reversion.set_comment("Updated saved aggregation view")

        return Response({
            'success': True,
            'saved_view': target.saved_aggregation_view
        })


class CompoundViewSet(ReversionMixin, viewsets.ModelViewSet):
    """
    CRUD operations for Compounds.

    Permission model:
    - READ: Anyone (including 'user' operating level)
    - CREATE: Requires 'contributor' or 'admin' operating level
    - UPDATE/DELETE: Requires 'admin' operating level only
    """
    queryset = Compound.objects.select_related('target', 'supplier', 'registered_by')
    permission_classes = [IsContributorCreateAdminUpdate]
    filter_backends = [DjangoFilterBackend, CompoundSearchFilter, filters.OrderingFilter]
    # Use dict format to enable __in lookup for reg_number (batch SMILES lookup)
    filterset_fields = {
        'target': ['exact'],
        'supplier': ['exact'],
        'stereo_comment': ['exact'],
        'reg_number': ['exact', 'in'],
    }
    # Note: search_fields is handled by CompoundSearchFilter which properly parses
    # NCL format identifiers (NCL-00018710, NCL-000187, etc.) and searches reg_number
    # by prefix matching, plus standard icontains on text fields
    search_fields = ['reg_number', 'smiles', 'supplier_ref', 'comments']
    ordering_fields = ['reg_number', 'registered_at', 'molecular_weight']
    ordering = ['-reg_number']

    def get_queryset(self):
        """Optimize queryset with batch count annotation for list view."""
        queryset = super().get_queryset()
        if self.action == 'list':
            # Annotate batch_count to avoid N+1 queries in list view
            return queryset.annotate(batch_count=Count('batches'))
        return queryset

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

    @action(detail=True, methods=['get'])
    def adjacent(self, request, pk=None):
        """
        Get adjacent compounds (previous and next) within the same target.

        Returns IDs and formatted_ids of the previous and next compounds,
        ordered by registration number (reg_number).

        Response:
        {
            "previous": {"id": "...", "formatted_id": "NCL-..."} | null,
            "next": {"id": "...", "formatted_id": "NCL-..."} | null
        }
        """
        compound = self.get_object()

        def format_result(result):
            """Convert query result to response format with formatted_id."""
            if result is None:
                return None
            return {
                'id': str(result['id']),
                'formatted_id': format_compound_id(result['reg_number']),
            }

        # Get previous compound (lower reg_number, same target)
        previous = (
            Compound.objects
            .filter(target=compound.target, reg_number__lt=compound.reg_number)
            .order_by('-reg_number')
            .values('id', 'reg_number')
            .first()
        )

        # Get next compound (higher reg_number, same target)
        next_compound = (
            Compound.objects
            .filter(target=compound.target, reg_number__gt=compound.reg_number)
            .order_by('reg_number')
            .values('id', 'reg_number')
            .first()
        )

        return Response({
            'previous': format_result(previous),
            'next': format_result(next_compound),
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

    @action(detail=False, methods=['get'])
    def resolve_by_smiles(self, request):
        """
        Resolve a compound by SMILES via InChI matching.

        Converts the input SMILES to InChI server-side and performs an exact
        lookup against the InChI field. This finds the canonical registry
        compound regardless of SMILES representation differences.

        Query parameters:
            smiles: SMILES string to resolve (required)

        Returns:
            JSON with found=true and compound data if matched,
            or found=false if no match.
        """
        if not RDKIT_AVAILABLE:
            return Response(
                {'error': 'Compound resolution not available - RDKit not installed'},
                status=status.HTTP_503_SERVICE_UNAVAILABLE
            )

        smiles_param = request.query_params.get('smiles')
        if not smiles_param:
            return Response(
                {'error': 'smiles parameter required'},
                status=status.HTTP_400_BAD_REQUEST
            )

        try:
            mol = Chem.MolFromSmiles(smiles_param)
            if mol is None:
                return Response(
                    {'error': 'Invalid SMILES string'},
                    status=status.HTTP_400_BAD_REQUEST
                )
            from rdkit.Chem import inchi as rdkit_inchi
            query_inchi = rdkit_inchi.MolToInchi(mol)
        except Exception as e:
            return Response(
                {'error': f'Failed to process SMILES: {str(e)}'},
                status=status.HTTP_400_BAD_REQUEST
            )

        if not query_inchi:
            return Response(
                {'error': 'Could not generate InChI from SMILES'},
                status=status.HTTP_400_BAD_REQUEST
            )

        compound = self.get_queryset().filter(inchi=query_inchi).first()

        if compound:
            serializer = CompoundDetailSerializer(compound)
            return Response({
                'found': True,
                'query_smiles': smiles_param,
                'query_inchi': query_inchi,
                'compound': serializer.data,
            })
        else:
            return Response({
                'found': False,
                'query_smiles': smiles_param,
                'query_inchi': query_inchi,
                'compound': None,
            })


class BatchViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Batches."""
    queryset = Batch.objects.select_related('compound', 'supplier')
    serializer_class = BatchSerializer
    permission_classes = [IsContributorOrReadOnly]
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
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ['batch']

    def get_serializer_class(self):
        """Use create serializer for POST (file upload), read serializer otherwise."""
        if self.action == 'create':
            return BatchQCFileCreateSerializer
        return BatchQCFileSerializer


class CompoundTemplateViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Compound Templates."""
    queryset = CompoundTemplate.objects.select_related('target')
    serializer_class = CompoundTemplateSerializer
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.SearchFilter]
    filterset_fields = ['target']
    search_fields = ['name']
