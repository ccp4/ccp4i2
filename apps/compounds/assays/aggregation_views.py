"""
Aggregation Views

REST API endpoints for compound assay data aggregation.
"""

from rest_framework import viewsets, status
from rest_framework.decorators import action
from rest_framework.parsers import JSONParser
from rest_framework.permissions import AllowAny
from rest_framework.response import Response

from compounds.registry.models import Target
from .models import Protocol
from .aggregation import (
    build_data_series_queryset,
    aggregate_compact,
    aggregate_medium,
    aggregate_long,
)


class AggregationViewSet(viewsets.ViewSet):
    """
    ViewSet for data aggregation operations.

    Endpoints:
        POST /aggregations/aggregate/ - Run aggregation query
        GET /aggregations/protocols/ - List available protocols
        GET /aggregations/targets/ - List available targets
    """

    # Allow unauthenticated access and bypass CSRF - aggregation is read-only
    authentication_classes = []
    permission_classes = [AllowAny]
    parser_classes = [JSONParser]

    @action(detail=False, methods=['post'])
    def aggregate(self, request):
        """
        Execute an aggregation query.

        Request body:
        {
            "predicates": {
                "targets": ["<uuid>", ...],
                "compounds": ["<uuid>", ...],
                "compound_search": "NCL-00026",
                "protocols": ["<uuid>", ...],
                "status": "valid"  // optional, defaults to "valid"
            },
            "output_format": "compact" | "medium" | "long",
            "aggregations": ["geomean", "count", "stdev", "list"],
            "group_by_batch": false  // optional, when true splits results by batch
        }

        Returns:
            Compact format: One row per compound (or compound/batch) with protocol columns
            Medium format: One row per compound-protocol (or compound-batch-protocol) pair
            Long format: One row per measurement (always includes batch info)

        When group_by_batch is true:
            - Compact/Medium: Separate rows are created for each batch of a compound
            - Long: Batch info is always included; the flag affects sorting/display hints
            - Rows include batch_id and batch_number fields
            - meta.group_by_batch indicates the grouping mode used
        """
        predicates = request.data.get('predicates', {})
        output_format = request.data.get('output_format', 'compact')
        aggregations = request.data.get('aggregations', ['geomean', 'count'])
        group_by_batch = request.data.get('group_by_batch', False)

        # Validate output format
        if output_format not in ('compact', 'medium', 'long'):
            return Response(
                {'error': 'output_format must be "compact", "medium", or "long"'},
                status=status.HTTP_400_BAD_REQUEST
            )

        # Validate aggregations
        valid_aggregations = {'geomean', 'count', 'stdev', 'list'}
        invalid = set(aggregations) - valid_aggregations
        if invalid:
            return Response(
                {'error': f'Invalid aggregations: {invalid}. Valid options: {valid_aggregations}'},
                status=status.HTTP_400_BAD_REQUEST
            )

        # Build and execute query
        queryset = build_data_series_queryset(predicates)

        if output_format == 'compact':
            result = aggregate_compact(queryset, aggregations, group_by_batch=group_by_batch)
        elif output_format == 'medium':
            result = aggregate_medium(queryset, aggregations, group_by_batch=group_by_batch)
        else:
            result = aggregate_long(queryset, aggregations, group_by_batch=group_by_batch)

        return Response(result)

    @action(detail=False, methods=['get'])
    def protocols(self, request):
        """
        List available protocols for the predicate builder.

        Query params:
            target: Optional target UUID to filter protocols
            search: Optional text search

        Returns:
            List of {id, name, analysis_method}
        """
        queryset = Protocol.objects.all()

        # Filter by target if provided
        target_id = request.query_params.get('target')
        if target_id:
            # Find protocols that have assays for this target
            queryset = queryset.filter(assays__target_id=target_id).distinct()

        # Search by name
        search = request.query_params.get('search')
        if search:
            queryset = queryset.filter(name__icontains=search)

        protocols = queryset.values('id', 'name', 'analysis_method')[:100]

        return Response(list(protocols))

    @action(detail=False, methods=['get'])
    def targets(self, request):
        """
        List available targets for the predicate builder.

        Query params:
            search: Optional text search

        Returns:
            List of {id, name}
        """
        queryset = Target.objects.all()

        # Search by name
        search = request.query_params.get('search')
        if search:
            queryset = queryset.filter(name__icontains=search)

        targets = queryset.values('id', 'name')[:100]

        return Response(list(targets))

    @action(detail=False, methods=['get'])
    def data_series(self, request):
        """
        Fetch detailed data series for a compound-protocol pair.

        Query params:
            compound: Compound UUID (required)
            protocol: Protocol UUID (required)
            status: Analysis status filter (optional, default: 'valid')

        Returns:
            List of data series with full dose-response data for charting.
        """
        from .models import DataSeries
        from .serializers import DataSeriesListSerializer

        compound_id = request.query_params.get('compound')
        protocol_id = request.query_params.get('protocol')
        status_filter = request.query_params.get('status', 'valid')

        if not compound_id or not protocol_id:
            return Response(
                {'error': 'Both compound and protocol parameters are required'},
                status=status.HTTP_400_BAD_REQUEST
            )

        queryset = DataSeries.objects.select_related(
            'assay',
            'assay__protocol',
            'compound',
            'analysis',
            'dilution_series',
        ).filter(
            compound_id=compound_id,
            assay__protocol_id=protocol_id,
        )

        if status_filter:
            queryset = queryset.filter(analysis__status=status_filter)

        # Order by assay date descending
        queryset = queryset.order_by('-assay__created_at')

        serializer = DataSeriesListSerializer(queryset, many=True)

        # Include compound and protocol info in response
        compound_info = None
        protocol_info = None
        if queryset.exists():
            first = queryset.first()
            if first.compound:
                compound_info = {
                    'id': str(first.compound.id),
                    'formatted_id': first.compound.formatted_id,
                    'smiles': first.compound.smiles or first.compound.rdkit_smiles,
                }
            protocol_info = {
                'id': str(first.assay.protocol.id),
                'name': first.assay.protocol.name,
            }

        return Response({
            'compound': compound_info,
            'protocol': protocol_info,
            'count': queryset.count(),
            'data_series': serializer.data,
        })
