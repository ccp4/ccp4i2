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
            "output_format": "compact" | "long",
            "aggregations": ["geomean", "count", "stdev", "list"]
        }

        Returns:
            Compact format: One row per compound with protocol columns
            Long format: One row per measurement
        """
        predicates = request.data.get('predicates', {})
        output_format = request.data.get('output_format', 'compact')
        aggregations = request.data.get('aggregations', ['geomean', 'count'])

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
            result = aggregate_compact(queryset, aggregations)
        elif output_format == 'medium':
            result = aggregate_medium(queryset, aggregations)
        else:
            result = aggregate_long(queryset, aggregations)

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
