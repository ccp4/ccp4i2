"""
Aggregation Views

REST API endpoints for compound assay data aggregation.
"""

from rest_framework import viewsets, status
from rest_framework.decorators import action
from rest_framework.parsers import JSONParser
from rest_framework.permissions import AllowAny
from rest_framework.response import Response

from compounds.registry.models import Target, MolecularPropertyThreshold
from compounds.registry.serializers import MolecularPropertyThresholdSerializer
from .models import Protocol
from .aggregation import (
    # Query mode decision
    should_use_compound_centric_mode,
    # Compound-centric functions (used when targets specified)
    build_compound_queryset,
    aggregate_compact_from_compounds,
    aggregate_medium_from_compounds,
    aggregate_long_from_compounds,
    # DataSeries-centric functions (used when only protocols specified)
    build_data_series_queryset,
    aggregate_compact,
    aggregate_medium,
    aggregate_long,
)

# Valid molecular properties that can be included in aggregation tables
VALID_MOLECULAR_PROPERTIES = {
    'molecular_weight',
    'heavy_atom_count',
    'hbd',
    'hba',
    'clogp',
    'tpsa',
    'rotatable_bonds',
    'fraction_sp3',
}


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
                "status": "valid"  // optional, defaults to "valid" (only for dataseries-centric mode)
            },
            "output_format": "compact" | "medium" | "long",
            "aggregations": ["geomean", "count", "stdev", "list"],
            "group_by_batch": false,  // optional, when true splits results by batch
            "include_tested_no_data": false,  // optional, only used in dataseries-centric mode
            "include_properties": ["molecular_weight", "clogp", ...]  // optional, molecular properties to include
        }

        Query Modes:
        ------------
        The API uses two distinct query strategies based on predicates:

        1. **Compound-centric mode** (when targets are specified):
           - Returns ALL compounds registered to the selected targets
           - Compounds with no data series still appear (with count=0)
           - If protocols are also specified, data is filtered to those protocols
           - Status filtering happens during aggregation (only valid KPIs in geomean/stdev)
           - include_tested_no_data is ignored (always true in this mode)

        2. **DataSeries-centric mode** (when only protocols specified, no targets):
           - Returns only compounds that have data for the specified protocols
           - Compounds without data for those protocols are excluded
           - Use include_tested_no_data=true to show compounds tested but with no valid KPI

        Returns:
            Compact format: One row per compound (or compound/batch) with protocol columns
            Medium format: One row per compound-protocol (or compound-batch-protocol) pair
            Long format: One row per measurement (always includes batch info)

        Data Counts:
            Each protocol entry includes:
            - count: Number of valid KPI values (used in geomean/stdev)
            - tested: Total DataSeries count
            - no_analysis: DataSeries where analysis is NULL
            - invalid: DataSeries where analysis.status='invalid'
            - unassigned: DataSeries where analysis.status='unassigned'

        When group_by_batch is true:
            - Compact/Medium: Separate rows are created for each batch of a compound
            - Long: Batch info is always included; the flag affects sorting/display hints
            - Rows include batch_id and batch_number fields
            - meta.group_by_batch indicates the grouping mode used

        When include_properties is provided:
            - Adds molecular property columns to each row
            - Valid properties: molecular_weight, heavy_atom_count, hbd, hba, clogp, tpsa, rotatable_bonds, fraction_sp3
            - Response includes 'property_thresholds' for RAG coloring
        """
        predicates = request.data.get('predicates', {})
        output_format = request.data.get('output_format', 'compact')
        aggregations = request.data.get('aggregations', ['geomean', 'count'])
        group_by_batch = request.data.get('group_by_batch', False)
        include_tested_no_data = request.data.get('include_tested_no_data', False)
        include_properties = request.data.get('include_properties', [])

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

        # Validate include_properties
        if include_properties:
            invalid_props = set(include_properties) - VALID_MOLECULAR_PROPERTIES
            if invalid_props:
                return Response(
                    {'error': f'Invalid properties: {invalid_props}. Valid options: {VALID_MOLECULAR_PROPERTIES}'},
                    status=status.HTTP_400_BAD_REQUEST
                )

        # Determine query mode based on predicates
        # Compound-centric: when targets are specified, show ALL compounds for those targets
        # DataSeries-centric: when only protocols specified, show compounds with data
        use_compound_centric = should_use_compound_centric_mode(predicates)
        protocol_ids = predicates.get('protocols', [])

        if use_compound_centric:
            # Compound-centric mode: ALL compounds for selected targets appear
            compound_queryset = build_compound_queryset(predicates)

            if output_format == 'compact':
                result = aggregate_compact_from_compounds(
                    compound_queryset, protocol_ids, aggregations,
                    group_by_batch=group_by_batch,
                    include_properties=include_properties,
                )
            elif output_format == 'medium':
                result = aggregate_medium_from_compounds(
                    compound_queryset, protocol_ids, aggregations,
                    group_by_batch=group_by_batch,
                    include_properties=include_properties,
                )
            else:
                result = aggregate_long_from_compounds(
                    compound_queryset, protocol_ids, aggregations,
                    group_by_batch=group_by_batch,
                    include_properties=include_properties,
                )
        else:
            # DataSeries-centric mode: only compounds with data for specified protocols
            queryset = build_data_series_queryset(predicates)

            if output_format == 'compact':
                result = aggregate_compact(
                    queryset, aggregations,
                    group_by_batch=group_by_batch,
                    include_tested_no_data=include_tested_no_data,
                    include_properties=include_properties,
                )
            elif output_format == 'medium':
                result = aggregate_medium(
                    queryset, aggregations,
                    group_by_batch=group_by_batch,
                    include_tested_no_data=include_tested_no_data,
                    include_properties=include_properties,
                )
            else:
                result = aggregate_long(
                    queryset, aggregations,
                    group_by_batch=group_by_batch,
                    include_tested_no_data=include_tested_no_data,
                    include_properties=include_properties,
                )

        # Include RAG thresholds if properties were requested
        if include_properties:
            thresholds = MolecularPropertyThreshold.objects.filter(
                property_name__in=include_properties,
                enabled=True
            )
            result['property_thresholds'] = MolecularPropertyThresholdSerializer(
                thresholds, many=True
            ).data

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
    def molecular_properties(self, request):
        """
        List available molecular properties and their RAG thresholds.

        Returns:
            {
                "properties": [
                    {
                        "name": "clogp",
                        "display_name": "cLogP",
                        "description": "Calculated LogP - lipophilicity indicator"
                    },
                    ...
                ],
                "thresholds": [
                    {
                        "property_name": "clogp",
                        "direction": "above",
                        "amber_threshold": 4,
                        "red_threshold": 5,
                        "enabled": true
                    },
                    ...
                ]
            }
        """
        properties = [
            {'name': 'molecular_weight', 'display_name': 'MW', 'description': 'Molecular weight (Da)'},
            {'name': 'heavy_atom_count', 'display_name': 'HAC', 'description': 'Heavy atom count (for ligand efficiency)'},
            {'name': 'hbd', 'display_name': 'HBD', 'description': 'Hydrogen bond donors'},
            {'name': 'hba', 'display_name': 'HBA', 'description': 'Hydrogen bond acceptors'},
            {'name': 'clogp', 'display_name': 'cLogP', 'description': 'Calculated LogP (lipophilicity)'},
            {'name': 'tpsa', 'display_name': 'TPSA', 'description': 'Topological polar surface area (Å²)'},
            {'name': 'rotatable_bonds', 'display_name': 'RotB', 'description': 'Rotatable bonds'},
            {'name': 'fraction_sp3', 'display_name': 'Fsp³', 'description': 'Fraction of sp³ carbons'},
        ]

        thresholds = MolecularPropertyThreshold.objects.all()
        threshold_data = MolecularPropertyThresholdSerializer(thresholds, many=True).data

        return Response({
            'properties': properties,
            'thresholds': threshold_data,
        })

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

        # Filter by analysis status
        # Options: 'valid', 'invalid', 'unassigned', 'no_analysis', or '' (all)
        if status_filter:
            if status_filter == 'no_analysis':
                queryset = queryset.filter(analysis__isnull=True)
            else:
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
