"""
Aggregation module for compound assay data.

Provides functions to aggregate KPI values across compounds and protocols,
supporting both compact (one row per compound) and long (one row per measurement)
output formats.
"""

import math
from collections import defaultdict
from typing import Any

from django.db.models import QuerySet

from .models import DataSeries, Protocol


def geometric_mean(values: list[float]) -> float | None:
    """
    Calculate the geometric mean of positive values.

    Uses the formula: 10^(sum(log10(values))/count)

    Args:
        values: List of positive numeric values

    Returns:
        Geometric mean, or None if no valid values
    """
    # Filter out non-positive values (can't take log of zero/negative)
    positive_values = [v for v in values if v is not None and v > 0]

    if not positive_values:
        return None

    log_sum = sum(math.log10(v) for v in positive_values)
    return math.pow(10, log_sum / len(positive_values))


def standard_deviation(values: list[float]) -> float | None:
    """
    Calculate the sample standard deviation.

    Args:
        values: List of numeric values

    Returns:
        Sample standard deviation, or None if fewer than 2 values
    """
    valid_values = [v for v in values if v is not None]

    if len(valid_values) < 2:
        return None

    n = len(valid_values)
    mean = sum(valid_values) / n
    variance = sum((v - mean) ** 2 for v in valid_values) / (n - 1)
    return math.sqrt(variance)


def aggregate_kpi_values(
    values: list[float],
    aggregations: list[str]
) -> dict[str, Any]:
    """
    Compute requested aggregation statistics for a list of KPI values.

    Args:
        values: List of KPI values
        aggregations: List of aggregation types to compute
                     ('geomean', 'count', 'stdev', 'list')

    Returns:
        Dictionary with requested aggregation results
    """
    result = {}
    valid_values = [v for v in values if v is not None]

    if 'geomean' in aggregations:
        result['geomean'] = geometric_mean(valid_values)

    if 'count' in aggregations:
        result['count'] = len(valid_values)

    if 'stdev' in aggregations:
        result['stdev'] = standard_deviation(valid_values)

    if 'list' in aggregations:
        # Format as comma-separated string with 2 decimal places
        formatted = [f"{v:.2f}" if isinstance(v, float) else str(v)
                     for v in valid_values]
        result['list'] = ', '.join(formatted)

    return result


def build_data_series_queryset(predicates: dict) -> QuerySet[DataSeries]:
    """
    Build an optimized queryset for data series based on filter predicates.

    Args:
        predicates: Dictionary of filter criteria:
            - targets: List of target UUIDs
            - compounds: List of compound UUIDs
            - compound_search: Text search for compound formatted_id
            - protocols: List of protocol UUIDs
            - status: Analysis status filter ('valid', 'invalid', 'unassigned')

    Returns:
        Filtered and optimized QuerySet
    """
    queryset = DataSeries.objects.select_related(
        'compound',
        'compound__target',
        'assay',
        'assay__protocol',
        'assay__target',
        'analysis',
        'dilution_series',
    )

    # Filter by targets (compound's registered target or assay's target)
    targets = predicates.get('targets', [])
    if targets:
        from django.db.models import Q
        queryset = queryset.filter(
            Q(compound__target_id__in=targets) |
            Q(assay__target_id__in=targets)
        )

    # Filter by specific compounds
    compounds = predicates.get('compounds', [])
    if compounds:
        queryset = queryset.filter(compound_id__in=compounds)

    # Filter by compound search (formatted_id pattern)
    # Supports comma-separated list of IDs (e.g., "NCL-00026123, NCL-00026124")
    compound_search = predicates.get('compound_search', '')
    if compound_search:
        import re
        # Find all NCL numbers in the search string (handles comma-separated lists)
        matches = re.findall(r'NCL-?(\d+)', compound_search, re.IGNORECASE)
        if matches:
            reg_numbers = [int(m) for m in matches]
            queryset = queryset.filter(compound__reg_number__in=reg_numbers)
        else:
            # Fall back to compound_name search (handles non-NCL formatted IDs)
            queryset = queryset.filter(compound_name__icontains=compound_search)

    # Filter by protocols
    protocols = predicates.get('protocols', [])
    if protocols:
        queryset = queryset.filter(assay__protocol_id__in=protocols)

    # Filter by analysis status (default to 'valid')
    status = predicates.get('status', 'valid')
    if status:
        queryset = queryset.filter(analysis__status=status)

    return queryset


def extract_kpi_value(data_series: DataSeries) -> float | None:
    """
    Extract the KPI value from a data series' analysis result.

    Args:
        data_series: DataSeries instance with loaded analysis

    Returns:
        KPI value or None if not available
    """
    if not data_series.analysis:
        return None

    results = data_series.analysis.results or {}
    kpi_key = results.get('KPI')

    if kpi_key and kpi_key in results:
        value = results[kpi_key]
        if isinstance(value, (int, float)):
            return float(value)

    return None


def aggregate_compact(
    queryset: QuerySet[DataSeries],
    aggregations: list[str]
) -> dict:
    """
    Aggregate data series into compact format (one row per compound).

    Each compound gets one row with columns for each protocol containing
    aggregated KPI statistics.

    Args:
        queryset: Filtered DataSeries queryset
        aggregations: List of aggregation functions to apply

    Returns:
        Dictionary with:
            - meta: Summary statistics
            - protocols: List of protocol info
            - data: List of compound rows with protocol aggregations
    """
    # Group data by compound and protocol
    # Structure: {compound_id: {protocol_id: [kpi_values]}}
    compound_protocol_values = defaultdict(lambda: defaultdict(list))
    compound_info = {}
    protocol_ids = set()

    for ds in queryset:
        if not ds.compound:
            continue

        compound_id = str(ds.compound.id)
        protocol_id = str(ds.assay.protocol_id)

        kpi_value = extract_kpi_value(ds)
        if kpi_value is not None:
            compound_protocol_values[compound_id][protocol_id].append(kpi_value)

        protocol_ids.add(ds.assay.protocol_id)

        # Store compound info (first occurrence)
        if compound_id not in compound_info:
            compound_info[compound_id] = {
                'compound_id': compound_id,
                'formatted_id': ds.compound.formatted_id,
                'smiles': ds.compound.smiles or ds.compound.rdkit_smiles,
                'target_name': ds.compound.target.name if ds.compound.target else None,
            }

    # Fetch protocol info
    protocols = list(Protocol.objects.filter(id__in=protocol_ids).values('id', 'name'))
    protocol_list = [{'id': str(p['id']), 'name': p['name']} for p in protocols]

    # Build result rows
    data = []
    total_measurements = 0

    for compound_id, protocol_values in compound_protocol_values.items():
        row = compound_info[compound_id].copy()
        row['protocols'] = {}

        for protocol_id, values in protocol_values.items():
            total_measurements += len(values)
            row['protocols'][protocol_id] = aggregate_kpi_values(values, aggregations)

        data.append(row)

    # Sort by formatted_id
    data.sort(key=lambda x: x.get('formatted_id', ''))

    return {
        'meta': {
            'compound_count': len(data),
            'protocol_count': len(protocol_list),
            'total_measurements': total_measurements,
        },
        'protocols': protocol_list,
        'data': data,
    }


def aggregate_medium(
    queryset: QuerySet[DataSeries],
    aggregations: list[str]
) -> dict:
    """
    Aggregate data series into medium format (one row per compound-protocol pair).

    Each row represents aggregated KPIs for a single compound tested with a
    single protocol.

    Args:
        queryset: Filtered DataSeries queryset
        aggregations: List of aggregation functions to apply

    Returns:
        Dictionary with:
            - meta: Summary statistics
            - data: List of compound-protocol rows with aggregations
    """
    # Group data by compound and protocol
    # Structure: {(compound_id, protocol_id): {'info': {...}, 'values': [...]}}
    compound_protocol_data = defaultdict(lambda: {'info': None, 'values': []})
    compound_ids = set()
    protocol_ids = set()

    for ds in queryset:
        if not ds.compound:
            continue

        compound_id = str(ds.compound.id)
        protocol_id = str(ds.assay.protocol_id)
        key = (compound_id, protocol_id)

        kpi_value = extract_kpi_value(ds)
        if kpi_value is not None:
            compound_protocol_data[key]['values'].append(kpi_value)

        compound_ids.add(ds.compound.id)
        protocol_ids.add(ds.assay.protocol_id)

        # Store info (first occurrence)
        if compound_protocol_data[key]['info'] is None:
            compound_protocol_data[key]['info'] = {
                'compound_id': compound_id,
                'formatted_id': ds.compound.formatted_id,
                'smiles': ds.compound.smiles or ds.compound.rdkit_smiles,
                'target_name': ds.compound.target.name if ds.compound.target else None,
                'protocol_id': protocol_id,
                'protocol_name': ds.assay.protocol.name,
            }

    # Build result rows
    data = []
    total_measurements = 0

    for key, item in compound_protocol_data.items():
        if item['info'] is None:
            continue

        values = item['values']
        total_measurements += len(values)

        row = item['info'].copy()
        row.update(aggregate_kpi_values(values, aggregations))
        data.append(row)

    # Sort by formatted_id, then protocol
    data.sort(key=lambda x: (x.get('formatted_id', ''), x.get('protocol_name', '')))

    return {
        'meta': {
            'compound_count': len(compound_ids),
            'protocol_count': len(protocol_ids),
            'total_measurements': total_measurements,
        },
        'data': data,
    }


def aggregate_long(
    queryset: QuerySet[DataSeries],
    aggregations: list[str]
) -> dict:
    """
    Aggregate data series into long format (one row per measurement).

    Each row represents a single KPI measurement with compound, protocol,
    and assay information.

    Args:
        queryset: Filtered DataSeries queryset
        aggregations: Not used in long format (included for API consistency)

    Returns:
        Dictionary with:
            - meta: Summary statistics
            - data: List of measurement rows
    """
    data = []
    compound_ids = set()
    protocol_ids = set()

    for ds in queryset:
        kpi_value = extract_kpi_value(ds)

        row = {
            'data_series_id': str(ds.id),
            'compound_id': str(ds.compound.id) if ds.compound else None,
            'formatted_id': ds.compound.formatted_id if ds.compound else None,
            'compound_name': ds.compound_name,
            'smiles': (ds.compound.smiles or ds.compound.rdkit_smiles) if ds.compound else None,
            'target_name': ds.compound.target.name if ds.compound and ds.compound.target else None,
            'protocol_id': str(ds.assay.protocol_id),
            'protocol_name': ds.assay.protocol.name,
            'assay_id': str(ds.assay.id),
            'assay_date': ds.assay.created_at.isoformat() if ds.assay.created_at else None,
            'kpi_value': kpi_value,
            'status': ds.analysis.status if ds.analysis else None,
        }
        data.append(row)

        if ds.compound:
            compound_ids.add(ds.compound.id)
        protocol_ids.add(ds.assay.protocol_id)

    # Sort by formatted_id, then protocol
    data.sort(key=lambda x: (x.get('formatted_id') or '', x.get('protocol_name') or ''))

    return {
        'meta': {
            'compound_count': len(compound_ids),
            'protocol_count': len(protocol_ids),
            'total_measurements': len(data),
        },
        'data': data,
    }
