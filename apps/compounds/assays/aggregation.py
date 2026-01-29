"""
Aggregation module for compound assay data.

Provides functions to aggregate KPI values across compounds and protocols,
supporting both compact (one row per compound) and long (one row per measurement)
output formats.
"""

import math
import re
from collections import defaultdict
from typing import Any

from django.db.models import QuerySet

from compounds.formatting import get_compound_pattern
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
            - group_by_batch: If True, batch info will be used for grouping

    Returns:
        Filtered and optimized QuerySet
    """
    queryset = DataSeries.objects.select_related(
        'compound',
        'compound__target',
        'batch',  # Include batch for batch-aware aggregation
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
    # Supports comma-separated list of IDs (e.g., "PREFIX-00026123, PREFIX-00026124")
    compound_search = predicates.get('compound_search', '')
    if compound_search:
        # Find all compound IDs in the search string (handles comma-separated lists)
        compound_pattern = get_compound_pattern(capturing=True)
        matches = compound_pattern.findall(compound_search)
        if matches:
            reg_numbers = [int(m) for m in matches]
            queryset = queryset.filter(compound__reg_number__in=reg_numbers)
        else:
            # Fall back to compound_name search (handles non-formatted IDs)
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


def get_kpi_unit(data_series: DataSeries) -> str | None:
    """
    Get KPI unit from data series.

    Checks multiple sources in priority order:
    1. Explicitly stored kpi_unit in analysis results (from imports)
    2. DilutionSeries concentration unit (for dose-response fitted data)
    3. Protocol's preferred_dilutions unit (fallback)

    Args:
        data_series: DataSeries instance with loaded analysis and dilution_series

    Returns:
        Unit string (e.g., 'nM', 'uM', 'mM') or None if not available
    """
    # Priority 1: Check for explicitly stored unit in analysis results
    if data_series.analysis and data_series.analysis.results:
        stored_unit = data_series.analysis.results.get('kpi_unit')
        if stored_unit:
            return stored_unit

    # Priority 2: From DilutionSeries (dose-response fitted data)
    if data_series.dilution_series:
        return data_series.dilution_series.unit

    # Priority 3: From protocol's preferred_dilutions (fallback)
    if data_series.assay and data_series.assay.protocol:
        if data_series.assay.protocol.preferred_dilutions:
            return data_series.assay.protocol.preferred_dilutions.unit

    return None


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


def extract_kpi_with_unit(data_series: DataSeries) -> tuple[float | None, str | None]:
    """
    Extract the KPI value and unit from a data series.

    Args:
        data_series: DataSeries instance with loaded analysis and dilution_series

    Returns:
        Tuple of (kpi_value, unit_string)
    """
    value = extract_kpi_value(data_series)
    unit = get_kpi_unit(data_series)
    return value, unit


def aggregate_compact(
    queryset: QuerySet[DataSeries],
    aggregations: list[str],
    group_by_batch: bool = False,
    include_tested_no_data: bool = False
) -> dict:
    """
    Aggregate data series into compact format (one row per compound or compound/batch).

    Each compound (or compound/batch when group_by_batch=True) gets one row with
    columns for each protocol containing aggregated KPI statistics.

    Args:
        queryset: Filtered DataSeries queryset
        aggregations: List of aggregation functions to apply
        group_by_batch: If True, create separate rows for each batch
        include_tested_no_data: If True, include compounds that were tested but
            have no valid KPI values (shown with count=0)

    Returns:
        Dictionary with:
            - meta: Summary statistics
            - protocols: List of protocol info (including kpi_unit)
            - data: List of compound rows with protocol aggregations
            - group_by_batch: Whether results are grouped by batch
    """
    # Group data by compound (and optionally batch) and protocol
    # Structure: {group_key: {protocol_id: [kpi_values]}}
    # group_key is (compound_id,) or (compound_id, batch_id)
    group_protocol_values = defaultdict(lambda: defaultdict(list))
    # Track which protocols each compound was tested on (regardless of KPI validity)
    tested_protocols: dict[tuple, set[str]] = defaultdict(set)
    group_info = {}
    protocol_ids = set()
    protocol_units: dict[str, str | None] = {}  # Track unit per protocol

    for ds in queryset:
        if not ds.compound:
            continue

        compound_id = str(ds.compound.id)
        protocol_id = str(ds.assay.protocol_id)

        # Determine group key based on group_by_batch setting
        if group_by_batch:
            batch_id = str(ds.batch_id) if ds.batch_id else None
            group_key = (compound_id, batch_id)
        else:
            group_key = (compound_id,)

        # Always track that this compound was tested on this protocol
        tested_protocols[group_key].add(protocol_id)

        kpi_value, kpi_unit = extract_kpi_with_unit(ds)
        if kpi_value is not None:
            group_protocol_values[group_key][protocol_id].append(kpi_value)

        # Track unit per protocol (use first encountered non-null unit)
        if protocol_id not in protocol_units and kpi_unit:
            protocol_units[protocol_id] = kpi_unit

        protocol_ids.add(ds.assay.protocol_id)

        # Store group info (first occurrence)
        if group_key not in group_info:
            info = {
                'compound_id': compound_id,
                'formatted_id': ds.compound.formatted_id,
                'smiles': ds.compound.smiles or ds.compound.rdkit_smiles,
                'target_name': ds.compound.target.name if ds.compound.target else None,
            }
            if group_by_batch:
                info['batch_id'] = str(ds.batch_id) if ds.batch_id else None
                info['batch_number'] = ds.batch.batch_number if ds.batch else None
            group_info[group_key] = info

    # Fetch protocol info and include kpi_unit
    protocols = list(Protocol.objects.filter(id__in=protocol_ids).values('id', 'name'))
    protocol_list = [
        {
            'id': str(p['id']),
            'name': p['name'],
            'kpi_unit': protocol_units.get(str(p['id'])),
        }
        for p in protocols
    ]

    # Build result rows
    data = []
    total_measurements = 0

    # Determine which groups to include:
    # - If include_tested_no_data=True, iterate over all tested compounds (group_info)
    # - Otherwise, only include compounds with at least one valid KPI value
    groups_to_include = group_info.keys() if include_tested_no_data else group_protocol_values.keys()

    for group_key in groups_to_include:
        if group_key not in group_info:
            continue

        row = group_info[group_key].copy()
        row['protocols'] = {}

        # Get protocols this compound was tested on
        protocols_tested = tested_protocols.get(group_key, set())

        # For include_tested_no_data mode, include all tested protocols (even with count=0)
        # Otherwise, only include protocols with valid KPI values
        if include_tested_no_data:
            protocols_to_include = protocols_tested
        else:
            protocols_to_include = set(group_protocol_values.get(group_key, {}).keys())

        for protocol_id in protocols_to_include:
            values = group_protocol_values.get(group_key, {}).get(protocol_id, [])
            total_measurements += len(values)
            row['protocols'][protocol_id] = aggregate_kpi_values(values, aggregations)

        # Only include row if it has at least one protocol
        if row['protocols']:
            data.append(row)

    # Sort by formatted_id, then batch_number if grouping by batch
    if group_by_batch:
        data.sort(key=lambda x: (x.get('formatted_id', ''), x.get('batch_number') or 0))
    else:
        data.sort(key=lambda x: x.get('formatted_id', ''))

    # Count unique compounds (may differ from row count when group_by_batch=True)
    unique_compounds = len(set(row['compound_id'] for row in data))

    return {
        'meta': {
            'compound_count': unique_compounds,
            'row_count': len(data),  # May differ when group_by_batch=True
            'protocol_count': len(protocol_list),
            'total_measurements': total_measurements,
            'group_by_batch': group_by_batch,
            'include_tested_no_data': include_tested_no_data,
        },
        'protocols': protocol_list,
        'data': data,
    }


def aggregate_medium(
    queryset: QuerySet[DataSeries],
    aggregations: list[str],
    group_by_batch: bool = False,
    include_tested_no_data: bool = False
) -> dict:
    """
    Aggregate data series into medium format (one row per compound-protocol pair,
    or compound-batch-protocol when group_by_batch=True).

    Each row represents aggregated KPIs for a single compound (or compound/batch)
    tested with a single protocol.

    Args:
        queryset: Filtered DataSeries queryset
        aggregations: List of aggregation functions to apply
        group_by_batch: If True, create separate rows for each batch
        include_tested_no_data: If True, include compounds that were tested but
            have no valid KPI values (shown with count=0)

    Returns:
        Dictionary with:
            - meta: Summary statistics
            - data: List of compound-protocol rows with aggregations (including kpi_unit)
    """
    # Group data by compound (and optionally batch) and protocol
    # Structure: {key: {'info': {...}, 'values': [], 'unit': ...}}
    group_protocol_data = defaultdict(lambda: {'info': None, 'values': [], 'unit': None})
    compound_ids = set()
    protocol_ids = set()

    for ds in queryset:
        if not ds.compound:
            continue

        compound_id = str(ds.compound.id)
        protocol_id = str(ds.assay.protocol_id)

        # Determine group key based on group_by_batch setting
        if group_by_batch:
            batch_id = str(ds.batch_id) if ds.batch_id else None
            key = (compound_id, batch_id, protocol_id)
        else:
            key = (compound_id, protocol_id)

        kpi_value, kpi_unit = extract_kpi_with_unit(ds)
        if kpi_value is not None:
            group_protocol_data[key]['values'].append(kpi_value)

        # Track unit (use first encountered non-null unit)
        if group_protocol_data[key]['unit'] is None and kpi_unit:
            group_protocol_data[key]['unit'] = kpi_unit

        compound_ids.add(ds.compound.id)
        protocol_ids.add(ds.assay.protocol_id)

        # Store info (first occurrence) - always track tested compounds
        if group_protocol_data[key]['info'] is None:
            info = {
                'compound_id': compound_id,
                'formatted_id': ds.compound.formatted_id,
                'smiles': ds.compound.smiles or ds.compound.rdkit_smiles,
                'target_name': ds.compound.target.name if ds.compound.target else None,
                'protocol_id': protocol_id,
                'protocol_name': ds.assay.protocol.name,
            }
            if group_by_batch:
                info['batch_id'] = str(ds.batch_id) if ds.batch_id else None
                info['batch_number'] = ds.batch.batch_number if ds.batch else None
            group_protocol_data[key]['info'] = info

    # Build result rows
    data = []
    total_measurements = 0

    for key, item in group_protocol_data.items():
        if item['info'] is None:
            continue

        values = item['values']

        # Skip rows with no valid KPI values unless include_tested_no_data is True
        if not values and not include_tested_no_data:
            continue

        total_measurements += len(values)

        row = item['info'].copy()
        row['kpi_unit'] = item['unit']  # Include unit in row
        row.update(aggregate_kpi_values(values, aggregations))
        data.append(row)

    # Sort by formatted_id, then batch_number (if applicable), then protocol
    if group_by_batch:
        data.sort(key=lambda x: (
            x.get('formatted_id', ''),
            x.get('batch_number') or 0,
            x.get('protocol_name', '')
        ))
    else:
        data.sort(key=lambda x: (x.get('formatted_id', ''), x.get('protocol_name', '')))

    return {
        'meta': {
            'compound_count': len(compound_ids),
            'row_count': len(data),  # May differ when group_by_batch=True
            'protocol_count': len(protocol_ids),
            'total_measurements': total_measurements,
            'group_by_batch': group_by_batch,
            'include_tested_no_data': include_tested_no_data,
        },
        'data': data,
    }


def aggregate_long(
    queryset: QuerySet[DataSeries],
    aggregations: list[str],
    group_by_batch: bool = False,
    include_tested_no_data: bool = False
) -> dict:
    """
    Aggregate data series into long format (one row per measurement).

    Each row represents a single KPI measurement with compound, protocol,
    and assay information. When group_by_batch=True, batch info is included.

    Args:
        queryset: Filtered DataSeries queryset
        aggregations: Not used in long format (included for API consistency)
        group_by_batch: If True, include batch information in rows
        include_tested_no_data: If True, include data series with no valid KPI values
            (long format already includes all rows, this is for API consistency)

    Returns:
        Dictionary with:
            - meta: Summary statistics
            - data: List of measurement rows (including kpi_unit and optionally batch)
    """
    data = []
    compound_ids = set()
    protocol_ids = set()

    for ds in queryset:
        kpi_value, kpi_unit = extract_kpi_with_unit(ds)

        # In long format, we show all rows unless filtering by valid KPI only
        # When include_tested_no_data=False, skip rows without valid KPI
        if kpi_value is None and not include_tested_no_data:
            # Still track the compound/protocol as tested
            if ds.compound:
                compound_ids.add(ds.compound.id)
            protocol_ids.add(ds.assay.protocol_id)
            continue

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
            'kpi_unit': kpi_unit,
            'status': ds.analysis.status if ds.analysis else None,
        }

        # Include batch info (always in long format, but flagged by group_by_batch for UI)
        row['batch_id'] = str(ds.batch_id) if ds.batch_id else None
        row['batch_number'] = ds.batch.batch_number if ds.batch else None

        data.append(row)

        if ds.compound:
            compound_ids.add(ds.compound.id)
        protocol_ids.add(ds.assay.protocol_id)

    # Sort by formatted_id, then batch_number, then protocol
    data.sort(key=lambda x: (
        x.get('formatted_id') or '',
        x.get('batch_number') or 0,
        x.get('protocol_name') or ''
    ))

    return {
        'meta': {
            'compound_count': len(compound_ids),
            'protocol_count': len(protocol_ids),
            'total_measurements': len(data),
            'group_by_batch': group_by_batch,
            'include_tested_no_data': include_tested_no_data,
        },
        'data': data,
    }
