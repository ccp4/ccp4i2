"""
Unit tests for the compound assay data aggregation module.

Tests validate:
1. Compound-centric mode (targets specified): ALL compounds for those targets appear
2. DataSeries-centric mode (only protocols): only compounds with data appear
3. Valid-only aggregation: only valid measurements in geomean/stdev
4. Count tracking: tested, no_analysis, invalid, unassigned counts
5. Protocol filtering in compound-centric mode
"""

import pytest
import uuid
from decimal import Decimal
from unittest.mock import MagicMock, patch

from compounds.assays.aggregation import (
    geometric_mean,
    standard_deviation,
    log_standard_deviation,
    aggregate_kpi_values,
    should_use_compound_centric_mode,
    extract_kpi_value,
)


class TestGeometricMean:
    """Test geometric mean calculation."""

    def test_single_value(self):
        """Single value should return itself."""
        result = geometric_mean([100.0])
        assert abs(result - 100.0) < 0.01

    def test_two_values(self):
        """Geometric mean of 10 and 1000 is 100."""
        result = geometric_mean([10.0, 1000.0])
        assert abs(result - 100.0) < 0.01

    def test_filters_non_positive(self):
        """Non-positive values should be filtered out."""
        result = geometric_mean([100.0, 0, -50, None])
        assert abs(result - 100.0) < 0.01

    def test_empty_returns_none(self):
        """Empty list returns None."""
        result = geometric_mean([])
        assert result is None

    def test_all_invalid_returns_none(self):
        """All invalid values returns None."""
        result = geometric_mean([0, -1, None])
        assert result is None


class TestStandardDeviation:
    """Test standard deviation calculation."""

    def test_two_values(self):
        """Sample stdev of [10, 20] is ~7.07."""
        result = standard_deviation([10.0, 20.0])
        # Sample stdev = sqrt(((10-15)^2 + (20-15)^2) / 1) = sqrt(50) = 7.07
        assert abs(result - 7.07) < 0.1

    def test_single_value_returns_none(self):
        """Single value can't have stdev."""
        result = standard_deviation([100.0])
        assert result is None

    def test_filters_none(self):
        """None values should be filtered."""
        result = standard_deviation([10.0, 20.0, None])
        assert abs(result - 7.07) < 0.1


class TestLogStandardDeviation:
    """Test log-space standard deviation."""

    def test_basic_calculation(self):
        """Test log-space stdev for pConc display."""
        # Values: 100nM and 10000nM
        # pConc: -log10(100e-9) = 7, -log10(10000e-9) = 5
        # stdev of [7, 5] in pConc space
        result = log_standard_deviation([100.0, 10000.0])
        # Mean = 6, stdev = sqrt(((7-6)^2 + (5-6)^2) / 1) = sqrt(2) = 1.414
        assert result is not None
        assert abs(result - 1.414) < 0.1

    def test_single_value_returns_none(self):
        """Single value can't have stdev."""
        result = log_standard_deviation([100.0])
        assert result is None


class TestAggregateKpiValues:
    """Test the aggregation function."""

    def test_geomean_aggregation(self):
        """Test geomean calculation."""
        result = aggregate_kpi_values([10.0, 1000.0], ['geomean'])
        assert 'geomean' in result
        assert abs(result['geomean'] - 100.0) < 0.01

    def test_count_aggregation(self):
        """Test count of valid values."""
        result = aggregate_kpi_values([10.0, 20.0, None], ['count'])
        assert result['count'] == 2

    def test_stdev_includes_log(self):
        """stdev aggregation should also include stdev_log."""
        result = aggregate_kpi_values([10.0, 20.0], ['stdev'])
        assert 'stdev' in result
        assert 'stdev_log' in result

    def test_list_aggregation(self):
        """Test list formatting."""
        result = aggregate_kpi_values([10.0, 20.5], ['list'])
        assert result['list'] == '10.00, 20.50'

    def test_multiple_aggregations(self):
        """Test multiple aggregations at once."""
        result = aggregate_kpi_values([10.0, 1000.0], ['geomean', 'count', 'stdev'])
        assert 'geomean' in result
        assert 'count' in result
        assert 'stdev' in result
        assert result['count'] == 2


class TestQueryModeDecision:
    """Test the query mode decision logic."""

    def test_targets_triggers_compound_centric(self):
        """When targets specified, use compound-centric mode."""
        predicates = {'targets': ['uuid-1']}
        assert should_use_compound_centric_mode(predicates) is True

    def test_empty_targets_uses_dataseries_centric(self):
        """When targets empty, use dataseries-centric mode."""
        predicates = {'targets': []}
        assert should_use_compound_centric_mode(predicates) is False

    def test_no_targets_uses_dataseries_centric(self):
        """When no targets key, use dataseries-centric mode."""
        predicates = {'protocols': ['uuid-1']}
        assert should_use_compound_centric_mode(predicates) is False

    def test_targets_with_protocols_still_compound_centric(self):
        """Targets + protocols still uses compound-centric."""
        predicates = {'targets': ['uuid-1'], 'protocols': ['uuid-2']}
        assert should_use_compound_centric_mode(predicates) is True

    def test_compounds_specified_uses_compound_centric(self):
        """When compounds explicitly specified, use compound-centric mode."""
        predicates = {'compounds': ['uuid-1', 'uuid-2']}
        assert should_use_compound_centric_mode(predicates) is True

    def test_compounds_with_targets_still_compound_centric(self):
        """Compounds + targets still uses compound-centric (compounds take priority)."""
        predicates = {'compounds': ['uuid-1'], 'targets': ['uuid-2']}
        assert should_use_compound_centric_mode(predicates) is True

    def test_compounds_with_protocols_uses_compound_centric(self):
        """Compounds + protocols uses compound-centric (compounds are definitive)."""
        predicates = {'compounds': ['uuid-1'], 'protocols': ['uuid-2']}
        assert should_use_compound_centric_mode(predicates) is True


class TestExtractKpiValue:
    """Test KPI value extraction from data series."""

    def test_valid_analysis_returns_value(self):
        """Valid analysis should return KPI value."""
        ds = MagicMock()
        ds.analysis = MagicMock()
        ds.analysis.status = 'valid'
        ds.analysis.results = {'KPI': 'EC50', 'EC50': 123.45}

        result = extract_kpi_value(ds, valid_only=True)
        assert result == 123.45

    def test_invalid_analysis_returns_none_when_valid_only(self):
        """Invalid analysis returns None when valid_only=True."""
        ds = MagicMock()
        ds.analysis = MagicMock()
        ds.analysis.status = 'invalid'
        ds.analysis.results = {'KPI': 'EC50', 'EC50': 123.45}

        result = extract_kpi_value(ds, valid_only=True)
        assert result is None

    def test_invalid_analysis_returns_value_when_not_valid_only(self):
        """Invalid analysis returns value when valid_only=False."""
        ds = MagicMock()
        ds.analysis = MagicMock()
        ds.analysis.status = 'invalid'
        ds.analysis.results = {'KPI': 'EC50', 'EC50': 123.45}

        result = extract_kpi_value(ds, valid_only=False)
        assert result == 123.45

    def test_no_analysis_returns_none(self):
        """No analysis returns None."""
        ds = MagicMock()
        ds.analysis = None

        result = extract_kpi_value(ds)
        assert result is None

    def test_unassigned_analysis_returns_none_when_valid_only(self):
        """Unassigned analysis returns None when valid_only=True."""
        ds = MagicMock()
        ds.analysis = MagicMock()
        ds.analysis.status = 'unassigned'
        ds.analysis.results = {'KPI': 'EC50', 'EC50': 123.45}

        result = extract_kpi_value(ds, valid_only=True)
        assert result is None


class TestAggregationBehaviorRequirements:
    """
    Tests documenting the required aggregation behavior.

    These tests verify the key requirements:
    1. When targets specified: ALL compounds for those targets appear
    2. When only protocols: only compounds with data appear
    3. Only valid KPIs go into geomean/stdev
    4. Counts distinguish valid, invalid, no_analysis, unassigned
    """

    def test_requirement_only_valid_in_geomean(self):
        """
        Requirement: Only valid measurements should go into geomean calculation.

        This is enforced by extract_kpi_value using valid_only=True by default.
        """
        # Valid analysis
        valid_ds = MagicMock()
        valid_ds.analysis = MagicMock()
        valid_ds.analysis.status = 'valid'
        valid_ds.analysis.results = {'KPI': 'EC50', 'EC50': 100.0}

        # Invalid analysis
        invalid_ds = MagicMock()
        invalid_ds.analysis = MagicMock()
        invalid_ds.analysis.status = 'invalid'
        invalid_ds.analysis.results = {'KPI': 'EC50', 'EC50': 1000.0}

        # Only valid should be extracted
        valid_value = extract_kpi_value(valid_ds)
        invalid_value = extract_kpi_value(invalid_ds)

        assert valid_value == 100.0
        assert invalid_value is None

        # Geomean should only include valid
        values = [v for v in [valid_value, invalid_value] if v is not None]
        geomean = geometric_mean(values)
        assert geomean == 100.0  # Not geometric_mean([100, 1000]) = 316

    def test_requirement_targets_triggers_compound_centric(self):
        """
        Requirement: If targets are specified, use compound-centric mode.

        This ensures ALL compounds registered to those targets appear,
        even those with no data series.
        """
        # Target specified
        assert should_use_compound_centric_mode({'targets': ['uuid-1']}) is True

        # No target specified
        assert should_use_compound_centric_mode({'protocols': ['uuid-1']}) is False
        assert should_use_compound_centric_mode({}) is False

    def test_requirement_compounds_are_definitive(self):
        """
        Requirement: If specific compounds are specified, they define the exact set.

        When compounds are explicitly specified:
        - They are the definitive set (highest priority)
        - Targets are ignored (no intersection)
        - All specified compounds appear, even with no data
        """
        # Compounds specified = compound-centric mode
        assert should_use_compound_centric_mode({'compounds': ['uuid-1']}) is True

        # Compounds + targets = still compound-centric (compounds take priority)
        assert should_use_compound_centric_mode({
            'compounds': ['uuid-1'],
            'targets': ['uuid-2']
        }) is True

        # Compounds + protocols = compound-centric (definitive compound list)
        assert should_use_compound_centric_mode({
            'compounds': ['uuid-1'],
            'protocols': ['uuid-2']
        }) is True

    def test_requirement_counts_are_tracked_separately(self):
        """
        Requirement: Track counts separately for:
        - count: valid KPIs (used in aggregations)
        - tested: total data series
        - no_analysis: data series with no analysis
        - invalid: data series with invalid analysis
        - unassigned: data series with unassigned analysis
        """
        # The aggregate_kpi_values function handles 'count' of valid values
        # The aggregate_compact/medium/long functions track tested/no_analysis/invalid/unassigned

        # Verify count only counts valid values
        result = aggregate_kpi_values([100.0, None, None], ['count'])
        assert result['count'] == 1  # Only 1 valid value

        # The other counts (tested, no_analysis, etc.) are computed in the
        # aggregate functions by examining ds.analysis status
