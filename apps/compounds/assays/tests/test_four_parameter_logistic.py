"""
Unit tests for the Four Parameter Logistic (4PL) curve fitting script.

Tests validate:
1. Correct IC50/Hill/top/bottom recovery from synthetic data
2. Return structure matches expected schema
3. Constraint modes (fix_hill, fix_top, fix_bottom, restrain_to_controls)
4. Edge cases and error handling
"""

import pytest
import numpy as np
from pathlib import Path
import sys

# Add the fitting_scripts directory to the path
fitting_scripts_dir = Path(__file__).parent.parent / 'fitting_scripts'
sys.path.insert(0, str(fitting_scripts_dir))

from four_parameter_logistic import fit, four_param_logistic


class TestFourParamLogisticFunction:
    """Test the underlying 4PL equation."""

    def test_at_zero_concentration(self):
        """At x→0, response should approach top."""
        # Very small x should give response near top
        result = four_param_logistic(1e-10, top=100, bottom=0, ic50=100, hill=1)
        assert result > 99.9, f"Expected ~100 at very low conc, got {result}"

    def test_at_infinite_concentration(self):
        """At x→∞, response should approach bottom."""
        result = four_param_logistic(1e10, top=100, bottom=0, ic50=100, hill=1)
        assert result < 0.1, f"Expected ~0 at very high conc, got {result}"

    def test_at_ic50(self):
        """At x=IC50, response should be midpoint."""
        result = four_param_logistic(100, top=100, bottom=0, ic50=100, hill=1)
        expected = 50  # Midpoint
        assert abs(result - expected) < 0.1, f"Expected {expected} at IC50, got {result}"

    def test_hill_slope_effect(self):
        """Higher Hill coefficient should give steeper curve."""
        x_values = np.array([50, 100, 200])

        # Hill=1 curve
        y_hill_1 = [four_param_logistic(x, 100, 0, 100, 1) for x in x_values]

        # Hill=2 curve (steeper)
        y_hill_2 = [four_param_logistic(x, 100, 0, 100, 2) for x in x_values]

        # At IC50, both should be 50
        assert abs(y_hill_1[1] - 50) < 0.1
        assert abs(y_hill_2[1] - 50) < 0.1

        # Away from IC50, Hill=2 should deviate more from 50
        # At x=50 (below IC50), Hill=2 should be higher (closer to top)
        assert y_hill_2[0] > y_hill_1[0], "Hill=2 should be steeper"


class TestFitFunction:
    """Test the fit() function with synthetic data."""

    def generate_synthetic_data(
        self,
        ic50: float = 100,
        hill: float = 1.0,
        top: float = 100,
        bottom: float = 0,
        noise_level: float = 0,
        n_points: int = 8
    ) -> dict:
        """Generate synthetic dose-response data."""
        # Log-spaced concentrations
        concentrations = np.logspace(
            np.log10(ic50 / 100),
            np.log10(ic50 * 100),
            n_points
        ).tolist()

        # Calculate responses using 4PL equation
        responses = [
            four_param_logistic(c, top, bottom, ic50, hill)
            for c in concentrations
        ]

        # Add noise if requested
        if noise_level > 0:
            np.random.seed(42)
            responses = [r + np.random.normal(0, noise_level) for r in responses]

        return {
            'concentrations': concentrations,
            'responses': responses,
            'controls': {'max': top, 'min': bottom},
            'parameters': {},
        }

    def test_basic_fit_returns_expected_structure(self):
        """Verify the return structure matches expected schema."""
        input_data = self.generate_synthetic_data()
        result = fit(input_data)

        # Check required fields
        assert 'fit_successful' in result
        assert 'kpi' in result
        assert result['kpi'] == 'ic50'

        # On success, should have these fields
        assert result['fit_successful'] is True
        assert 'ic50' in result
        assert 'hill_slope' in result
        assert 'top' in result
        assert 'bottom' in result
        assert 'r_squared' in result
        assert 'curve_points' in result
        assert 'flags' in result
        assert isinstance(result['flags'], list)

        # curve_points should be list of [x, y] pairs
        assert len(result['curve_points']) > 0
        assert len(result['curve_points'][0]) == 2

    def test_recovers_known_ic50(self):
        """Fit should recover known IC50 from clean synthetic data."""
        true_ic50 = 100
        input_data = self.generate_synthetic_data(ic50=true_ic50)
        result = fit(input_data)

        assert result['fit_successful'] is True
        fitted_ic50 = result['ic50']

        # Should be within 5% of true value
        error = abs(fitted_ic50 - true_ic50) / true_ic50
        assert error < 0.05, f"IC50 error {error:.1%} exceeds 5% (fitted={fitted_ic50}, true={true_ic50})"

    def test_recovers_known_hill(self):
        """Fit should recover known Hill coefficient from clean synthetic data."""
        true_hill = 1.5
        input_data = self.generate_synthetic_data(hill=true_hill)
        result = fit(input_data)

        assert result['fit_successful'] is True
        fitted_hill = result['hill_slope']

        # Should be within 10% of true value
        error = abs(fitted_hill - true_hill) / true_hill
        assert error < 0.10, f"Hill error {error:.1%} exceeds 10% (fitted={fitted_hill}, true={true_hill})"

    def test_recovers_known_top_bottom(self):
        """Fit should recover known top/bottom from clean synthetic data."""
        true_top, true_bottom = 95, 5
        input_data = self.generate_synthetic_data(top=true_top, bottom=true_bottom)
        result = fit(input_data)

        assert result['fit_successful'] is True

        # Should be within 5% of dynamic range
        dynamic_range = true_top - true_bottom
        top_error = abs(result['top'] - true_top) / dynamic_range
        bottom_error = abs(result['bottom'] - true_bottom) / dynamic_range

        assert top_error < 0.05, f"Top error {top_error:.1%} exceeds 5%"
        assert bottom_error < 0.05, f"Bottom error {bottom_error:.1%} exceeds 5%"

    def test_r_squared_high_for_clean_data(self):
        """R² should be very high (>0.99) for noise-free synthetic data."""
        input_data = self.generate_synthetic_data(noise_level=0)
        result = fit(input_data)

        assert result['fit_successful'] is True
        assert result['r_squared'] > 0.99, f"R² should be >0.99 for clean data, got {result['r_squared']}"

    def test_handles_noisy_data(self):
        """Fit should still work with moderate noise."""
        input_data = self.generate_synthetic_data(noise_level=5)
        result = fit(input_data)

        assert result['fit_successful'] is True
        # R² will be lower with noise but should still be reasonable
        assert result['r_squared'] > 0.8, f"R² should be >0.8 even with noise, got {result['r_squared']}"


class TestConstraintModes:
    """Test the constraint modes (fix_hill, fix_top, fix_bottom, restrain_to_controls)."""

    def generate_data(self) -> dict:
        """Generate test data."""
        concentrations = [10, 30, 100, 300, 1000, 3000, 10000, 30000]
        # Synthetic data with IC50=1000, Hill=1, top=100, bottom=5
        responses = [
            four_param_logistic(c, top=100, bottom=5, ic50=1000, hill=1)
            for c in concentrations
        ]
        return {
            'concentrations': concentrations,
            'responses': responses,
            'controls': {'max': 100, 'min': 5},
            'parameters': {},
        }

    def test_fix_hill_to_one(self):
        """When fix_hill=1.0, fitted Hill should be exactly 1.0."""
        input_data = self.generate_data()
        input_data['parameters'] = {'fix_hill': 1.0}
        result = fit(input_data)

        assert result['fit_successful'] is True
        assert result['hill_slope'] == 1.0, f"Hill should be exactly 1.0, got {result['hill_slope']}"

    def test_fix_top_to_control(self):
        """When fix_top=True, fitted top should equal max control."""
        input_data = self.generate_data()
        input_data['parameters'] = {'fix_top': True}
        result = fit(input_data)

        assert result['fit_successful'] is True
        assert result['top'] == 100, f"Top should be exactly 100, got {result['top']}"

    def test_fix_bottom_to_control(self):
        """When fix_bottom=True, fitted bottom should equal min control."""
        input_data = self.generate_data()
        input_data['parameters'] = {'fix_bottom': True}
        result = fit(input_data)

        assert result['fit_successful'] is True
        assert result['bottom'] == 5, f"Bottom should be exactly 5, got {result['bottom']}"

    def test_fix_both_asymptotes(self):
        """When both fix_top and fix_bottom are True."""
        input_data = self.generate_data()
        input_data['parameters'] = {'fix_top': True, 'fix_bottom': True}
        result = fit(input_data)

        assert result['fit_successful'] is True
        assert result['top'] == 100
        assert result['bottom'] == 5

    def test_restrain_to_controls(self):
        """restrain_to_controls should produce fit close to controls without hard-fixing."""
        input_data = self.generate_data()
        input_data['parameters'] = {'restrain_to_controls': True}
        result = fit(input_data)

        assert result['fit_successful'] is True
        assert result.get('restraint_applied') is True

        # Fitted values should be close to controls but not necessarily exact
        assert abs(result['top'] - 100) < 5, "Top should be close to max control"
        assert abs(result['bottom'] - 5) < 5, "Bottom should be close to min control"

    def test_restrain_skipped_for_fixed_params(self):
        """When fix_top=True, restraint should not add pseudo point for top."""
        input_data = self.generate_data()
        input_data['parameters'] = {
            'fix_top': True,
            'restrain_to_controls': True,
        }
        result = fit(input_data)

        assert result['fit_successful'] is True
        # Top is hard-fixed, bottom is restrained
        assert result['top'] == 100
        # restraint_applied should still be True because bottom was restrained
        assert result.get('restraint_applied') is True

    def test_restrain_no_effect_when_both_fixed(self):
        """When both asymptotes are fixed, restraint has no effect."""
        input_data = self.generate_data()
        input_data['parameters'] = {
            'fix_top': True,
            'fix_bottom': True,
            'restrain_to_controls': True,
        }
        result = fit(input_data)

        assert result['fit_successful'] is True
        # Both hard-fixed, so no restraint pseudo points added
        assert result.get('restraint_applied') is False


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_insufficient_data_points(self):
        """Should fail with fewer than 4 data points."""
        input_data = {
            'concentrations': [10, 100, 1000],
            'responses': [90, 50, 10],
            'controls': {},
            'parameters': {},
        }
        result = fit(input_data)

        assert result['fit_successful'] is False
        assert 'insufficient_data' in result['flags']

    def test_mismatched_lengths(self):
        """Should fail when concentrations and responses have different lengths."""
        input_data = {
            'concentrations': [10, 100, 1000, 10000],
            'responses': [90, 50, 10],  # One fewer
            'controls': {},
            'parameters': {},
        }
        result = fit(input_data)

        assert result['fit_successful'] is False
        assert 'data_mismatch' in result['flags']

    def test_handles_nan_values(self):
        """Should filter out NaN values and still fit."""
        input_data = {
            'concentrations': [10, 100, float('nan'), 1000, 3000, 10000, 30000, 100000],
            'responses': [95, 80, 50, 40, 20, 10, 5, 3],
            'controls': {'max': 100, 'min': 0},
            'parameters': {},
        }
        result = fit(input_data)

        # Should succeed after filtering NaN
        assert result['fit_successful'] is True

    def test_missing_controls_for_fix_top(self):
        """Should flag when fix_top=True but no max control provided."""
        input_data = {
            'concentrations': [10, 100, 1000, 10000, 100000],
            'responses': [95, 80, 50, 20, 5],
            'controls': {},  # No controls
            'parameters': {'fix_top': True},
        }
        result = fit(input_data)

        # Should still succeed but flag the issue
        assert 'no_max_control_for_fix_top' in result['flags']

    def test_missing_controls_for_restrain(self):
        """Should flag when restrain requested but controls missing."""
        input_data = {
            'concentrations': [10, 100, 1000, 10000, 100000],
            'responses': [95, 80, 50, 20, 5],
            'controls': {},  # No controls
            'parameters': {'restrain_to_controls': True},
        }
        result = fit(input_data)

        assert result['fit_successful'] is True
        assert 'restraint_no_max_control' in result['flags']
        assert 'restraint_no_min_control' in result['flags']
        assert result.get('restraint_applied') is False

    def test_flags_ic50_extrapolation(self):
        """Should flag when IC50 is outside data range."""
        # Data that doesn't span the IC50 well
        input_data = {
            'concentrations': [1, 3, 10, 30],  # All below IC50
            'responses': [98, 95, 90, 80],  # Haven't reached inflection
            'controls': {'max': 100, 'min': 0},
            'parameters': {},
        }
        result = fit(input_data)

        if result['fit_successful']:
            assert 'ic50_extrapolated' in result['flags']

    def test_flags_poor_fit(self):
        """Should flag when R² is low."""
        # Random-ish data that won't fit well
        input_data = {
            'concentrations': [10, 100, 1000, 10000, 100000],
            'responses': [50, 80, 20, 90, 30],  # No clear sigmoid pattern
            'controls': {'max': 100, 'min': 0},
            'parameters': {},
        }
        result = fit(input_data)

        if result['fit_successful']:
            # Should have poor fit flag if R² < 0.8
            if result['r_squared'] < 0.8:
                assert 'poor_fit' in result['flags']


class TestEndPercent:
    """Test the end_percent calculation."""

    def test_end_percent_at_full_inhibition(self):
        """end_percent should be ~0% when highest conc achieves full inhibition."""
        # Data where highest concentration reaches bottom
        concentrations = [10, 100, 1000, 10000]
        responses = [
            four_param_logistic(c, top=100, bottom=0, ic50=100, hill=1)
            for c in concentrations
        ]

        input_data = {
            'concentrations': concentrations,
            'responses': responses,
            'controls': {'max': 100, 'min': 0},
            'parameters': {},
        }
        result = fit(input_data)

        assert result['fit_successful'] is True
        assert 'end_percent' in result
        # At 10000 nM with IC50=100, should be nearly 0%
        assert result['end_percent'] < 5

    def test_end_percent_display_format(self):
        """end_percent_display should be formatted string."""
        input_data = {
            'concentrations': [10, 100, 1000, 10000],
            'responses': [95, 50, 10, 2],
            'controls': {'max': 100, 'min': 0},
            'parameters': {},
        }
        result = fit(input_data)

        assert result['fit_successful'] is True
        assert 'end_percent_display' in result
        assert isinstance(result['end_percent_display'], str)
        assert '%' in result['end_percent_display']
        assert 'at' in result['end_percent_display']
