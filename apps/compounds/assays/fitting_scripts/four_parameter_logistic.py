"""
Four Parameter Logistic (4PL) Curve Fitting

Standard sigmoidal dose-response curve fitting using the Hill-Langmuir equation.

This script implements the fit() interface for the FittingMethod model.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr


def four_param_logistic(x, top, bottom, ic50, hill):
    """
    Four parameter logistic function.

    Args:
        x: Concentration values
        top: Maximum response (asymptote at high concentrations)
        bottom: Minimum response (asymptote at low concentrations)
        ic50: Concentration at 50% response
        hill: Hill coefficient (slope)

    Returns:
        Response values
    """
    return bottom + (top - bottom) / (1 + (x / ic50) ** hill)


def fit(input_data: dict) -> dict:
    """
    Fit dose-response data to a 4-parameter logistic curve.

    Args:
        input_data: {
            "concentrations": [10000, 3333, 1111, ...],  # in order from high to low
            "responses": [95.2, 87.1, 62.3, ...],        # corresponding responses
            "controls": {"max": 100.0, "min": 2.3},      # optional control values
            "parameters": {                              # optional fitting parameters
                "fix_hill": null,      # fixed Hill coefficient (or null for free)
                "fix_top": null,       # fixed top asymptote (or null for free)
                "fix_bottom": null,    # fixed bottom asymptote (or null for free)
            }
        }

    Returns:
        {
            "ic50": 234.5,
            "hill_slope": 1.2,
            "top": 98.7,
            "bottom": 3.1,
            "r_squared": 0.994,
            "curve_points": [[x1, y1], ...],  # for plotting
            "flags": [],                       # warning flags
            "kpi": "ic50",                     # which result is the primary KPI
            "fit_successful": true,
        }
    """
    concentrations = np.array(input_data.get('concentrations', []))
    responses = np.array(input_data.get('responses', []))
    controls = input_data.get('controls', {})
    parameters = input_data.get('parameters', {})

    flags = []

    # Validate input
    if len(concentrations) < 4:
        return {
            'fit_successful': False,
            'flags': ['insufficient_data'],
            'kpi': 'ic50',
            'error': 'Need at least 4 data points for 4PL fitting',
        }

    if len(concentrations) != len(responses):
        return {
            'fit_successful': False,
            'flags': ['data_mismatch'],
            'kpi': 'ic50',
            'error': 'Concentrations and responses must have same length',
        }

    # Remove any NaN or infinite values
    valid_mask = np.isfinite(concentrations) & np.isfinite(responses) & (concentrations > 0)
    concentrations = concentrations[valid_mask]
    responses = responses[valid_mask]

    if len(concentrations) < 4:
        return {
            'fit_successful': False,
            'flags': ['insufficient_valid_data'],
            'kpi': 'ic50',
            'error': 'Less than 4 valid data points after filtering',
        }

    # Initial parameter estimates
    response_range = responses.max() - responses.min()

    # Use control values if available, otherwise estimate from data
    top_init = controls.get('max', responses.max())
    bottom_init = controls.get('min', responses.min())

    # Estimate IC50 as geometric mean of concentration range
    ic50_init = np.sqrt(concentrations.min() * concentrations.max())

    # Hill coefficient typically around 1
    hill_init = 1.0

    # Set up bounds and initial guess
    p0 = [top_init, bottom_init, ic50_init, hill_init]

    # Dynamic bounds based on actual data range (handles both percent and raw values)
    data_min = min(responses.min(), bottom_init if bottom_init else responses.min())
    data_max = max(responses.max(), top_init if top_init else responses.max())
    data_range = data_max - data_min

    # Allow some headroom beyond observed data
    bounds_lower = [data_min - data_range * 0.5, data_min - data_range * 0.5, concentrations.min() / 100, 0.1]
    bounds_upper = [data_max + data_range * 0.5, data_max + data_range * 0.5, concentrations.max() * 100, 10]

    # Handle fixed parameters
    fix_hill = parameters.get('fix_hill')
    fix_top = parameters.get('fix_top')
    fix_bottom = parameters.get('fix_bottom')

    try:
        if fix_hill is not None or fix_top is not None or fix_bottom is not None:
            # Create a modified function with fixed parameters
            def fit_func(x, *args):
                idx = 0
                top = fix_top if fix_top is not None else args[idx]
                if fix_top is None:
                    idx += 1
                bottom = fix_bottom if fix_bottom is not None else args[idx]
                if fix_bottom is None:
                    idx += 1
                ic50 = args[idx]
                idx += 1
                hill = fix_hill if fix_hill is not None else args[idx]
                return four_param_logistic(x, top, bottom, ic50, hill)

            # Adjust p0 and bounds for free parameters only
            p0_adj = []
            bounds_lower_adj = []
            bounds_upper_adj = []

            if fix_top is None:
                p0_adj.append(top_init)
                bounds_lower_adj.append(data_min - data_range * 0.5)
                bounds_upper_adj.append(data_max + data_range * 0.5)
            if fix_bottom is None:
                p0_adj.append(bottom_init)
                bounds_lower_adj.append(data_min - data_range * 0.5)
                bounds_upper_adj.append(data_max + data_range * 0.5)

            p0_adj.append(ic50_init)  # IC50 always free
            bounds_lower_adj.append(concentrations.min() / 100)
            bounds_upper_adj.append(concentrations.max() * 100)

            if fix_hill is None:
                p0_adj.append(hill_init)
                bounds_lower_adj.append(0.1)
                bounds_upper_adj.append(10)

            popt, pcov = curve_fit(
                fit_func, concentrations, responses,
                p0=p0_adj,
                bounds=(bounds_lower_adj, bounds_upper_adj),
                maxfev=5000
            )

            # Reconstruct full parameter set
            idx = 0
            top_fit = fix_top if fix_top is not None else popt[idx]
            if fix_top is None:
                idx += 1
            bottom_fit = fix_bottom if fix_bottom is not None else popt[idx]
            if fix_bottom is None:
                idx += 1
            ic50_fit = popt[idx]
            idx += 1
            hill_fit = fix_hill if fix_hill is not None else popt[idx]
        else:
            # Standard 4-parameter fit
            popt, pcov = curve_fit(
                four_param_logistic, concentrations, responses,
                p0=p0,
                bounds=(bounds_lower, bounds_upper),
                maxfev=5000
            )
            top_fit, bottom_fit, ic50_fit, hill_fit = popt

    except RuntimeError as e:
        flags.append('fit_failed')
        return {
            'fit_successful': False,
            'flags': flags,
            'kpi': 'ic50',
            'error': f'Curve fitting failed: {str(e)}',
        }
    except Exception as e:
        flags.append('fit_error')
        return {
            'fit_successful': False,
            'flags': flags,
            'kpi': 'ic50',
            'error': f'Unexpected error: {str(e)}',
        }

    # Calculate R-squared
    y_pred = four_param_logistic(concentrations, top_fit, bottom_fit, ic50_fit, hill_fit)
    ss_res = np.sum((responses - y_pred) ** 2)
    ss_tot = np.sum((responses - np.mean(responses)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    # Generate smooth curve for plotting
    x_smooth = np.logspace(
        np.log10(concentrations.min() / 2),
        np.log10(concentrations.max() * 2),
        100
    )
    y_smooth = four_param_logistic(x_smooth, top_fit, bottom_fit, ic50_fit, hill_fit)
    curve_points = [[float(x), float(y)] for x, y in zip(x_smooth, y_smooth)]

    # Quality checks
    if r_squared < 0.8:
        flags.append('poor_fit')

    if ic50_fit < concentrations.min() or ic50_fit > concentrations.max():
        flags.append('ic50_extrapolated')

    if abs(top_fit - bottom_fit) < 10:
        flags.append('low_dynamic_range')

    if hill_fit < 0.3 or hill_fit > 5:
        flags.append('unusual_hill_slope')

    # Check for incomplete curve
    response_at_highest = four_param_logistic(concentrations.max(), top_fit, bottom_fit, ic50_fit, hill_fit)
    response_at_lowest = four_param_logistic(concentrations.min(), top_fit, bottom_fit, ic50_fit, hill_fit)

    if abs(response_at_highest - top_fit) > 0.2 * abs(top_fit - bottom_fit):
        flags.append('incomplete_top')

    if abs(response_at_lowest - bottom_fit) > 0.2 * abs(top_fit - bottom_fit):
        flags.append('incomplete_bottom')

    return {
        'ic50': float(ic50_fit),
        'hill_slope': float(hill_fit),
        'top': float(top_fit),
        'bottom': float(bottom_fit),
        'r_squared': float(r_squared),
        'curve_points': curve_points,
        'flags': flags,
        'kpi': 'ic50',
        'fit_successful': True,
    }
