"""
Tight-Binding Competition Analysis using the Wang Equation

For competitive binding assays where inhibitor concentrations approach
or exceed the protein concentration, standard IC50 analysis is invalid.
This script implements the Wang equation to calculate true Ki values.

Reference:
Wang, Z-X. (1995) "An exact mathematical expression for describing
competitive binding of two different ligands to a protein molecule"
FEBS Letters 360, 111-114.

Script interface for the FittingMethod model.
"""

import numpy as np
from scipy.optimize import curve_fit
from typing import Tuple, Optional


def solve_wang_cubic(
    P_total: float,
    L_total: float,
    I_total: float,
    Kd_L: float,
    Ki: float
) -> Tuple[float, float]:
    """
    Solve the Wang cubic equation for free ligand and inhibitor concentrations.

    The competitive binding equilibrium leads to a cubic equation in [L]_free.
    This function finds the physically meaningful root.

    Args:
        P_total: Total protein concentration
        L_total: Total labeled ligand concentration
        I_total: Total inhibitor concentration
        Kd_L: Dissociation constant of labeled ligand
        Ki: Inhibition constant of inhibitor

    Returns:
        Tuple of (L_free, I_free) concentrations
    """
    # Handle edge case of no inhibitor
    if I_total <= 0:
        # Simple binding equation: [L]_free from quadratic
        # [P]_total * [L]_free / (Kd_L + [L]_free) + [L]_free = [L]_total
        # Rearranging: [L]_free^2 + (Kd_L + P_total - L_total)*[L]_free - Kd_L*L_total = 0
        a = 1
        b = Kd_L + P_total - L_total
        c = -Kd_L * L_total
        discriminant = b**2 - 4*a*c
        L_free = (-b + np.sqrt(discriminant)) / (2*a)
        return max(0, min(L_free, L_total)), 0.0

    # Coefficients for cubic in [L]_free derived from Wang equation
    # The competitive binding mass balance equations lead to:
    # [L]^3 + a[L]^2 + b[L] + c = 0
    a = Kd_L + Ki + L_total + I_total - P_total
    b = Ki * (L_total - P_total) + Kd_L * (I_total - P_total) + Kd_L * Ki
    c = -Kd_L * Ki * L_total

    # Solve cubic using numpy.roots (returns all 3 roots)
    coeffs = [1, a, b, c]
    roots = np.roots(coeffs)

    # Find the physically meaningful root (real, positive, <= L_total)
    L_free = None
    valid_roots = []
    for root in roots:
        if np.isreal(root):
            val = np.real(root)
            if 0 < val <= L_total * 1.001:  # Small tolerance for numerical error
                valid_roots.append(val)

    if valid_roots:
        # Take the smallest valid root (most physical)
        L_free = min(valid_roots)
    else:
        # Fallback: use smallest positive real root
        real_positive = [np.real(r) for r in roots if np.isreal(r) and np.real(r) > 0]
        if real_positive:
            L_free = min(real_positive)
        else:
            # Last resort fallback
            L_free = L_total * 0.5

    L_free = max(0, min(L_free, L_total))

    # Calculate I_free from mass balance
    # [PL] = P_total * L_free / (Kd_L + L_free * (1 + I_free/Ki))
    # This requires solving another equation, but we can approximate:
    # At equilibrium, I_free = I_total - [PI]
    # [PI] = P_total * I_free / (Ki + I_free * (1 + L_free/Kd_L))

    # Simplified: assume most inhibitor is free at moderate concentrations
    # More accurate calculation using mass balance
    PL = P_total * L_free / (Kd_L + L_free) if (Kd_L + L_free) > 0 else 0
    P_free = P_total - PL

    # From I equilibrium: [PI] = P_free * I_total / (Ki + I_total) approximately
    # This is simplified; full solution would need another iteration
    if Ki + I_total > 0:
        PI = P_free * I_total / (Ki + I_total)
        I_free = I_total - PI
    else:
        I_free = I_total

    return L_free, max(0, I_free)


def fractional_binding(
    I_total: float,
    P_total: float,
    L_total: float,
    Kd_L: float,
    Ki: float
) -> float:
    """
    Calculate fractional binding of labeled ligand at given inhibitor concentration.

    Returns Y = [PL] / [PL]_0 where [PL]_0 is binding without inhibitor.
    """
    # Baseline binding (no inhibitor)
    L_free_0, _ = solve_wang_cubic(P_total, L_total, 0, Kd_L, Ki)
    PL_0 = P_total * L_free_0 / (Kd_L + L_free_0) if (Kd_L + L_free_0) > 0 else P_total

    if I_total <= 0:
        return 1.0

    # With inhibitor
    L_free, _ = solve_wang_cubic(P_total, L_total, I_total, Kd_L, Ki)
    PL = P_total * L_free / (Kd_L + L_free) if (Kd_L + L_free) > 0 else 0

    return PL / PL_0 if PL_0 > 0 else 0


def wang_response_curve(
    I_concentrations: np.ndarray,
    Ki: float,
    P_total: float,
    L_total: float,
    Kd_L: float,
    top: float,
    bottom: float
) -> np.ndarray:
    """
    Calculate expected response at each inhibitor concentration using Wang model.

    Args:
        I_concentrations: Array of inhibitor concentrations
        Ki: Inhibition constant
        P_total: Total protein concentration
        L_total: Total labeled ligand concentration
        Kd_L: Dissociation constant of labeled ligand
        top: Response at zero inhibitor (maximum signal)
        bottom: Response at saturating inhibitor (minimum signal)

    Returns:
        Array of response values
    """
    responses = []
    for I_total in I_concentrations:
        frac = fractional_binding(float(I_total), P_total, L_total, Kd_L, Ki)
        # Response scales with fractional binding
        response = bottom + (top - bottom) * frac
        responses.append(response)
    return np.array(responses)


def fit(input_data: dict) -> dict:
    """
    Fit dose-response data using the Wang tight-binding equation.

    Args:
        input_data: {
            "concentrations": [10000, 3333, ...],  # Inhibitor concentrations (nM)
            "responses": [95.2, 87.1, ...],         # Corresponding responses
            "controls": {"max": 100.0, "min": 2.3}, # Control values
            "parameters": {
                "protein_conc": 50.0,    # [P]_total in nM - REQUIRED
                "ligand_conc": 10.0,     # [L]_total in nM - REQUIRED
                "ligand_kd": 5.0,        # Kd of labeled ligand in nM - REQUIRED
                "fix_top": null,         # Optional: fix top asymptote
                "fix_bottom": null,      # Optional: fix bottom asymptote
            }
        }

    Returns:
        {
            "ki": 15.3,              # Fitted Ki value (primary KPI)
            "ic50_apparent": 234.5,  # Apparent IC50 for reference
            "top": 98.7,
            "bottom": 3.1,
            "r_squared": 0.994,
            "curve_points": [[x1, y1], ...],
            "flags": [],
            "kpi": "ki",
            "fit_successful": true,
            "tight_binding_params": {
                "protein_conc": 50.0,
                "ligand_conc": 10.0,
                "ligand_kd": 5.0,
            }
        }
    """
    concentrations = np.array(input_data.get('concentrations', []))
    responses = np.array(input_data.get('responses', []))
    controls = input_data.get('controls', {})
    parameters = input_data.get('parameters', {})

    flags = []

    # Validate required tight-binding parameters
    P_total = parameters.get('protein_conc')
    L_total = parameters.get('ligand_conc')
    Kd_L = parameters.get('ligand_kd')

    missing_params = []
    if P_total is None:
        missing_params.append('protein_conc')
    if L_total is None:
        missing_params.append('ligand_conc')
    if Kd_L is None:
        missing_params.append('ligand_kd')

    if missing_params:
        return {
            'fit_successful': False,
            'flags': ['missing_tight_binding_params'],
            'kpi': 'ki',
            'error': f'Missing required parameters: {", ".join(missing_params)}. '
                     f'Configure protein_conc, ligand_conc, and ligand_kd in protocol fitting_parameters.',
        }

    # Convert to float and validate
    try:
        P_total = float(P_total)
        L_total = float(L_total)
        Kd_L = float(Kd_L)
    except (TypeError, ValueError) as e:
        return {
            'fit_successful': False,
            'flags': ['invalid_tight_binding_params'],
            'kpi': 'ki',
            'error': f'Invalid parameter values: {e}',
        }

    if P_total <= 0 or L_total <= 0 or Kd_L <= 0:
        return {
            'fit_successful': False,
            'flags': ['invalid_tight_binding_params'],
            'kpi': 'ki',
            'error': 'protein_conc, ligand_conc, and ligand_kd must be positive values',
        }

    # Validate input data
    if len(concentrations) < 4:
        return {
            'fit_successful': False,
            'flags': ['insufficient_data'],
            'kpi': 'ki',
            'error': 'Need at least 4 data points for fitting',
        }

    if len(concentrations) != len(responses):
        return {
            'fit_successful': False,
            'flags': ['data_mismatch'],
            'kpi': 'ki',
            'error': 'Concentrations and responses must have same length',
        }

    # Remove invalid values
    valid_mask = np.isfinite(concentrations) & np.isfinite(responses) & (concentrations > 0)
    concentrations = concentrations[valid_mask]
    responses = responses[valid_mask]

    if len(concentrations) < 4:
        return {
            'fit_successful': False,
            'flags': ['insufficient_valid_data'],
            'kpi': 'ki',
            'error': 'Less than 4 valid data points after filtering',
        }

    # Initial estimates
    resp_min = float(responses.min())
    resp_max = float(responses.max())

    # Get control values - handle potentially mislabeled controls
    # For inhibition assays: 'max' should be high signal (uninhibited), 'min' should be low signal (inhibited)
    ctrl_max = controls.get('max')
    ctrl_min = controls.get('min')

    # If controls are provided but inverted (max < min), swap them
    if ctrl_max is not None and ctrl_min is not None and ctrl_max < ctrl_min:
        ctrl_max, ctrl_min = ctrl_min, ctrl_max

    # Use controls if available, otherwise use data extremes
    top_init = ctrl_max if ctrl_max is not None else resp_max
    bottom_init = ctrl_min if ctrl_min is not None else resp_min

    # Final safety check: ensure top > bottom
    if top_init < bottom_init:
        top_init, bottom_init = bottom_init, top_init

    # Estimate initial Ki from apparent IC50
    # Rough IC50 estimate at midpoint response
    mid_response = (top_init + bottom_init) / 2
    idx_mid = np.argmin(np.abs(responses - mid_response))
    ki_init = float(concentrations[idx_mid])

    # Set bounds with generous margins
    conc_min = float(concentrations.min())
    conc_max = float(concentrations.max())

    # Ki bounds: very wide range from very low to very high
    ki_lower = conc_min / 10000
    ki_upper = conc_max * 10000

    # Clamp ki_init to be strictly within bounds
    ki_init = max(ki_lower * 2, min(ki_init, ki_upper / 2))

    # Response bounds: allow fitting beyond observed range
    resp_range = max(resp_max - resp_min, 10)  # Minimum range of 10

    # Top should be at or above max response
    top_lower = resp_min - resp_range
    top_upper = resp_max + resp_range * 2

    # Bottom should be at or below min response
    bottom_lower = resp_min - resp_range * 2
    bottom_upper = resp_max + resp_range

    # Clamp initial values to be strictly within bounds
    top_init = max(top_lower + 1, min(float(top_init), top_upper - 1))
    bottom_init = max(bottom_lower + 1, min(float(bottom_init), bottom_upper - 1))

    # Ensure top_init > bottom_init after clamping
    if top_init <= bottom_init:
        mid = (top_init + bottom_init) / 2
        top_init = mid + 5
        bottom_init = mid - 5

    # Handle fixed parameters
    fix_top = parameters.get('fix_top')
    fix_bottom = parameters.get('fix_bottom')

    # Initialize for error reporting
    p0 = None
    bounds_lower = None
    bounds_upper = None

    # Pre-validate that all bounds are properly ordered
    def validate_and_fix_bounds(p0_list, lower_list, upper_list, names):
        """Ensure bounds are valid and p0 is strictly within them."""
        fixed_p0 = []
        fixed_lower = []
        fixed_upper = []

        for i, (p, lo, hi, name) in enumerate(zip(p0_list, lower_list, upper_list, names)):
            # Ensure lower < upper
            if lo >= hi:
                # Widen the bounds
                mid = (lo + hi) / 2
                spread = max(abs(mid) * 0.5, 100)
                lo = mid - spread
                hi = mid + spread

            # Ensure p is strictly within [lo, hi]
            margin = (hi - lo) * 0.01  # 1% margin
            if p <= lo:
                p = lo + margin
            elif p >= hi:
                p = hi - margin

            fixed_p0.append(float(p))
            fixed_lower.append(float(lo))
            fixed_upper.append(float(hi))

        return fixed_p0, fixed_lower, fixed_upper

    try:
        if fix_top is not None or fix_bottom is not None:
            # Fit with some parameters fixed
            def fit_func(I, *args):
                idx = 0
                ki = args[idx]
                idx += 1
                top = float(fix_top) if fix_top is not None else args[idx]
                if fix_top is None:
                    idx += 1
                bottom = float(fix_bottom) if fix_bottom is not None else args[idx]
                return wang_response_curve(I, ki, P_total, L_total, Kd_L, top, bottom)

            # Build p0 and bounds
            p0 = [ki_init]
            bounds_lower = [ki_lower]
            bounds_upper = [ki_upper]

            if fix_top is None:
                p0.append(top_init)
                bounds_lower.append(top_lower)
                bounds_upper.append(top_upper)
            if fix_bottom is None:
                p0.append(bottom_init)
                bounds_lower.append(bottom_lower)
                bounds_upper.append(bottom_upper)

            # Validate bounds before fitting
            param_names = ['ki'] + (['top'] if fix_top is None else []) + (['bottom'] if fix_bottom is None else [])
            p0, bounds_lower, bounds_upper = validate_and_fix_bounds(
                p0, bounds_lower, bounds_upper, param_names
            )

            popt, pcov = curve_fit(
                fit_func, concentrations, responses,
                p0=p0,
                bounds=(bounds_lower, bounds_upper),
                maxfev=10000
            )

            # Extract results
            idx = 0
            ki_fit = popt[idx]
            idx += 1
            top_fit = float(fix_top) if fix_top is not None else popt[idx]
            if fix_top is None:
                idx += 1
            bottom_fit = float(fix_bottom) if fix_bottom is not None else popt[idx]
        else:
            # Full 3-parameter fit: Ki, top, bottom
            def fit_func(I, ki, top, bottom):
                return wang_response_curve(I, ki, P_total, L_total, Kd_L, top, bottom)

            p0 = [ki_init, top_init, bottom_init]
            bounds_lower = [ki_lower, top_lower, bottom_lower]
            bounds_upper = [ki_upper, top_upper, bottom_upper]

            # Validate bounds before fitting
            p0, bounds_lower, bounds_upper = validate_and_fix_bounds(
                p0, bounds_lower, bounds_upper, ['ki', 'top', 'bottom']
            )

            popt, pcov = curve_fit(
                fit_func, concentrations, responses,
                p0=p0,
                bounds=(bounds_lower, bounds_upper),
                maxfev=10000
            )
            ki_fit, top_fit, bottom_fit = popt

    except RuntimeError as e:
        flags.append('fit_failed')
        return {
            'fit_successful': False,
            'flags': flags,
            'kpi': 'ki',
            'error': f'Curve fitting failed: {str(e)}',
            'debug': {
                'p0': p0,
                'bounds_lower': bounds_lower,
                'bounds_upper': bounds_upper,
                'conc_range': [conc_min, conc_max],
                'resp_range': [resp_min, resp_max],
            }
        }
    except Exception as e:
        flags.append('fit_error')
        return {
            'fit_successful': False,
            'flags': flags,
            'kpi': 'ki',
            'error': f'Unexpected error: {str(e)}. For Wang fit with [P]={P_total} (nM), Kd={Kd_L} (nM), [L]={L_total} (nM)',
            'debug': {
                'p0': p0 if 'p0' in dir() else None,
                'bounds_lower': bounds_lower if 'bounds_lower' in dir() else None,
                'bounds_upper': bounds_upper if 'bounds_upper' in dir() else None,
                'conc_range': [conc_min, conc_max],
                'resp_range': [resp_min, resp_max],
            }
        }

    # Calculate R-squared
    y_pred = wang_response_curve(concentrations, ki_fit, P_total, L_total, Kd_L, top_fit, bottom_fit)
    ss_res = np.sum((responses - y_pred) ** 2)
    ss_tot = np.sum((responses - np.mean(responses)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    # Calculate apparent IC50 for reference (concentration at 50% effect)
    ic50_apparent = None
    try:
        target_response = (top_fit + bottom_fit) / 2
        # Binary search for IC50
        low = float(concentrations.min()) / 100
        high = float(concentrations.max()) * 100
        for _ in range(50):
            mid = np.sqrt(low * high)
            resp = wang_response_curve(
                np.array([mid]), ki_fit, P_total, L_total, Kd_L, top_fit, bottom_fit
            )[0]
            if resp > target_response:
                low = mid
            else:
                high = mid
        ic50_apparent = float(np.sqrt(low * high))
    except Exception:
        pass

    # Generate curve points for plotting
    x_smooth = np.logspace(
        np.log10(float(concentrations.min()) / 2),
        np.log10(float(concentrations.max()) * 2),
        100
    )
    y_smooth = wang_response_curve(x_smooth, ki_fit, P_total, L_total, Kd_L, top_fit, bottom_fit)
    curve_points = [[float(x), float(y)] for x, y in zip(x_smooth, y_smooth)]

    # Quality checks
    if r_squared < 0.8:
        flags.append('poor_fit')

    if ki_fit < float(concentrations.min()) / 10 or ki_fit > float(concentrations.max()) * 10:
        flags.append('ki_extrapolated')

    if abs(top_fit - bottom_fit) < 10:
        flags.append('low_dynamic_range')

    # Check tight-binding condition
    if ki_fit > P_total * 10:
        flags.append('weak_binding_use_standard_analysis')

    return {
        'ki': float(ki_fit),
        'ic50_apparent': ic50_apparent,
        'top': float(top_fit),
        'bottom': float(bottom_fit),
        'r_squared': float(r_squared),
        'curve_points': curve_points,
        'flags': flags,
        'kpi': 'ki',
        'fit_successful': True,
        'tight_binding_params': {
            'protein_conc': P_total,
            'ligand_conc': L_total,
            'ligand_kd': Kd_L,
        }
    }
