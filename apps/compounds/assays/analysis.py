"""
Assay Analysis Module

Runs curve fitting analysis on data series using FittingMethod scripts.
"""

import importlib.util
import logging
from typing import Optional

from .models import DataSeries, AnalysisResult, FittingMethod, Assay

logger = logging.getLogger(__name__)

# Default flags that cause automatic invalidation
# These are used when a protocol doesn't specify custom validation_rules
DEFAULT_INVALIDATING_FLAGS = ['poor_fit', 'insufficient_data', 'unusual_hill_slope']

# All available flags from fitting scripts (for reference and UI)
AVAILABLE_FLAGS = [
    {'id': 'poor_fit', 'label': 'Poor Fit', 'description': 'R² < 0.8'},
    {'id': 'insufficient_data', 'label': 'Insufficient Data', 'description': 'Fewer than 4 data points'},
    {'id': 'unusual_hill_slope', 'label': 'Unusual Hill Slope', 'description': 'Hill coefficient outside 0.3–5.0'},
    {'id': 'incomplete_top', 'label': 'Incomplete Top', 'description': 'Top asymptote not reached'},
    {'id': 'incomplete_bottom', 'label': 'Incomplete Bottom', 'description': 'Bottom asymptote not reached'},
    {'id': 'ic50_extrapolated', 'label': 'IC50 Extrapolated', 'description': 'IC50 outside measured concentration range'},
    {'id': 'low_dynamic_range', 'label': 'Low Dynamic Range', 'description': '|top - bottom| < 10'},
    {'id': 'ki_extrapolated', 'label': 'Ki Extrapolated', 'description': 'Ki outside reasonable range (Wang fitting)'},
    {'id': 'weak_binding_use_standard_analysis', 'label': 'Weak Binding', 'description': 'Ki > 10× protein concentration (Wang fitting)'},
]


def _create_or_update_analysis(data_series: 'DataSeries', status: str, results: dict) -> 'AnalysisResult':
    """Helper to create or update an AnalysisResult for a data series."""
    if data_series.analysis:
        analysis = data_series.analysis
        analysis.status = status
        analysis.results = results
        analysis.save()
    else:
        analysis = AnalysisResult.objects.create(
            status=status,
            results=results,
        )
        data_series.analysis = analysis
        data_series.save()
    return analysis


def get_builtin_fitting_function(method_slug: str = 'four-parameter-logistic'):
    """Load a built-in fitting function from fitting_scripts directory."""
    from . import fitting_scripts
    from .fitting_scripts import four_parameter_logistic
    return four_parameter_logistic.fit


def execute_fitting_script(script_code: str, input_data: dict) -> dict:
    """
    Execute a FittingMethod script on input data.

    Args:
        script_code: Python source code containing a fit() function
        input_data: Dictionary with concentrations, responses, controls, parameters

    Returns:
        Result dictionary from the fit() function
    """
    # Create a namespace with common imports available for fitting scripts
    # Include __builtins__ so import statements in scripts work
    import builtins
    import typing
    import numpy as np
    import scipy
    namespace = {
        '__builtins__': builtins,
        'np': np,
        'numpy': np,
        'scipy': scipy,
        'typing': typing,
        'Tuple': typing.Tuple,
        'Optional': typing.Optional,
    }
    # Add scipy.optimize for curve_fit
    try:
        from scipy.optimize import curve_fit, minimize
        namespace['curve_fit'] = curve_fit
        namespace['minimize'] = minimize
    except ImportError:
        pass

    # Execute the script to define the fit() function
    try:
        exec(script_code, namespace)
    except Exception as e:
        logger.error(f"Failed to execute fitting script: {e}")
        return {
            'fit_successful': False,
            'error': f'Script execution error: {str(e)}',
            'flags': ['script_error'],
            'kpi': 'ic50',
        }

    # Get the fit function
    fit_func = namespace.get('fit')
    if not callable(fit_func):
        return {
            'fit_successful': False,
            'error': 'Script does not define a fit() function',
            'flags': ['no_fit_function'],
            'kpi': 'ic50',
        }

    # Call the fit function
    try:
        result = fit_func(input_data)
        return result
    except Exception as e:
        logger.error(f"Fitting function raised exception: {e}")
        return {
            'fit_successful': False,
            'error': f'Fit execution error: {str(e)}',
            'flags': ['fit_exception'],
            'kpi': 'ic50',
        }


def analyse_data_series(data_series: DataSeries, fitting_method: Optional[FittingMethod] = None) -> AnalysisResult:
    """
    Run curve fitting analysis on a single data series.

    Args:
        data_series: DataSeries object with extracted_data
        fitting_method: FittingMethod to use (or None for default 4PL)

    Returns:
        AnalysisResult object (saved to database)
    """
    # Get concentrations from dilution series or assay protocol
    concentrations = []
    if data_series.dilution_series:
        concentrations = data_series.dilution_series.concentrations or []
    elif data_series.assay.protocol.preferred_dilutions:
        concentrations = data_series.assay.protocol.preferred_dilutions.concentrations or []

    # Get response data - format is [min_control, data1, ..., dataN, max_control]
    raw_data = data_series.extracted_data
    if not isinstance(raw_data, list):
        raw_data = []

    # Extract controls from first/last positions and responses from the middle
    # Data format: [high_signal_control, data1, ..., dataN, low_signal_control]
    # For inhibition assays:
    #   - First position = uninhibited control = HIGH signal → top asymptote for fitting
    #   - Last position = fully inhibited control = LOW signal → bottom asymptote for fitting
    controls = {}
    responses = []
    if len(raw_data) >= 3:
        # First element is high signal control (top of curve)
        if raw_data[0] is not None:
            controls['max'] = raw_data[0]
        # Last element is low signal control (bottom of curve)
        if raw_data[-1] is not None:
            controls['min'] = raw_data[-1]
        # Middle elements are the dose-response data
        responses = raw_data[1:-1]
    elif len(raw_data) > 0:
        # No controls embedded, treat all as responses
        responses = raw_data

    # Validate concentrations are available - do NOT auto-generate
    if not concentrations:
        error_msg = (
            f"No dilution series available for {data_series.compound_name}. "
            f"Set dilution_series on the data series or preferred_dilutions on the protocol."
        )
        logger.error(error_msg)
        # Create failed analysis result
        analysis = _create_or_update_analysis(
            data_series,
            status='invalid',
            results={
                'fit_successful': False,
                'error': error_msg,
                'flags': ['no_dilution_series'],
                'KPI': 'EC50',
            }
        )
        return analysis

    # Validate concentration/response length match
    if len(concentrations) != len(responses):
        error_msg = (
            f"Concentration/response length mismatch for {data_series.compound_name}: "
            f"{len(concentrations)} concentrations vs {len(responses)} responses."
        )
        logger.error(error_msg)
        analysis = _create_or_update_analysis(
            data_series,
            status='invalid',
            results={
                'fit_successful': False,
                'error': error_msg,
                'flags': ['data_length_mismatch'],
                'KPI': 'EC50',
            }
        )
        return analysis

    # Build input data for fitting
    input_data = {
        'concentrations': concentrations,
        'responses': responses,
        'controls': controls,
        'parameters': data_series.assay.protocol.fitting_parameters or {},
    }

    # Apply skip_points filter
    skip_points = data_series.skip_points or []
    if skip_points and concentrations and responses:
        filtered_conc = []
        filtered_resp = []
        for i, (c, r) in enumerate(zip(concentrations, responses)):
            if i not in skip_points:
                filtered_conc.append(c)
                filtered_resp.append(r)
        input_data['concentrations'] = filtered_conc
        input_data['responses'] = filtered_resp

    # Run fitting
    if fitting_method and fitting_method.script:
        result = execute_fitting_script(fitting_method.script, input_data)
    else:
        # Use built-in 4PL fitting
        try:
            fit_func = get_builtin_fitting_function()
            result = fit_func(input_data)
        except Exception as e:
            logger.error(f"Built-in fitting failed: {e}")
            result = {
                'fit_successful': False,
                'error': str(e),
                'flags': ['builtin_error'],
                'kpi': 'ic50',
            }

    # Determine status based on result and protocol validation rules
    # Get validation rules from protocol's fitting_parameters
    fitting_params = data_series.assay.protocol.fitting_parameters or {}
    validation_rules = fitting_params.get('validation_rules', {})
    invalidating_flags = validation_rules.get('invalidating_flags', DEFAULT_INVALIDATING_FLAGS)

    if result.get('fit_successful', False):
        flags = result.get('flags', [])
        # Mark as invalid if any flag is in the invalidating list
        if any(flag in invalidating_flags for flag in flags):
            status = 'invalid'
        else:
            status = 'valid'
    else:
        status = 'invalid'

    # Build results dict for storage
    # Concentrations come from dilution_series, not stored here
    kpi_value = result.get('kpi', 'ic50')
    kpi_key = kpi_value.upper() if kpi_value else 'IC50'

    # Get the primary KPI value from the fitting result
    # The fitting script returns lowercase keys (e.g., 'ic50', 'ki')
    primary_kpi_value = result.get(kpi_value) if kpi_value else result.get('ic50')

    # Determine algorithm identifier for curve drawing and parameter extraction
    # Use fitting_method slug if available, otherwise 'four-parameter-logistic' for built-in
    algorithm = fitting_method.slug if fitting_method else 'four-parameter-logistic'

    stored_results = {
        # Store the KPI value under the uppercased key that matches 'KPI' pointer
        kpi_key: primary_kpi_value,
        'IC50_apparent': result.get('ic50_apparent'),
        'Hill': result.get('hill_slope'),
        'minVal': result.get('bottom'),
        'maxVal': result.get('top'),
        'r_squared': result.get('r_squared'),
        'KPI': kpi_key,
        'flags': result.get('flags', []),
        'curve_points': result.get('curve_points'),
        'fit_successful': result.get('fit_successful', False),
        'error': result.get('error'),
        'tight_binding_params': result.get('tight_binding_params'),
        'algorithm': algorithm,  # Identifies fitting algorithm for curve drawing
        # End percent: percentage of inhibition at highest concentration
        'end_percent': result.get('end_percent'),
        'end_percent_display': result.get('end_percent_display'),
    }

    return _create_or_update_analysis(data_series, status, stored_results)


def analyse_assay(assay: Assay) -> dict:
    """
    Run analysis on all data series in an assay.

    Args:
        assay: Assay object to analyse

    Returns:
        Summary dict with success/failure counts
    """
    # Get fitting method from protocol
    fitting_method = assay.protocol.get_effective_fitting_method()

    data_series_list = assay.data_series.all()
    results = {
        'total': data_series_list.count(),
        'successful': 0,
        'failed': 0,
        'valid': 0,
        'invalid': 0,
        'details': [],
    }

    for ds in data_series_list:
        try:
            analysis = analyse_data_series(ds, fitting_method)
            results['successful'] += 1
            if analysis.status == 'valid':
                results['valid'] += 1
            else:
                results['invalid'] += 1
            results['details'].append({
                'data_series_id': str(ds.id),
                'compound_name': ds.compound_name,
                'status': analysis.status,
                'kpi_value': analysis.kpi_value,
            })
        except Exception as e:
            logger.error(f"Failed to analyse data series {ds.id}: {e}")
            results['failed'] += 1
            results['details'].append({
                'data_series_id': str(ds.id),
                'compound_name': ds.compound_name,
                'error': str(e),
            })

    return results

