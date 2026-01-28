"""
Quality Control (QC) metrics for assay plates.

Calculates standard HTS quality metrics from control wells:
- Z-prime (Z') and Robust Z-prime
- Signal-to-Background (S/B) and Signal-to-Noise (S/N)
- Coefficient of Variation (CV%) for controls
- SSMD (Strictly Standardized Mean Difference)

Also aggregates curve fitting quality from AnalysisResult data.
"""

from dataclasses import dataclass
from statistics import mean, stdev, median
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .models import Assay


def _median_absolute_deviation(values: list[float]) -> float:
    """Calculate Median Absolute Deviation (MAD)."""
    med = median(values)
    return median([abs(v - med) for v in values])


@dataclass
class QCMetrics:
    """Container for assay QC metrics."""

    # Plate-level metrics (from controls)
    z_prime: float | None = None
    robust_z_prime: float | None = None
    signal_to_background: float | None = None
    signal_to_noise: float | None = None
    ssmd: float | None = None
    cv_high_controls: float | None = None
    cv_low_controls: float | None = None

    # Control statistics
    high_control_mean: float | None = None
    high_control_stdev: float | None = None
    low_control_mean: float | None = None
    low_control_stdev: float | None = None
    n_high_controls: int = 0
    n_low_controls: int = 0

    # Curve fitting quality (from AnalysisResult)
    total_curves: int = 0
    valid_curves: int = 0
    invalid_curves: int = 0
    unassigned_curves: int = 0
    curves_with_good_fit: int = 0  # R² > 0.8
    curves_with_flags: int = 0

    # Overall assessment
    overall_quality: str = 'unknown'  # excellent, good, acceptable, poor, insufficient_data
    issues: list[str] | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            # Plate-level metrics
            'z_prime': self.z_prime,
            'robust_z_prime': self.robust_z_prime,
            'signal_to_background': self.signal_to_background,
            'signal_to_noise': self.signal_to_noise,
            'ssmd': self.ssmd,
            'cv_high_controls': self.cv_high_controls,
            'cv_low_controls': self.cv_low_controls,

            # Control statistics
            'high_control_mean': self.high_control_mean,
            'high_control_stdev': self.high_control_stdev,
            'low_control_mean': self.low_control_mean,
            'low_control_stdev': self.low_control_stdev,
            'n_high_controls': self.n_high_controls,
            'n_low_controls': self.n_low_controls,

            # Curve fitting quality
            'total_curves': self.total_curves,
            'valid_curves': self.valid_curves,
            'invalid_curves': self.invalid_curves,
            'unassigned_curves': self.unassigned_curves,
            'curves_with_good_fit': self.curves_with_good_fit,
            'curves_with_flags': self.curves_with_flags,

            # Derived percentages
            'percent_valid': (
                round(100 * self.valid_curves / self.total_curves, 1)
                if self.total_curves > 0 else None
            ),
            'percent_good_fit': (
                round(100 * self.curves_with_good_fit / self.total_curves, 1)
                if self.total_curves > 0 else None
            ),

            # Overall assessment
            'overall_quality': self.overall_quality,
            'issues': self.issues or [],
        }


def calculate_qc_metrics(assay: 'Assay') -> QCMetrics:
    """
    Calculate comprehensive QC metrics for an assay.

    Extracts control values from all data series and computes standard
    HTS quality metrics. Also aggregates curve fitting statistics.

    Args:
        assay: The Assay instance to analyze

    Returns:
        QCMetrics dataclass with all calculated values
    """
    metrics = QCMetrics()
    issues = []

    # Collect control values and curve stats from all data series
    high_controls = []
    low_controls = []

    data_series_qs = assay.data_series.select_related('analysis').all()

    for ds in data_series_qs:
        # Extract controls from extracted_data
        # Format: [high_signal_control, data1, ..., dataN, low_signal_control]
        data = ds.extracted_data
        if data and isinstance(data, list) and len(data) >= 3:
            high_val = data[0]
            low_val = data[-1]

            # Only include numeric values
            if isinstance(high_val, (int, float)) and high_val is not None:
                high_controls.append(float(high_val))
            if isinstance(low_val, (int, float)) and low_val is not None:
                low_controls.append(float(low_val))

        # Aggregate curve fitting stats
        metrics.total_curves += 1

        if ds.analysis:
            status = ds.analysis.status
            if status == 'valid':
                metrics.valid_curves += 1
            elif status == 'invalid':
                metrics.invalid_curves += 1
            else:
                metrics.unassigned_curves += 1

            results = ds.analysis.results or {}

            # Check R² for good fit
            r_squared = results.get('r_squared')
            if r_squared is not None and r_squared > 0.8:
                metrics.curves_with_good_fit += 1

            # Check for flags
            flags = results.get('flags', [])
            if flags:
                metrics.curves_with_flags += 1
        else:
            metrics.unassigned_curves += 1

    # Store control counts
    metrics.n_high_controls = len(high_controls)
    metrics.n_low_controls = len(low_controls)

    # Calculate plate-level metrics if we have enough controls
    if len(high_controls) >= 3 and len(low_controls) >= 3:
        # Basic statistics
        mean_high = mean(high_controls)
        mean_low = mean(low_controls)
        stdev_high = stdev(high_controls)
        stdev_low = stdev(low_controls)

        metrics.high_control_mean = round(mean_high, 3)
        metrics.high_control_stdev = round(stdev_high, 3)
        metrics.low_control_mean = round(mean_low, 3)
        metrics.low_control_stdev = round(stdev_low, 3)

        # Coefficient of Variation (CV%)
        if mean_high != 0:
            metrics.cv_high_controls = round(100 * stdev_high / abs(mean_high), 2)
        if mean_low != 0:
            metrics.cv_low_controls = round(100 * stdev_low / abs(mean_low), 2)

        # Signal separation
        signal_diff = abs(mean_high - mean_low)

        if signal_diff > 0:
            # Z-prime: Z' = 1 - 3(σ_h + σ_l) / |μ_h - μ_l|
            metrics.z_prime = round(
                1 - (3 * (stdev_high + stdev_low) / signal_diff),
                3
            )

            # Robust Z-prime using median and MAD
            # Scale factor 1.4826 converts MAD to stdev equivalent for normal distribution
            median_high = median(high_controls)
            median_low = median(low_controls)
            mad_high = _median_absolute_deviation(high_controls)
            mad_low = _median_absolute_deviation(low_controls)
            robust_signal_diff = abs(median_high - median_low)

            if robust_signal_diff > 0:
                metrics.robust_z_prime = round(
                    1 - (3 * 1.4826 * (mad_high + mad_low) / robust_signal_diff),
                    3
                )

            # Signal-to-Background (S/B)
            if mean_low != 0:
                metrics.signal_to_background = round(mean_high / mean_low, 2)

            # Signal-to-Noise (S/N): (μ_h - μ_l) / σ_l
            if stdev_low > 0:
                metrics.signal_to_noise = round(signal_diff / stdev_low, 2)

            # SSMD: (μ_h - μ_l) / √(σ_h² + σ_l²)
            combined_var = stdev_high**2 + stdev_low**2
            if combined_var > 0:
                metrics.ssmd = round(signal_diff / (combined_var ** 0.5), 2)
    else:
        if len(high_controls) < 3:
            issues.append(f"Insufficient high controls ({len(high_controls)} < 3)")
        if len(low_controls) < 3:
            issues.append(f"Insufficient low controls ({len(low_controls)} < 3)")

    # Flag CV issues
    if metrics.cv_high_controls is not None and metrics.cv_high_controls > 20:
        issues.append(f"High control CV is elevated ({metrics.cv_high_controls}% > 20%)")
    if metrics.cv_low_controls is not None and metrics.cv_low_controls > 20:
        issues.append(f"Low control CV is elevated ({metrics.cv_low_controls}% > 20%)")

    # Flag low valid curve rate
    if metrics.total_curves > 0:
        valid_pct = 100 * metrics.valid_curves / metrics.total_curves
        if valid_pct < 50:
            issues.append(f"Low valid curve rate ({valid_pct:.0f}%)")

    # Determine overall quality based on Z-prime
    if metrics.z_prime is not None:
        if metrics.z_prime >= 0.5:
            metrics.overall_quality = 'excellent'
        elif metrics.z_prime >= 0.3:
            metrics.overall_quality = 'good'
        elif metrics.z_prime >= 0:
            metrics.overall_quality = 'acceptable'
        else:
            metrics.overall_quality = 'poor'
            issues.append(f"Z-prime is negative ({metrics.z_prime})")
    elif metrics.total_curves > 0:
        # No plate controls but we have curve data
        metrics.overall_quality = 'no_plate_controls'
    else:
        metrics.overall_quality = 'insufficient_data'

    metrics.issues = issues
    return metrics
