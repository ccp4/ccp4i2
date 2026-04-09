'use client';

import {
  Box,
  Typography,
  Grid2 as Grid,
  Skeleton,
  Alert,
  Chip,
  Accordion,
  AccordionSummary,
  AccordionDetails,
} from '@mui/material';
import {
  CheckCircleOutline,
  WarningAmber,
  ErrorOutline,
  InfoOutlined,
  TrendingUp,
  Science,
  ExpandMore,
} from '@mui/icons-material';
import { useCompoundsApi } from '@/lib/compounds/api';

interface QCMetrics {
  // Plate-level metrics
  z_prime: number | null;
  robust_z_prime: number | null;
  signal_to_background: number | null;
  signal_to_noise: number | null;
  ssmd: number | null;
  cv_high_controls: number | null;
  cv_low_controls: number | null;

  // Control statistics
  high_control_mean: number | null;
  high_control_stdev: number | null;
  low_control_mean: number | null;
  low_control_stdev: number | null;
  n_high_controls: number;
  n_low_controls: number;

  // Curve fitting quality
  total_curves: number;
  valid_curves: number;
  invalid_curves: number;
  unassigned_curves: number;
  curves_with_good_fit: number;
  curves_with_flags: number;

  // Derived
  percent_valid: number | null;
  percent_good_fit: number | null;

  // Overall
  overall_quality: string;
  issues: string[];
}

interface QCPanelProps {
  assayId: string;
}

// RAG (Red-Amber-Green) color scheme for QC metrics
const RAG_COLORS = {
  green: { bg: '#e8f5e9', border: '#4caf50', text: '#2e7d32' },   // Good
  amber: { bg: '#fff3e0', border: '#ff9800', text: '#e65100' },   // Warning
  red: { bg: '#ffebee', border: '#f44336', text: '#c62828' },     // Poor
  neutral: { bg: '#f5f5f5', border: '#9e9e9e', text: '#424242' }, // No data
};

function MetricCard({
  label,
  value,
  description,
  quality,
  suffix,
}: {
  label: string;
  value: string | number | null;
  description?: string;
  quality?: 'good' | 'warning' | 'error' | 'neutral';
  suffix?: string;
}) {
  const getColors = () => {
    switch (quality) {
      case 'good':
        return RAG_COLORS.green;
      case 'warning':
        return RAG_COLORS.amber;
      case 'error':
        return RAG_COLORS.red;
      default:
        return RAG_COLORS.neutral;
    }
  };

  const colors = getColors();

  return (
    <Box
      sx={{
        p: 1.5,
        bgcolor: colors.bg,
        borderRadius: 1,
        borderLeft: `4px solid ${colors.border}`,
        minWidth: 100,
        height: '100%',
      }}
    >
      <Typography
        variant="subtitle2"
        fontWeight={600}
        sx={{ color: 'text.primary', mb: 0.5 }}
      >
        {label}
      </Typography>
      <Typography
        variant="h5"
        fontWeight={700}
        fontFamily="monospace"
        sx={{ color: colors.text, lineHeight: 1.2 }}
      >
        {value !== null && value !== undefined ? (
          <>
            {typeof value === 'number' ? value.toFixed(2) : value}
            {suffix && (
              <Typography
                component="span"
                variant="body2"
                sx={{ color: colors.text, opacity: 0.8, ml: 0.25 }}
              >
                {suffix}
              </Typography>
            )}
          </>
        ) : (
          '-'
        )}
      </Typography>
      {description && (
        <Typography
          variant="caption"
          sx={{ color: 'text.secondary', display: 'block', mt: 0.5, lineHeight: 1.2 }}
        >
          {description}
        </Typography>
      )}
    </Box>
  );
}

function QualityBadge({ quality }: { quality: string }) {
  const config: Record<string, { color: 'success' | 'info' | 'warning' | 'error' | 'default'; icon: React.ReactNode; label: string }> = {
    excellent: { color: 'success', icon: <CheckCircleOutline fontSize="small" />, label: 'Excellent' },
    good: { color: 'success', icon: <CheckCircleOutline fontSize="small" />, label: 'Good' },
    acceptable: { color: 'info', icon: <InfoOutlined fontSize="small" />, label: 'Acceptable' },
    poor: { color: 'error', icon: <ErrorOutline fontSize="small" />, label: 'Poor' },
    no_plate_controls: { color: 'warning', icon: <WarningAmber fontSize="small" />, label: 'No Plate Controls' },
    insufficient_data: { color: 'default', icon: <InfoOutlined fontSize="small" />, label: 'Insufficient Data' },
  };

  const { color, icon, label } = config[quality] || config.insufficient_data;

  return (
    <Chip
      icon={icon as React.ReactElement}
      label={label}
      color={color}
      size="small"
      sx={{ fontWeight: 500 }}
    />
  );
}

function getZPrimeQuality(value: number | null): 'good' | 'warning' | 'error' | 'neutral' {
  if (value === null) return 'neutral';
  if (value >= 0.5) return 'good';
  if (value >= 0) return 'warning';
  return 'error';
}

function getCVQuality(value: number | null): 'good' | 'warning' | 'error' | 'neutral' {
  if (value === null) return 'neutral';
  if (value <= 10) return 'good';
  if (value <= 20) return 'warning';
  return 'error';
}

function getPercentQuality(value: number | null, threshold: number = 70): 'good' | 'warning' | 'error' | 'neutral' {
  if (value === null) return 'neutral';
  if (value >= threshold) return 'good';
  if (value >= 50) return 'warning';
  return 'error';
}

export function QCPanel({ assayId }: QCPanelProps) {
  const api = useCompoundsApi();
  const { data: metrics, isLoading, error } = api.get<QCMetrics>(`assays/${assayId}/qc/`);

  if (isLoading) {
    return (
      <Accordion sx={{ mb: 3 }} disabled>
        <AccordionSummary>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
            <Skeleton variant="circular" width={24} height={24} />
            <Skeleton variant="text" width={120} />
            <Skeleton variant="rounded" width={80} height={24} sx={{ ml: 2 }} />
          </Box>
        </AccordionSummary>
      </Accordion>
    );
  }

  if (error || !metrics) {
    return null; // Silently fail - QC is supplementary info
  }

  const hasPlateControls = metrics.n_high_controls >= 3 && metrics.n_low_controls >= 3;

  return (
    <Accordion defaultExpanded={false} sx={{ mb: 3 }}>
      <AccordionSummary
        expandIcon={<ExpandMore />}
        sx={{ '& .MuiAccordionSummary-content': { alignItems: 'center' } }}
      >
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, flex: 1 }}>
          <Science sx={{ color: 'info.main' }} />
          <Typography variant="h6">Quality Control</Typography>
          <Box sx={{ ml: 2 }}>
            <QualityBadge quality={metrics.overall_quality} />
          </Box>
          {metrics.z_prime !== null && (
            <Typography
              variant="body2"
              color="text.secondary"
              sx={{ ml: 2, display: { xs: 'none', sm: 'block' } }}
            >
              Z&apos; = {metrics.z_prime.toFixed(2)}
            </Typography>
          )}
        </Box>
      </AccordionSummary>
      <AccordionDetails>
        {/* Plate-level metrics */}
        {hasPlateControls && (
          <>
            <Typography variant="subtitle2" color="text.secondary" sx={{ mb: 1 }}>
              Plate Metrics (from {metrics.n_high_controls} high / {metrics.n_low_controls} low controls)
            </Typography>
            <Grid container spacing={1.5} sx={{ mb: 2 }}>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="Z-prime"
                  value={metrics.z_prime}
                  quality={getZPrimeQuality(metrics.z_prime)}
                  description=">0.5 excellent, 0-0.5 marginal"
                />
              </Grid>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="Robust Z'"
                  value={metrics.robust_z_prime}
                  quality={getZPrimeQuality(metrics.robust_z_prime)}
                  description="Outlier-resistant (median/MAD)"
                />
              </Grid>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="Signal/Background"
                  value={metrics.signal_to_background}
                  description="High ÷ Low control mean"
                />
              </Grid>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="Signal/Noise"
                  value={metrics.signal_to_noise}
                  description="Separation ÷ Low stdev"
                />
              </Grid>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="CV% High"
                  value={metrics.cv_high_controls}
                  quality={getCVQuality(metrics.cv_high_controls)}
                  suffix="%"
                  description="High control variation (<10% good)"
                />
              </Grid>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="CV% Low"
                  value={metrics.cv_low_controls}
                  quality={getCVQuality(metrics.cv_low_controls)}
                  suffix="%"
                  description="Low control variation (<10% good)"
                />
              </Grid>
            </Grid>
          </>
        )}

        {/* Curve fitting metrics */}
        {metrics.total_curves > 0 && (
          <>
            <Typography variant="subtitle2" color="text.secondary" sx={{ mb: 1 }}>
              <TrendingUp sx={{ fontSize: 16, verticalAlign: 'text-bottom', mr: 0.5 }} />
              Curve Fitting ({metrics.total_curves} curves)
            </Typography>
            <Grid container spacing={1.5}>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="Valid Curves"
                  value={metrics.percent_valid}
                  quality={getPercentQuality(metrics.percent_valid)}
                  suffix="%"
                  description={`${metrics.valid_curves} of ${metrics.total_curves} curves`}
                />
              </Grid>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="Good Fit (R²)"
                  value={metrics.percent_good_fit}
                  quality={getPercentQuality(metrics.percent_good_fit, 80)}
                  suffix="%"
                  description={`${metrics.curves_with_good_fit} curves with R² > 0.8`}
                />
              </Grid>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="Invalid"
                  value={metrics.invalid_curves}
                  quality={metrics.invalid_curves > 0 ? 'warning' : 'good'}
                  description="Curves marked invalid"
                />
              </Grid>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="Flagged"
                  value={metrics.curves_with_flags}
                  quality={metrics.curves_with_flags > 0 ? 'warning' : 'good'}
                  description="Curves with QC flags"
                />
              </Grid>
              <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                <MetricCard
                  label="Unassigned"
                  value={metrics.unassigned_curves}
                  description="Awaiting review"
                />
              </Grid>
            </Grid>
          </>
        )}

        {/* Issues */}
        {metrics.issues && metrics.issues.length > 0 && (
          <Box sx={{ mt: 2 }}>
            {metrics.issues.map((issue, idx) => (
              <Alert key={idx} severity="warning" sx={{ mb: idx < metrics.issues.length - 1 ? 1 : 0 }}>
                {issue}
              </Alert>
            ))}
          </Box>
        )}
      </AccordionDetails>
    </Accordion>
  );
}
