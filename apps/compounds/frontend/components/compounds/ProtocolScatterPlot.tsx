'use client';

import { useState, useMemo, useRef, useEffect, useCallback } from 'react';
import {
  Box,
  Button,
  Dialog,
  DialogTitle,
  DialogContent,
  Typography,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
  IconButton,
  Tooltip,
  Chip,
  Switch,
  FormControlLabel,
  Paper,
  Popper,
  Fade,
  Skeleton,
  Divider,
} from '@mui/material';
import {
  BubbleChart,
  Close,
  Refresh,
  TrendingUp,
} from '@mui/icons-material';
import { useRDKit } from '@/lib/compounds/rdkit-context';
import {
  Chart as ChartJS,
  LinearScale,
  LogarithmicScale,
  PointElement,
  LineElement,
  Tooltip as ChartTooltip,
  Legend,
  ChartOptions,
  ChartData,
} from 'chart.js';
import { Scatter } from 'react-chartjs-2';
import { useRouter } from 'next/navigation';
import { CompactRow, ProtocolInfo } from '@/types/compounds/aggregation';

// Register Chart.js components
ChartJS.register(
  LinearScale,
  LogarithmicScale,
  PointElement,
  LineElement,
  ChartTooltip,
  Legend
);

/**
 * Linear regression result interface
 */
interface RegressionResult {
  slope: number;
  intercept: number;
  rSquared: number;
  // For display: in log-log space, y = a * x^b where a = 10^intercept, b = slope
  // For display: in linear space, y = slope * x + intercept
}

/**
 * Calculate linear regression, optionally in log-space for each axis
 * When both axes are logarithmic, we fit log(y) = slope * log(x) + intercept
 * This corresponds to y = 10^intercept * x^slope (power law)
 */
function calculateRegression(
  points: Array<{ x: number; y: number }>,
  xLog: boolean,
  yLog: boolean
): RegressionResult | null {
  if (points.length < 2) return null;

  // Transform values based on scale type
  const transformedPoints = points
    .filter((p) => p.x > 0 && p.y > 0) // Must be positive for log
    .map((p) => ({
      x: xLog ? Math.log10(p.x) : p.x,
      y: yLog ? Math.log10(p.y) : p.y,
    }));

  if (transformedPoints.length < 2) return null;

  const n = transformedPoints.length;
  const sumX = transformedPoints.reduce((acc, p) => acc + p.x, 0);
  const sumY = transformedPoints.reduce((acc, p) => acc + p.y, 0);
  const sumXY = transformedPoints.reduce((acc, p) => acc + p.x * p.y, 0);
  const sumXX = transformedPoints.reduce((acc, p) => acc + p.x * p.x, 0);
  const sumYY = transformedPoints.reduce((acc, p) => acc + p.y * p.y, 0);

  const denominator = n * sumXX - sumX * sumX;
  if (Math.abs(denominator) < 1e-10) return null;

  const slope = (n * sumXY - sumX * sumY) / denominator;
  const intercept = (sumY - slope * sumX) / n;

  // Calculate R² (coefficient of determination)
  const meanY = sumY / n;
  const ssTotal = sumYY - n * meanY * meanY;
  const ssResidual = transformedPoints.reduce((acc, p) => {
    const predicted = slope * p.x + intercept;
    return acc + (p.y - predicted) ** 2;
  }, 0);

  const rSquared = ssTotal > 0 ? 1 - ssResidual / ssTotal : 0;

  return { slope, intercept, rSquared };
}

interface ProtocolScatterPlotProps {
  /** Compact aggregation data rows */
  data: CompactRow[];
  /** Available protocols with geomean values */
  protocols: ProtocolInfo[];
}

interface AxisConfig {
  protocolId: string | null;
  scale: 'linear' | 'logarithmic';
  min: string;
  max: string;
}

/**
 * Scatter plot component for comparing geomean values between two protocols.
 * Each point represents a compound, with coordinates being the geomean values
 * for the selected X and Y protocols.
 */
export function ProtocolScatterPlot({ data, protocols }: ProtocolScatterPlotProps) {
  const router = useRouter();
  const chartRef = useRef<ChartJS<'scatter'>>(null);

  const [open, setOpen] = useState(false);
  const [xAxis, setXAxis] = useState<AxisConfig>({
    protocolId: protocols.length > 0 ? protocols[0].id : null,
    scale: 'logarithmic',
    min: '',
    max: '',
  });
  const [yAxis, setYAxis] = useState<AxisConfig>({
    protocolId: protocols.length > 1 ? protocols[1].id : null,
    scale: 'logarithmic',
    min: '',
    max: '',
  });
  const [showRegression, setShowRegression] = useState(false);

  // Custom tooltip state for showing molecule structure
  const [hoveredPoint, setHoveredPoint] = useState<{
    x: number;
    y: number;
    compound_id: string;
    formatted_id: string;
    smiles?: string;
    target_name?: string;
    screenX: number;
    screenY: number;
  } | null>(null);
  const tooltipAnchorRef = useRef<HTMLDivElement>(null);
  const { rdkitModule, isLoading: rdkitLoading } = useRDKit();

  // Generate molecule SVG for hovered point
  const moleculeSvgUrl = useMemo(() => {
    if (!hoveredPoint?.smiles || !rdkitModule) return null;
    try {
      const mol = rdkitModule.get_mol(hoveredPoint.smiles);
      if (!mol) return null;
      const svg = mol.get_svg(180, 180);
      mol.delete();
      const blob = new Blob([svg], { type: 'image/svg+xml' });
      return URL.createObjectURL(blob);
    } catch {
      return null;
    }
  }, [hoveredPoint?.smiles, rdkitModule]);

  // Clean up SVG URL when it changes
  useEffect(() => {
    return () => {
      if (moleculeSvgUrl) {
        URL.revokeObjectURL(moleculeSvgUrl);
      }
    };
  }, [moleculeSvgUrl]);

  // Clear hovered point when dialog closes or protocols change
  useEffect(() => {
    if (!open) {
      setHoveredPoint(null);
    }
  }, [open]);

  useEffect(() => {
    setHoveredPoint(null);
  }, [xAxis.protocolId, yAxis.protocolId]);

  // Update axis selections when protocols change
  useEffect(() => {
    if (protocols.length > 0 && !xAxis.protocolId) {
      setXAxis((prev) => ({ ...prev, protocolId: protocols[0].id }));
    }
    if (protocols.length > 1 && !yAxis.protocolId) {
      setYAxis((prev) => ({ ...prev, protocolId: protocols[1].id }));
    }
  }, [protocols, xAxis.protocolId, yAxis.protocolId]);

  // Get protocol names for labels
  const xProtocolName = useMemo(
    () => protocols.find((p) => p.id === xAxis.protocolId)?.name || 'Select Protocol',
    [protocols, xAxis.protocolId]
  );
  const yProtocolName = useMemo(
    () => protocols.find((p) => p.id === yAxis.protocolId)?.name || 'Select Protocol',
    [protocols, yAxis.protocolId]
  );

  // Build scatter plot data
  const chartData: ChartData<'scatter'> = useMemo(() => {
    if (!xAxis.protocolId || !yAxis.protocolId) {
      return { datasets: [] };
    }

    const points: Array<{
      x: number;
      y: number;
      compound_id: string;
      formatted_id: string;
      smiles?: string;
      target_name?: string;
    }> = [];

    for (const row of data) {
      const xValue = row.protocols[xAxis.protocolId]?.geomean;
      const yValue = row.protocols[yAxis.protocolId]?.geomean;

      // Only include points that have both values
      if (xValue != null && yValue != null && xValue > 0 && yValue > 0) {
        points.push({
          x: xValue,
          y: yValue,
          compound_id: row.compound_id,
          formatted_id: row.formatted_id,
          smiles: row.smiles ?? undefined,
          target_name: row.target_name ?? undefined,
        });
      }
    }

    return {
      datasets: [
        {
          label: 'Compounds',
          data: points,
          backgroundColor: 'rgba(25, 118, 210, 0.6)',
          borderColor: 'rgba(25, 118, 210, 1)',
          borderWidth: 1,
          pointRadius: 6,
          pointHoverRadius: 8,
          pointHoverBackgroundColor: 'rgba(25, 118, 210, 0.8)',
        },
      ],
      // Store raw points for regression calculation
      _rawPoints: points,
    };
  }, [data, xAxis.protocolId, yAxis.protocolId]);

  // Calculate regression based on current scale settings
  const regression = useMemo(() => {
    const points = (chartData as any)._rawPoints as Array<{ x: number; y: number }> || [];
    if (points.length < 2) return null;
    return calculateRegression(
      points,
      xAxis.scale === 'logarithmic',
      yAxis.scale === 'logarithmic'
    );
  }, [chartData, xAxis.scale, yAxis.scale]);

  // Calculate axis bounds
  const axisBounds = useMemo(() => {
    const points = chartData.datasets[0]?.data as Array<{ x: number; y: number }> || [];
    if (points.length === 0) {
      return { xMin: 0.01, xMax: 1000, yMin: 0.01, yMax: 1000 };
    }

    const xValues = points.map((p) => p.x);
    const yValues = points.map((p) => p.y);

    const xMin = Math.min(...xValues);
    const xMax = Math.max(...xValues);
    const yMin = Math.min(...yValues);
    const yMax = Math.max(...yValues);

    // Add padding (10% on each side for log scale)
    const xPadding = xAxis.scale === 'logarithmic' ? Math.pow(10, 0.1) : (xMax - xMin) * 0.1;
    const yPadding = yAxis.scale === 'logarithmic' ? Math.pow(10, 0.1) : (yMax - yMin) * 0.1;

    return {
      xMin: xAxis.scale === 'logarithmic' ? xMin / xPadding : xMin - xPadding,
      xMax: xAxis.scale === 'logarithmic' ? xMax * xPadding : xMax + xPadding,
      yMin: yAxis.scale === 'logarithmic' ? yMin / yPadding : yMin - yPadding,
      yMax: yAxis.scale === 'logarithmic' ? yMax * yPadding : yMax + yPadding,
    };
  }, [chartData, xAxis.scale, yAxis.scale]);

  // Generate regression line points for display
  const regressionLineData = useMemo(() => {
    if (!showRegression || !regression) return null;

    const xLog = xAxis.scale === 'logarithmic';
    const yLog = yAxis.scale === 'logarithmic';

    // Get effective axis bounds
    const xMin = xAxis.min ? parseFloat(xAxis.min) : axisBounds.xMin;
    const xMax = xAxis.max ? parseFloat(xAxis.max) : axisBounds.xMax;

    // Generate points for the regression line
    // Use more points for log scale to ensure smooth curve appearance
    const numPoints = 100;
    const linePoints: Array<{ x: number; y: number }> = [];

    for (let i = 0; i <= numPoints; i++) {
      let x: number;
      if (xLog) {
        // Logarithmic spacing
        const logMin = Math.log10(Math.max(xMin, 1e-10));
        const logMax = Math.log10(Math.max(xMax, 1e-10));
        x = Math.pow(10, logMin + (i / numPoints) * (logMax - logMin));
      } else {
        // Linear spacing
        x = xMin + (i / numPoints) * (xMax - xMin);
      }

      // Calculate y using the regression equation in the appropriate space
      let y: number;
      if (xLog && yLog) {
        // log(y) = slope * log(x) + intercept
        y = Math.pow(10, regression.slope * Math.log10(x) + regression.intercept);
      } else if (xLog && !yLog) {
        // y = slope * log(x) + intercept
        y = regression.slope * Math.log10(x) + regression.intercept;
      } else if (!xLog && yLog) {
        // log(y) = slope * x + intercept
        y = Math.pow(10, regression.slope * x + regression.intercept);
      } else {
        // y = slope * x + intercept
        y = regression.slope * x + regression.intercept;
      }

      // Only include points within reasonable bounds
      if (y > 0 && isFinite(y)) {
        linePoints.push({ x, y });
      }
    }

    return linePoints;
  }, [showRegression, regression, xAxis, yAxis, axisBounds]);

  // Final chart data including regression line
  const finalChartData: ChartData<'scatter'> = useMemo(() => {
    const datasets = [...chartData.datasets];

    if (showRegression && regressionLineData && regressionLineData.length > 0) {
      datasets.push({
        label: 'Regression Line',
        data: regressionLineData,
        backgroundColor: 'transparent',
        borderColor: 'rgba(211, 47, 47, 0.8)',
        borderWidth: 2,
        pointRadius: 0,
        pointHoverRadius: 0,
        showLine: true,
        tension: 0,
        order: 1, // Draw behind scatter points
      } as any);
    }

    return { datasets };
  }, [chartData, showRegression, regressionLineData]);

  // Chart options
  const chartOptions: ChartOptions<'scatter'> = useMemo(
    () => ({
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          display: false,
        },
        tooltip: {
          enabled: false, // Disable built-in tooltip, we use custom
          external: (context) => {
            const { tooltip, chart } = context;

            // Hide tooltip when not visible or hovering regression line
            if (tooltip.opacity === 0 || !tooltip.dataPoints || tooltip.dataPoints.length === 0) {
              setHoveredPoint(null);
              return;
            }

            const dataPoint = tooltip.dataPoints[0];
            // Skip regression line points
            if (dataPoint.dataset.label === 'Regression Line') {
              setHoveredPoint(null);
              return;
            }

            const point = dataPoint.raw as {
              x: number;
              y: number;
              compound_id: string;
              formatted_id: string;
              smiles?: string;
              target_name?: string;
            };

            // Skip if no formatted_id (regression line points)
            if (!point.formatted_id) {
              setHoveredPoint(null);
              return;
            }

            // Get screen position from canvas position
            const canvas = chart.canvas;
            const rect = canvas.getBoundingClientRect();
            const screenX = rect.left + tooltip.caretX;
            const screenY = rect.top + tooltip.caretY;

            setHoveredPoint({
              ...point,
              screenX,
              screenY,
            });
          },
        },
      },
      scales: {
        x: {
          type: xAxis.scale,
          title: {
            display: true,
            text: xProtocolName,
            font: { weight: 'bold' },
          },
          min: xAxis.min ? parseFloat(xAxis.min) : axisBounds.xMin,
          max: xAxis.max ? parseFloat(xAxis.max) : axisBounds.xMax,
          ticks: {
            callback: (value) => {
              if (typeof value === 'number') {
                return value < 0.01 || value > 10000 ? value.toExponential(0) : value;
              }
              return value;
            },
          },
        },
        y: {
          type: yAxis.scale,
          title: {
            display: true,
            text: yProtocolName,
            font: { weight: 'bold' },
          },
          min: yAxis.min ? parseFloat(yAxis.min) : axisBounds.yMin,
          max: yAxis.max ? parseFloat(yAxis.max) : axisBounds.yMax,
          ticks: {
            callback: (value) => {
              if (typeof value === 'number') {
                return value < 0.01 || value > 10000 ? value.toExponential(0) : value;
              }
              return value;
            },
          },
        },
      },
      onClick: (_, elements) => {
        if (elements.length > 0) {
          const element = elements[0];
          // Only navigate for scatter points (dataset index 0), not regression line
          if (element.datasetIndex === 0) {
            const point = (finalChartData.datasets[0].data as any[])[element.index];
            if (point?.compound_id) {
              router.push(`/registry/compounds/${point.compound_id}`);
            }
          }
        }
      },
    }),
    [xAxis, yAxis, xProtocolName, yProtocolName, axisBounds, finalChartData, router, setHoveredPoint]
  );

  // Reset axis ranges to auto
  const handleResetRanges = () => {
    setXAxis((prev) => ({ ...prev, min: '', max: '' }));
    setYAxis((prev) => ({ ...prev, min: '', max: '' }));
  };

  // Count points that will be shown
  const pointCount = chartData.datasets[0]?.data?.length || 0;

  // Only show if we have at least 2 protocols with geomean data
  if (protocols.length < 2) {
    return null;
  }

  return (
    <>
      <Button
        variant="outlined"
        startIcon={<BubbleChart />}
        onClick={() => setOpen(true)}
        sx={{ mb: 2 }}
      >
        Compare Protocols
        {pointCount > 0 && (
          <Chip
            label={`${pointCount} compounds`}
            size="small"
            variant="outlined"
            sx={{ ml: 1 }}
          />
        )}
      </Button>

      <Dialog
        open={open}
        onClose={() => setOpen(false)}
        maxWidth="lg"
        fullWidth
        PaperProps={{ sx: { height: '80vh' } }}
      >
        <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
            <BubbleChart color="primary" />
            <Typography variant="h6">Protocol Comparison Plot</Typography>
          </Box>
          <IconButton onClick={() => setOpen(false)} size="small">
            <Close />
          </IconButton>
        </DialogTitle>

        <DialogContent sx={{ display: 'flex', flexDirection: 'column' }}>
          {/* Axis controls */}
          <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap', mb: 2 }}>
            {/* X-Axis controls */}
            <Box sx={{ flex: '1 1 300px', minWidth: 280 }}>
              <Typography variant="caption" color="text.secondary" gutterBottom>
                X-Axis
              </Typography>
              <Box sx={{ display: 'flex', gap: 1, alignItems: 'flex-start' }}>
                <FormControl size="small" sx={{ minWidth: 150 }}>
                  <InputLabel>Protocol</InputLabel>
                  <Select
                    value={xAxis.protocolId || ''}
                    label="Protocol"
                    onChange={(e) =>
                      setXAxis((prev) => ({ ...prev, protocolId: e.target.value }))
                    }
                  >
                    {protocols.map((p) => (
                      <MenuItem key={p.id} value={p.id}>
                        {p.name}
                      </MenuItem>
                    ))}
                  </Select>
                </FormControl>
                <ToggleButtonGroup
                  value={xAxis.scale}
                  exclusive
                  onChange={(_, value) =>
                    value && setXAxis((prev) => ({ ...prev, scale: value }))
                  }
                  size="small"
                >
                  <ToggleButton value="logarithmic">Log</ToggleButton>
                  <ToggleButton value="linear">Lin</ToggleButton>
                </ToggleButtonGroup>
              </Box>
              <Box sx={{ display: 'flex', gap: 1, mt: 1 }}>
                <TextField
                  label="Min"
                  size="small"
                  type="number"
                  value={xAxis.min}
                  onChange={(e) =>
                    setXAxis((prev) => ({ ...prev, min: e.target.value }))
                  }
                  sx={{ width: 80 }}
                  inputProps={{ step: 'any' }}
                />
                <TextField
                  label="Max"
                  size="small"
                  type="number"
                  value={xAxis.max}
                  onChange={(e) =>
                    setXAxis((prev) => ({ ...prev, max: e.target.value }))
                  }
                  sx={{ width: 80 }}
                  inputProps={{ step: 'any' }}
                />
              </Box>
            </Box>

            {/* Y-Axis controls */}
            <Box sx={{ flex: '1 1 300px', minWidth: 280 }}>
              <Typography variant="caption" color="text.secondary" gutterBottom>
                Y-Axis
              </Typography>
              <Box sx={{ display: 'flex', gap: 1, alignItems: 'flex-start' }}>
                <FormControl size="small" sx={{ minWidth: 150 }}>
                  <InputLabel>Protocol</InputLabel>
                  <Select
                    value={yAxis.protocolId || ''}
                    label="Protocol"
                    onChange={(e) =>
                      setYAxis((prev) => ({ ...prev, protocolId: e.target.value }))
                    }
                  >
                    {protocols.map((p) => (
                      <MenuItem key={p.id} value={p.id}>
                        {p.name}
                      </MenuItem>
                    ))}
                  </Select>
                </FormControl>
                <ToggleButtonGroup
                  value={yAxis.scale}
                  exclusive
                  onChange={(_, value) =>
                    value && setYAxis((prev) => ({ ...prev, scale: value }))
                  }
                  size="small"
                >
                  <ToggleButton value="logarithmic">Log</ToggleButton>
                  <ToggleButton value="linear">Lin</ToggleButton>
                </ToggleButtonGroup>
              </Box>
              <Box sx={{ display: 'flex', gap: 1, mt: 1 }}>
                <TextField
                  label="Min"
                  size="small"
                  type="number"
                  value={yAxis.min}
                  onChange={(e) =>
                    setYAxis((prev) => ({ ...prev, min: e.target.value }))
                  }
                  sx={{ width: 80 }}
                  inputProps={{ step: 'any' }}
                />
                <TextField
                  label="Max"
                  size="small"
                  type="number"
                  value={yAxis.max}
                  onChange={(e) =>
                    setYAxis((prev) => ({ ...prev, max: e.target.value }))
                  }
                  sx={{ width: 80 }}
                  inputProps={{ step: 'any' }}
                />
              </Box>
            </Box>

            {/* Reset button */}
            <Box sx={{ display: 'flex', alignItems: 'flex-end', pb: 0.5 }}>
              <Tooltip title="Reset axis ranges to auto">
                <IconButton onClick={handleResetRanges} size="small">
                  <Refresh />
                </IconButton>
              </Tooltip>
            </Box>

            {/* Regression controls */}
            <Box sx={{ display: 'flex', alignItems: 'flex-end', pb: 0.5, ml: 'auto' }}>
              <FormControlLabel
                control={
                  <Switch
                    checked={showRegression}
                    onChange={(e) => setShowRegression(e.target.checked)}
                    size="small"
                  />
                }
                label={
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                    <TrendingUp fontSize="small" />
                    <Typography variant="body2">Regression</Typography>
                  </Box>
                }
              />
            </Box>
          </Box>

          {/* Regression statistics */}
          {showRegression && regression && (
            <Box
              sx={{
                display: 'flex',
                gap: 2,
                mb: 1,
                p: 1,
                bgcolor: 'grey.50',
                borderRadius: 1,
                flexWrap: 'wrap',
              }}
            >
              <Typography variant="body2">
                <strong>R²:</strong> {regression.rSquared.toFixed(4)}
              </Typography>
              <Typography variant="body2">
                <strong>
                  {xAxis.scale === 'logarithmic' && yAxis.scale === 'logarithmic'
                    ? 'Exponent (b):'
                    : 'Slope:'}
                </strong>{' '}
                {regression.slope.toFixed(4)}
              </Typography>
              <Typography variant="body2">
                <strong>
                  {xAxis.scale === 'logarithmic' && yAxis.scale === 'logarithmic'
                    ? 'Coefficient (a):'
                    : 'Intercept:'}
                </strong>{' '}
                {xAxis.scale === 'logarithmic' && yAxis.scale === 'logarithmic'
                  ? Math.pow(10, regression.intercept).toExponential(3)
                  : regression.intercept.toFixed(4)}
              </Typography>
              <Typography variant="caption" color="text.secondary" sx={{ width: '100%' }}>
                {xAxis.scale === 'logarithmic' && yAxis.scale === 'logarithmic'
                  ? 'Power law fit: y = a × x^b'
                  : xAxis.scale === 'logarithmic'
                    ? 'Semi-log fit: y = slope × log₁₀(x) + intercept'
                    : yAxis.scale === 'logarithmic'
                      ? 'Semi-log fit: log₁₀(y) = slope × x + intercept'
                      : 'Linear fit: y = slope × x + intercept'}
              </Typography>
            </Box>
          )}

          {/* Chart */}
          <Box sx={{ flex: 1, minHeight: 0 }}>
            {xAxis.protocolId && yAxis.protocolId ? (
              pointCount > 0 ? (
                <Box sx={{ height: '100%' }}>
                  <Scatter ref={chartRef} data={finalChartData} options={chartOptions} />
                </Box>
              ) : (
                <Box
                  sx={{
                    height: '100%',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    bgcolor: 'grey.50',
                    borderRadius: 1,
                  }}
                >
                  <Typography color="text.secondary">
                    No compounds have geomean values for both selected protocols
                  </Typography>
                </Box>
              )
            ) : (
              <Box
                sx={{
                  height: '100%',
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center',
                  bgcolor: 'grey.50',
                  borderRadius: 1,
                }}
              >
                <Typography color="text.secondary">
                  Select protocols for both X and Y axes
                </Typography>
              </Box>
            )}
          </Box>

          <Typography variant="caption" color="text.secondary" sx={{ mt: 1 }}>
            Click on a point to view compound details
          </Typography>

          {/* Hidden anchor for tooltip positioning */}
          <Box
            ref={tooltipAnchorRef}
            sx={{
              position: 'fixed',
              left: hoveredPoint?.screenX ?? 0,
              top: hoveredPoint?.screenY ?? 0,
              width: 1,
              height: 1,
              pointerEvents: 'none',
            }}
          />

          {/* Custom tooltip with molecule structure */}
          <Popper
            open={!!hoveredPoint}
            anchorEl={tooltipAnchorRef.current}
            placement="right-start"
            transition
            sx={{ zIndex: 1400, pointerEvents: 'none' }}
            modifiers={[
              {
                name: 'flip',
                enabled: true,
                options: {
                  fallbackPlacements: ['left-start', 'bottom', 'top'],
                },
              },
              {
                name: 'preventOverflow',
                enabled: true,
                options: {
                  boundary: 'viewport',
                  padding: 8,
                },
              },
              {
                name: 'offset',
                options: {
                  offset: [0, 10],
                },
              },
            ]}
          >
            {({ TransitionProps }) => (
              <Fade {...TransitionProps} timeout={150}>
                <Paper
                  elevation={8}
                  sx={{
                    p: 1.5,
                    maxWidth: 280,
                    bgcolor: 'background.paper',
                    border: 1,
                    borderColor: 'divider',
                  }}
                >
                  {hoveredPoint && (
                    <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
                      {/* Compound ID and target */}
                      <Box>
                        <Typography variant="subtitle2" fontWeight="bold">
                          {hoveredPoint.formatted_id}
                        </Typography>
                        {hoveredPoint.target_name && (
                          <Typography variant="caption" color="text.secondary">
                            {hoveredPoint.target_name}
                          </Typography>
                        )}
                      </Box>

                      {/* Molecule structure */}
                      {hoveredPoint.smiles && (
                        <>
                          <Divider />
                          <Box
                            sx={{
                              width: 180,
                              height: 180,
                              display: 'flex',
                              alignItems: 'center',
                              justifyContent: 'center',
                              bgcolor: 'grey.50',
                              borderRadius: 1,
                              mx: 'auto',
                            }}
                          >
                            {rdkitLoading ? (
                              <Skeleton variant="rectangular" width={180} height={180} />
                            ) : moleculeSvgUrl ? (
                              <img
                                src={moleculeSvgUrl}
                                alt={`Structure: ${hoveredPoint.formatted_id}`}
                                style={{ width: 180, height: 180, objectFit: 'contain' }}
                              />
                            ) : (
                              <Typography variant="caption" color="text.secondary">
                                Unable to render
                              </Typography>
                            )}
                          </Box>
                        </>
                      )}

                      {/* Values */}
                      <Divider />
                      <Box sx={{ display: 'flex', flexDirection: 'column', gap: 0.25 }}>
                        <Typography variant="body2">
                          <strong>{xProtocolName}:</strong> {hoveredPoint.x.toFixed(2)}
                        </Typography>
                        <Typography variant="body2">
                          <strong>{yProtocolName}:</strong> {hoveredPoint.y.toFixed(2)}
                        </Typography>
                      </Box>
                    </Box>
                  )}
                </Paper>
              </Fade>
            )}
          </Popper>
        </DialogContent>
      </Dialog>
    </>
  );
}
