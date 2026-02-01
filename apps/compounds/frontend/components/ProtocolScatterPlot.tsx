'use client';

import { useState, useMemo, useRef, useEffect } from 'react';
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
} from '@mui/material';
import {
  BubbleChart,
  Close,
  Refresh,
} from '@mui/icons-material';
import {
  Chart as ChartJS,
  LinearScale,
  LogarithmicScale,
  PointElement,
  Tooltip as ChartTooltip,
  Legend,
  ChartOptions,
  ChartData,
} from 'chart.js';
import { Scatter } from 'react-chartjs-2';
import { useRouter } from 'next/navigation';
import { CompactRow, ProtocolInfo } from '@/types/aggregation';

// Register Chart.js components
ChartJS.register(
  LinearScale,
  LogarithmicScale,
  PointElement,
  ChartTooltip,
  Legend
);

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
    };
  }, [data, xAxis.protocolId, yAxis.protocolId]);

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
          callbacks: {
            label: (context) => {
              const point = context.raw as {
                x: number;
                y: number;
                formatted_id: string;
                target_name?: string;
              };
              return [
                point.formatted_id,
                point.target_name ? `Target: ${point.target_name}` : '',
                `${xProtocolName}: ${point.x.toFixed(2)}`,
                `${yProtocolName}: ${point.y.toFixed(2)}`,
              ].filter(Boolean);
            },
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
          const dataIndex = elements[0].index;
          const point = (chartData.datasets[0].data as any[])[dataIndex];
          if (point?.compound_id) {
            router.push(`/registry/compounds/${point.compound_id}`);
          }
        }
      },
    }),
    [xAxis, yAxis, xProtocolName, yProtocolName, axisBounds, chartData, router]
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
          </Box>

          {/* Chart */}
          <Box sx={{ flex: 1, minHeight: 0 }}>
            {xAxis.protocolId && yAxis.protocolId ? (
              pointCount > 0 ? (
                <Box sx={{ height: '100%' }}>
                  <Scatter ref={chartRef} data={chartData} options={chartOptions} />
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
        </DialogContent>
      </Dialog>
    </>
  );
}
