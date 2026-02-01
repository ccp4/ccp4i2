'use client';

import { useMemo } from 'react';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  LogarithmicScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend,
  ChartOptions,
  ChartData,
} from 'chart.js';
import { Scatter } from 'react-chartjs-2';
import { Box, Typography, Paper } from '@mui/material';

// Register Chart.js components
ChartJS.register(
  CategoryScale,
  LinearScale,
  LogarithmicScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend
);

/**
 * Four-parameter logistic (4PL) equation for dose-response curve fitting
 * y = minVal + (maxVal - minVal) / (1 + (x/EC50)^Hill)
 *
 * This matches the Python fitting implementation where:
 * - minVal (bottom) is approached at high concentrations
 * - maxVal (top) is approached at low concentrations
 */
function hillLangmuir(
  x: number,
  ec50: number,
  hill: number,
  minVal: number,
  maxVal: number
): number {
  if (x <= 0 || ec50 <= 0) return maxVal;
  return minVal + (maxVal - minVal) / (1 + Math.pow(x / ec50, hill));
}

/**
 * Generate fitted curve points for smooth line
 */
function generateFittedCurve(
  concentrations: number[],
  ec50: number,
  hill: number,
  minVal: number,
  maxVal: number,
  numPoints: number = 100
): { x: number; y: number }[] {
  if (!concentrations.length) return [];

  const minConc = Math.min(...concentrations.filter(c => c > 0));
  const maxConc = Math.max(...concentrations);

  // Generate logarithmically spaced points
  const logMin = Math.log10(minConc / 10);
  const logMax = Math.log10(maxConc * 10);
  const step = (logMax - logMin) / numPoints;

  const points: { x: number; y: number }[] = [];
  for (let i = 0; i <= numPoints; i++) {
    const x = Math.pow(10, logMin + i * step);
    const y = hillLangmuir(x, ec50, hill, minVal, maxVal);
    points.push({ x, y });
  }

  return points;
}

export interface DoseResponseData {
  concentrations: number[];
  responses: number[];
  unit?: string;
}

export interface FitParameters {
  ec50?: number | null;
  hill?: number | null;
  minVal?: number | null;
  maxVal?: number | null;
  status?: string;
}

interface DoseResponseChartProps {
  data: DoseResponseData;
  fit?: FitParameters;
  title?: string;
  compoundName?: string;
  width?: number;
  height?: number;
  showLegend?: boolean;
  skipPoints?: number[];
}

export function DoseResponseChart({
  data,
  fit,
  title,
  compoundName,
  width = 400,
  height = 300,
  showLegend = true,
  skipPoints = [],
}: DoseResponseChartProps) {
  const chartData = useMemo<ChartData<'scatter'>>(() => {
    const { concentrations, responses } = data;

    // Separate included and skipped points
    const includedPoints: { x: number; y: number }[] = [];
    const skippedPoints: { x: number; y: number }[] = [];

    concentrations.forEach((conc, idx) => {
      if (conc > 0 && responses[idx] !== undefined) {
        const point = { x: conc, y: responses[idx] };
        if (Array.isArray(skipPoints) && skipPoints.includes(idx)) {
          skippedPoints.push(point);
        } else {
          includedPoints.push(point);
        }
      }
    });

    const datasets: ChartData<'scatter'>['datasets'] = [
      {
        label: 'Data Points',
        data: includedPoints,
        backgroundColor: 'rgba(70, 130, 180, 0.8)',
        borderColor: 'rgba(70, 130, 180, 1)',
        pointRadius: 6,
        pointHoverRadius: 8,
      },
    ];

    // Add skipped points if any
    if (skippedPoints.length > 0) {
      datasets.push({
        label: 'Excluded Points',
        data: skippedPoints,
        backgroundColor: 'rgba(200, 200, 200, 0.6)',
        borderColor: 'rgba(150, 150, 150, 1)',
        pointRadius: 5,
        pointHoverRadius: 7,
        pointStyle: 'crossRot',
      });
    }

    // Add fitted curve if we have valid numeric parameters
    const hasValidFitParams =
      fit?.ec50 != null &&
      fit.ec50 > 0 &&
      fit.hill != null &&
      typeof fit.hill === 'number' &&
      fit.minVal != null &&
      typeof fit.minVal === 'number' &&
      fit.maxVal != null &&
      typeof fit.maxVal === 'number';

    if (hasValidFitParams) {
      // Normalize min/max: for standard dose-response, minVal should be < maxVal
      // Legacy data may have these inverted - detect and correct
      let normalizedMin = fit!.minVal!;
      let normalizedMax = fit!.maxVal!;
      if (normalizedMin > normalizedMax) {
        // Swap if inverted (legacy data with opposite interpretation)
        [normalizedMin, normalizedMax] = [normalizedMax, normalizedMin];
      }

      const curvePoints = generateFittedCurve(
        concentrations,
        fit!.ec50!,
        fit!.hill!,
        normalizedMin,
        normalizedMax
      );

      // Use different color for invalid fits
      const curveColor = fit?.status === 'valid'
        ? 'rgba(220, 53, 69, 1)'
        : 'rgba(150, 150, 150, 0.7)';

      // Use type: 'line' for the curve dataset in mixed chart
      datasets.push({
        type: 'line' as const,
        label: `Fitted Curve (EC50: ${fit!.ec50!.toExponential(2)})`,
        data: curvePoints,
        backgroundColor: 'transparent',
        borderColor: curveColor,
        borderWidth: 2,
        pointRadius: 0,
        fill: false,
      } as any);

      // Add EC50 marker line
      const ec50Y = hillLangmuir(fit!.ec50!, fit!.ec50!, fit!.hill!, normalizedMin, normalizedMax);
      datasets.push({
        type: 'line' as const,
        label: 'EC50',
        data: [
          { x: fit!.ec50!, y: normalizedMin },
          { x: fit!.ec50!, y: ec50Y },
        ],
        backgroundColor: 'transparent',
        borderColor: curveColor.replace('1)', '0.5)'),
        borderWidth: 1,
        borderDash: [5, 5],
        pointRadius: 0,
        fill: false,
      } as any);
    }

    return { datasets };
  }, [data, fit, skipPoints]);

  const options = useMemo<ChartOptions<'scatter'>>(() => ({
    responsive: true,
    maintainAspectRatio: false,
    animation: false,
    plugins: {
      legend: {
        display: showLegend,
        position: 'bottom' as const,
        labels: {
          usePointStyle: true,
          filter: (item) => !item.text?.includes('EC50'),
        },
      },
      title: {
        display: !!title || !!compoundName,
        text: title || compoundName || '',
        font: { size: 14 },
      },
      tooltip: {
        callbacks: {
          label: (context) => {
            const point = context.raw as { x: number; y: number };
            return `[${point.x.toExponential(2)}]: ${point.y.toFixed(1)}`;
          },
        },
      },
    },
    scales: {
      x: {
        type: 'logarithmic' as const,
        title: {
          display: true,
          text: `Concentration (${data.unit || 'nM'})`,
        },
        ticks: {
          callback: (value) => {
            const num = Number(value);
            if (num >= 1000) return `${num / 1000}k`;
            if (num >= 1) return num.toString();
            return num.toExponential(0);
          },
        },
      },
      y: {
        type: 'linear' as const,
        title: {
          display: true,
          text: 'Response',
        },
      },
    },
  }), [title, compoundName, showLegend, data.unit]);

  if (!data.concentrations.length || !data.responses.length) {
    return (
      <Paper sx={{ p: 2, width, height, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
        <Typography color="text.secondary">No data available</Typography>
      </Paper>
    );
  }

  return (
    <Box sx={{ width, height }}>
      <Scatter data={chartData} options={options} />
    </Box>
  );
}

/**
 * Compact chart for table thumbnails
 */
interface DoseResponseThumbProps {
  data: DoseResponseData;
  fit?: FitParameters;
  size?: number;
}

export function DoseResponseThumb({ data, fit, size = 120 }: DoseResponseThumbProps) {
  const chartData = useMemo<ChartData<'scatter'>>(() => {
    const { concentrations, responses } = data;

    const points = concentrations
      .map((conc, idx) => ({ x: conc, y: responses[idx] }))
      .filter(p => p.x > 0 && p.y !== undefined);

    const datasets: ChartData<'scatter'>['datasets'] = [
      {
        label: 'Data',
        data: points,
        backgroundColor: 'rgba(70, 130, 180, 0.8)',
        borderColor: 'rgba(70, 130, 180, 1)',
        pointRadius: 3,
      },
    ];

    // Add fitted curve if we have valid numeric parameters
    const hasValidFitParams =
      fit?.ec50 != null &&
      fit.ec50 > 0 &&
      fit.hill != null &&
      typeof fit.hill === 'number' &&
      fit.minVal != null &&
      typeof fit.minVal === 'number' &&
      fit.maxVal != null &&
      typeof fit.maxVal === 'number';

    if (hasValidFitParams) {
      // Normalize min/max: for standard dose-response, minVal should be < maxVal
      // Legacy data may have these inverted - detect and correct
      let normalizedMin = fit!.minVal!;
      let normalizedMax = fit!.maxVal!;
      if (normalizedMin > normalizedMax) {
        [normalizedMin, normalizedMax] = [normalizedMax, normalizedMin];
      }

      const curvePoints = generateFittedCurve(
        concentrations,
        fit!.ec50!,
        fit!.hill!,
        normalizedMin,
        normalizedMax,
        50
      );

      // Use red for valid fits, dark gray for invalid
      const curveColor = fit?.status === 'valid'
        ? 'rgba(220, 53, 69, 1)'
        : 'rgba(100, 100, 100, 0.9)';

      // For scatter charts, use type: 'line' for the curve dataset
      datasets.push({
        type: 'line' as const,
        label: 'Fit',
        data: curvePoints,
        backgroundColor: 'transparent',
        borderColor: curveColor,
        borderWidth: 2,
        pointRadius: 0,
        fill: false,
      } as any);
    }

    return { datasets };
  }, [data, fit]);

  // Calculate axis bounds from data
  const axisBounds = useMemo(() => {
    const { concentrations, responses } = data;
    const validConcs = concentrations.filter(c => c > 0);
    const validResponses = responses.filter(r => r !== undefined && r !== null);

    return {
      xMin: validConcs.length ? Math.min(...validConcs) : 1,
      xMax: validConcs.length ? Math.max(...validConcs) : 1000,
      yMin: validResponses.length ? Math.min(...validResponses) : 0,
      yMax: validResponses.length ? Math.max(...validResponses) : 100,
    };
  }, [data]);

  const options = useMemo<ChartOptions<'scatter'>>(() => ({
    responsive: true,
    maintainAspectRatio: false,
    animation: false,
    layout: {
      padding: { left: 2, right: 2, top: 2, bottom: 2 },
    },
    plugins: {
      legend: { display: false },
      title: { display: false },
      tooltip: { enabled: false },
    },
    scales: {
      x: {
        type: 'logarithmic' as const,
        display: true,
        grid: { display: false },
        border: { display: true, color: 'rgba(0,0,0,0.1)' },
        ticks: {
          display: true,
          font: { size: 8 },
          color: 'rgba(0,0,0,0.5)',
          maxTicksLimit: 2,
          callback: (value) => {
            const num = Number(value);
            if (num >= 1000) return `${(num / 1000).toFixed(0)}k`;
            if (num >= 1) return num.toFixed(0);
            return num.toExponential(0);
          },
        },
      },
      y: {
        type: 'linear' as const,
        display: true,
        grid: { display: false },
        border: { display: true, color: 'rgba(0,0,0,0.1)' },
        ticks: {
          display: true,
          font: { size: 8 },
          color: 'rgba(0,0,0,0.5)',
          maxTicksLimit: 2,
          callback: (value) => {
            const num = Number(value);
            return num.toFixed(0);
          },
        },
      },
    },
    elements: {
      point: { radius: 2 },
    },
  }), []);

  if (!data.concentrations.length) {
    return (
      <Box
        sx={{
          width: size,
          height: size,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          bgcolor: 'grey.100',
          borderRadius: 1,
        }}
      >
        <Typography variant="caption" color="text.secondary">
          No data
        </Typography>
      </Box>
    );
  }

  return (
    <Box sx={{ width: size, height: size }}>
      <Scatter data={chartData} options={options} />
    </Box>
  );
}
