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
 * Hill-Langmuir equation for dose-response curve fitting
 * y = minVal + (maxVal - minVal) / (1 + (EC50/x)^Hill)
 */
function hillLangmuir(
  x: number,
  ec50: number,
  hill: number,
  minVal: number,
  maxVal: number
): number {
  if (x <= 0 || ec50 <= 0) return minVal;
  return minVal + (maxVal - minVal) / (1 + Math.pow(ec50 / x, hill));
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

    // Add fitted curve if we have valid parameters
    if (
      fit?.ec50 &&
      fit.ec50 > 0 &&
      fit.hill !== undefined &&
      fit.hill !== null &&
      fit.minVal !== undefined &&
      fit.minVal !== null &&
      fit.maxVal !== undefined &&
      fit.maxVal !== null &&
      fit.status === 'valid'
    ) {
      const curvePoints = generateFittedCurve(
        concentrations,
        fit.ec50,
        fit.hill,
        fit.minVal,
        fit.maxVal
      );

      datasets.push({
        label: `Fitted Curve (EC50: ${fit.ec50.toExponential(2)})`,
        data: curvePoints,
        backgroundColor: 'transparent',
        borderColor: 'rgba(220, 53, 69, 1)',
        borderWidth: 2,
        pointRadius: 0,
        showLine: true,
        tension: 0.4,
      });

      // Add EC50 marker line
      const ec50Y = hillLangmuir(fit.ec50, fit.ec50, fit.hill, fit.minVal, fit.maxVal);
      datasets.push({
        label: 'EC50',
        data: [
          { x: fit.ec50, y: fit.minVal },
          { x: fit.ec50, y: ec50Y },
        ],
        backgroundColor: 'transparent',
        borderColor: 'rgba(220, 53, 69, 0.5)',
        borderWidth: 1,
        borderDash: [5, 5],
        pointRadius: 0,
        showLine: true,
      });
    }

    return { datasets };
  }, [data, fit, skipPoints]);

  const options = useMemo<ChartOptions<'scatter'>>(() => ({
    responsive: true,
    maintainAspectRatio: false,
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
        data: points,
        backgroundColor: 'rgba(70, 130, 180, 0.8)',
        borderColor: 'rgba(70, 130, 180, 1)',
        pointRadius: 3,
      },
    ];

    // Add fitted curve
    if (
      fit?.ec50 &&
      fit.ec50 > 0 &&
      fit.hill !== undefined &&
      fit.minVal !== undefined &&
      fit.maxVal !== undefined &&
      fit.status === 'valid'
    ) {
      const curvePoints = generateFittedCurve(
        concentrations,
        fit.ec50,
        fit.hill!,
        fit.minVal!,
        fit.maxVal!,
        50
      );

      datasets.push({
        data: curvePoints,
        backgroundColor: 'transparent',
        borderColor: 'rgba(220, 53, 69, 1)',
        borderWidth: 1.5,
        pointRadius: 0,
        showLine: true,
        tension: 0.4,
      });
    }

    return { datasets };
  }, [data, fit]);

  const options = useMemo<ChartOptions<'scatter'>>(() => ({
    responsive: true,
    maintainAspectRatio: false,
    plugins: {
      legend: { display: false },
      title: { display: false },
      tooltip: { enabled: false },
    },
    scales: {
      x: {
        type: 'logarithmic' as const,
        display: false,
      },
      y: {
        type: 'linear' as const,
        display: false,
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
