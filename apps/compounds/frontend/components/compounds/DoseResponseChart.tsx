/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
'use client';

import { useMemo } from 'react';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  LogarithmicScale,
  PointElement,
  LineElement,
  LineController,
  Title,
  Tooltip,
  Legend,
  ChartOptions,
  ChartData,
} from 'chart.js';
import { Scatter } from 'react-chartjs-2';
import { Box, Typography, Paper, useTheme } from '@mui/material';

// Register Chart.js components
// LineController is needed for mixed chart datasets with type: 'line'
ChartJS.register(
  CategoryScale,
  LinearScale,
  LogarithmicScale,
  PointElement,
  LineElement,
  LineController,
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
  /** Observed low-signal (min) control response value */
  minControlResponse?: number | null;
  /** Observed high-signal (max) control response value */
  maxControlResponse?: number | null;
}

export interface FitParameters {
  ec50?: number | null;
  hill?: number | null;
  minVal?: number | null;
  maxVal?: number | null;
  status?: string;
  /** Pre-generated curve points from backend fitting script [[x1,y1], [x2,y2], ...] */
  curvePoints?: [number, number][] | null;
  /** The KPI name for display (e.g., 'EC50', 'IC50', 'Ki') */
  kpiName?: string | null;
  /** The fitting algorithm used (e.g., 'four-parameter-logistic', 'tight-binding-wang') */
  algorithm?: string | null;
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
  const theme = useTheme();
  const isDarkMode = theme.palette.mode === 'dark';

  // Theme-aware colors
  const colors = useMemo(() => ({
    dataPoint: isDarkMode ? 'rgba(100, 181, 246, 0.9)' : 'rgba(70, 130, 180, 0.8)',
    dataPointBorder: isDarkMode ? 'rgba(100, 181, 246, 1)' : 'rgba(70, 130, 180, 1)',
    skippedPoint: isDarkMode ? 'rgba(120, 120, 120, 0.6)' : 'rgba(200, 200, 200, 0.6)',
    skippedPointBorder: isDarkMode ? 'rgba(100, 100, 100, 1)' : 'rgba(150, 150, 150, 1)',
    validCurve: isDarkMode ? 'rgba(239, 83, 80, 1)' : 'rgba(220, 53, 69, 1)',
    invalidCurve: isDarkMode ? 'rgba(120, 120, 120, 0.7)' : 'rgba(150, 150, 150, 0.7)',
    textColor: theme.palette.text.primary,
    secondaryText: theme.palette.text.secondary,
    gridColor: isDarkMode ? 'rgba(255, 255, 255, 0.1)' : 'rgba(0, 0, 0, 0.1)',
  }), [isDarkMode, theme.palette]);

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
        backgroundColor: colors.dataPoint,
        borderColor: colors.dataPointBorder,
        pointRadius: 6,
        pointHoverRadius: 8,
      },
    ];

    // Add skipped points if any
    if (skippedPoints.length > 0) {
      datasets.push({
        label: 'Excluded Points',
        data: skippedPoints,
        backgroundColor: colors.skippedPoint,
        borderColor: colors.skippedPointBorder,
        pointRadius: 5,
        pointHoverRadius: 7,
        pointStyle: 'crossRot',
      });
    }

    // Determine KPI name for display (use provided name or default to EC50)
    const kpiName = fit?.kpiName || 'EC50';
    const kpiValue = fit?.ec50;

    // Add fitted curve - prefer backend curve_points for algorithm-agnostic rendering
    const hasBackendCurve = fit?.curvePoints && Array.isArray(fit.curvePoints) && fit.curvePoints.length > 0;

    // Fallback: generate curve using 4PL equation if no backend curve and we have valid 4PL params
    const hasValid4PLParams =
      !hasBackendCurve &&
      fit?.ec50 != null &&
      fit.ec50 > 0 &&
      fit.hill != null &&
      typeof fit.hill === 'number' &&
      fit.minVal != null &&
      typeof fit.minVal === 'number' &&
      fit.maxVal != null &&
      typeof fit.maxVal === 'number';

    if (hasBackendCurve || hasValid4PLParams) {
      let curvePoints: { x: number; y: number }[];

      if (hasBackendCurve) {
        // Use backend-generated curve points (works for all algorithms)
        curvePoints = fit!.curvePoints!.map(([x, y]) => ({ x, y }));
      } else {
        // Fallback: generate using 4PL equation
        // Normalize min/max: for standard dose-response, minVal should be < maxVal
        let normalizedMin = fit!.minVal!;
        let normalizedMax = fit!.maxVal!;
        if (normalizedMin > normalizedMax) {
          [normalizedMin, normalizedMax] = [normalizedMax, normalizedMin];
        }

        curvePoints = generateFittedCurve(
          concentrations,
          fit!.ec50!,
          fit!.hill!,
          normalizedMin,
          normalizedMax
        );
      }

      // Use different color for invalid fits
      const curveColor = fit?.status === 'valid'
        ? colors.validCurve
        : colors.invalidCurve;

      // Format KPI value for legend
      const kpiDisplay = kpiValue != null && kpiValue > 0
        ? `${kpiName}: ${kpiValue.toExponential(2)}`
        : kpiName;

      // Use type: 'line' for the curve dataset in mixed chart
      datasets.push({
        type: 'line' as const,
        label: `Fitted Curve (${kpiDisplay})`,
        data: curvePoints,
        backgroundColor: 'transparent',
        borderColor: curveColor,
        borderWidth: 2,
        pointRadius: 0,
        fill: false,
      } as any);

      // Add KPI marker line (only if we have the value and 4PL-style params for the calculation)
      if (kpiValue != null && kpiValue > 0 && fit?.hill != null && fit?.minVal != null && fit?.maxVal != null) {
        let normalizedMin = fit.minVal;
        let normalizedMax = fit.maxVal;
        if (normalizedMin > normalizedMax) {
          [normalizedMin, normalizedMax] = [normalizedMax, normalizedMin];
        }

        const kpiY = hillLangmuir(kpiValue, kpiValue, fit.hill, normalizedMin, normalizedMax);
        datasets.push({
          type: 'line' as const,
          label: kpiName,
          data: [
            { x: kpiValue, y: normalizedMin },
            { x: kpiValue, y: kpiY },
          ],
          backgroundColor: 'transparent',
          borderColor: curveColor.replace('1)', '0.5)'),
          borderWidth: 1,
          borderDash: [5, 5],
          pointRadius: 0,
          fill: false,
        } as any);
      }
    }

    return { datasets };
  }, [data, fit, skipPoints, colors]);

  // Calculate y-axis bounds from data points, fitted asymptotes, AND control values.
  // Including controls ensures that flat/inactive data looks flat in the plot
  // rather than being visually exaggerated by a poor fit's narrow range.
  const yAxisBounds = useMemo(() => {
    const { responses } = data;
    const validResponses = responses.filter(r => r !== undefined && r !== null) as number[];

    if (validResponses.length === 0) {
      return { min: 0, max: 100 };
    }

    let yMin = Math.min(...validResponses);
    let yMax = Math.max(...validResponses);

    // Include fitted asymptotes if available
    if (fit?.minVal != null && typeof fit.minVal === 'number') {
      yMin = Math.min(yMin, fit.minVal);
    }
    if (fit?.maxVal != null && typeof fit.maxVal === 'number') {
      yMax = Math.max(yMax, fit.maxVal);
    }

    // Include observed control values so the axis reflects the full assay range
    if (data.minControlResponse != null) {
      yMin = Math.min(yMin, data.minControlResponse);
      yMax = Math.max(yMax, data.minControlResponse);
    }
    if (data.maxControlResponse != null) {
      yMin = Math.min(yMin, data.maxControlResponse);
      yMax = Math.max(yMax, data.maxControlResponse);
    }

    // Add 5% padding for visual clarity
    const range = yMax - yMin;
    const padding = range * 0.05;

    return {
      min: yMin - padding,
      max: yMax + padding,
    };
  }, [data, fit]);

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
          color: colors.textColor,
          // Hide the KPI marker line from legend (it's redundant with the curve label)
          filter: (item) => {
            const kpiName = fit?.kpiName || 'EC50';
            return item.text !== kpiName;
          },
        },
      },
      title: {
        display: !!title || !!compoundName,
        text: title || compoundName || '',
        font: { size: 14 },
        color: colors.textColor,
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
          color: colors.textColor,
        },
        ticks: {
          color: colors.secondaryText,
          callback: (value) => {
            const num = Number(value);
            if (num >= 1000) return `${num / 1000}k`;
            if (num >= 1) return num.toString();
            return num.toExponential(0);
          },
        },
        grid: {
          color: colors.gridColor,
        },
      },
      y: {
        type: 'linear' as const,
        min: yAxisBounds.min,
        max: yAxisBounds.max,
        title: {
          display: true,
          text: 'Response',
          color: colors.textColor,
        },
        ticks: {
          color: colors.secondaryText,
        },
        grid: {
          color: colors.gridColor,
        },
      },
    },
  }), [title, compoundName, showLegend, data.unit, colors, yAxisBounds]);

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

// =============================================================================
// Aggregated Dose-Response Chart (% Inhibition with Error Bars)
// =============================================================================

export interface NormalizedSeries {
  concentrations: number[];
  percentInhibition: number[];
}

/**
 * Normalize a single data series to % inhibition using its fitted min/max.
 * 0% = no inhibition (response at maxVal), 100% = full inhibition (response at minVal).
 */
export function normalizeToPercentInhibition(
  responses: number[],
  minVal: number,
  maxVal: number
): number[] {
  const range = maxVal - minVal;
  if (range === 0) return responses.map(() => 0);
  return responses.map(r => 100 * (maxVal - r) / range);
}

/**
 * Compute mean and standard deviation for an array of numbers
 */
function computeStats(values: number[]): { mean: number; sem: number; n: number } {
  const n = values.length;
  if (n === 0) return { mean: 0, sem: 0, n: 0 };
  const mean = values.reduce((a, b) => a + b, 0) / n;
  if (n === 1) return { mean, sem: 0, n: 1 };
  const variance = values.reduce((sum, v) => sum + (v - mean) ** 2, 0) / (n - 1);
  const stdev = Math.sqrt(variance);
  return { mean, sem: stdev / Math.sqrt(n), n };
}

/**
 * Chart.js plugin that draws vertical error bars on scatter points.
 * Expects dataset.errorBars: { low: number; high: number }[] parallel to data.
 */
const errorBarPlugin = {
  id: 'errorBars',
  afterDatasetsDraw(chart: any) {
    const { ctx } = chart;
    chart.data.datasets.forEach((dataset: any, datasetIndex: number) => {
      if (!dataset.errorBars) return;
      const meta = chart.getDatasetMeta(datasetIndex);
      if (meta.hidden) return;

      ctx.save();
      ctx.strokeStyle = dataset.errorBarColor || dataset.borderColor || 'rgba(0,0,0,0.5)';
      ctx.lineWidth = 1.5;
      const capWidth = 4;

      meta.data.forEach((point: any, index: number) => {
        const errorBar = dataset.errorBars[index];
        if (!errorBar) return;

        const yScale = chart.scales.y;
        const yLow = yScale.getPixelForValue(errorBar.low);
        const yHigh = yScale.getPixelForValue(errorBar.high);
        const x = point.x;

        // Vertical line
        ctx.beginPath();
        ctx.moveTo(x, yLow);
        ctx.lineTo(x, yHigh);
        ctx.stroke();

        // Top cap
        ctx.beginPath();
        ctx.moveTo(x - capWidth, yHigh);
        ctx.lineTo(x + capWidth, yHigh);
        ctx.stroke();

        // Bottom cap
        ctx.beginPath();
        ctx.moveTo(x - capWidth, yLow);
        ctx.lineTo(x + capWidth, yLow);
        ctx.stroke();
      });

      ctx.restore();
    });
  },
};

interface AggregatedDoseResponseChartProps {
  /** All data series to aggregate (should share the same dilution series) */
  seriesData: {
    concentrations: number[];
    responses: number[];
    minVal: number;
    maxVal: number;
  }[];
  unit?: string;
  title?: string;
  width?: number;
  height?: number;
}

export function AggregatedDoseResponseChart({
  seriesData,
  unit = 'nM',
  title,
  width = 430,
  height = 350,
}: AggregatedDoseResponseChartProps) {
  const theme = useTheme();
  const isDarkMode = theme.palette.mode === 'dark';

  const colors = useMemo(() => ({
    meanPoint: isDarkMode ? 'rgba(100, 181, 246, 0.9)' : 'rgba(70, 130, 180, 0.8)',
    meanPointBorder: isDarkMode ? 'rgba(100, 181, 246, 1)' : 'rgba(70, 130, 180, 1)',
    errorBar: isDarkMode ? 'rgba(100, 181, 246, 0.6)' : 'rgba(70, 130, 180, 0.5)',
    individualLine: isDarkMode ? 'rgba(255, 255, 255, 0.12)' : 'rgba(0, 0, 0, 0.08)',
    textColor: theme.palette.text.primary,
    secondaryText: theme.palette.text.secondary,
    gridColor: isDarkMode ? 'rgba(255, 255, 255, 0.1)' : 'rgba(0, 0, 0, 0.1)',
  }), [isDarkMode, theme.palette]);

  // Normalize all series and compute per-concentration stats
  const { concentrations, stats, normalizedSeries } = useMemo(() => {
    if (seriesData.length === 0) {
      return { concentrations: [], stats: [], normalizedSeries: [] };
    }

    // Use concentrations from the first series (all should be the same)
    const concs = seriesData[0].concentrations;

    // Normalize each series
    const normalized = seriesData.map(s =>
      normalizeToPercentInhibition(s.responses, s.minVal, s.maxVal)
    );

    // Compute stats at each concentration
    const perConcStats = concs.map((_, concIdx) => {
      const values = normalized
        .map(series => series[concIdx])
        .filter(v => v !== undefined && !isNaN(v));
      return computeStats(values);
    });

    return { concentrations: concs, stats: perConcStats, normalizedSeries: normalized };
  }, [seriesData]);

  const chartData = useMemo<ChartData<'scatter'>>(() => {
    if (concentrations.length === 0) return { datasets: [] };

    const meanPoints = concentrations
      .map((conc, idx) => ({ x: conc, y: stats[idx].mean }))
      .filter(p => p.x > 0);

    const errorBars = concentrations
      .map((conc, idx) => {
        if (conc <= 0) return null;
        const { mean, sem } = stats[idx];
        return { low: mean - sem, high: mean + sem };
      })
      .filter(Boolean);

    const datasets: ChartData<'scatter'>['datasets'] = [];

    // Individual series as faint lines for context
    normalizedSeries.forEach((series, seriesIdx) => {
      const points = concentrations
        .map((conc, idx) => ({ x: conc, y: series[idx] }))
        .filter(p => p.x > 0 && !isNaN(p.y));

      datasets.push({
        type: 'line' as const,
        label: seriesIdx === 0 ? 'Individual Series' : `_series_${seriesIdx}`,
        data: points,
        backgroundColor: 'transparent',
        borderColor: colors.individualLine,
        borderWidth: 1,
        pointRadius: 0,
        fill: false,
      } as any);
    });

    // Mean points with error bars
    datasets.push({
      label: `Mean (n=${seriesData.length})`,
      data: meanPoints,
      backgroundColor: colors.meanPoint,
      borderColor: colors.meanPointBorder,
      pointRadius: 6,
      pointHoverRadius: 8,
      errorBars,
      errorBarColor: colors.errorBar,
    } as any);

    return { datasets };
  }, [concentrations, stats, normalizedSeries, seriesData.length, colors]);

  const options = useMemo<ChartOptions<'scatter'>>(() => ({
    responsive: true,
    maintainAspectRatio: false,
    animation: false,
    plugins: {
      legend: {
        display: true,
        position: 'bottom' as const,
        labels: {
          usePointStyle: true,
          color: colors.textColor,
          filter: (item) => !item.text?.startsWith('_'),
        },
      },
      title: {
        display: !!title,
        text: title || '',
        font: { size: 14 },
        color: colors.textColor,
      },
      tooltip: {
        callbacks: {
          label: (context) => {
            const point = context.raw as { x: number; y: number };
            const idx = context.dataIndex;
            const stat = stats[idx];
            if (stat && stat.n > 1) {
              return `[${point.x.toExponential(2)}]: ${point.y.toFixed(1)}% ± ${stat.sem.toFixed(1)}%`;
            }
            return `[${point.x.toExponential(2)}]: ${point.y.toFixed(1)}%`;
          },
        },
      },
    },
    scales: {
      x: {
        type: 'logarithmic' as const,
        title: {
          display: true,
          text: `Concentration (${unit})`,
          color: colors.textColor,
        },
        ticks: {
          color: colors.secondaryText,
          callback: (value) => {
            const num = Number(value);
            if (num >= 1000) return `${num / 1000}k`;
            if (num >= 1) return num.toString();
            return num.toExponential(0);
          },
        },
        grid: { color: colors.gridColor },
      },
      y: {
        type: 'linear' as const,
        title: {
          display: true,
          text: '% Inhibition',
          color: colors.textColor,
        },
        ticks: { color: colors.secondaryText },
        grid: { color: colors.gridColor },
      },
    },
  }), [title, unit, colors, stats]);

  if (concentrations.length === 0) {
    return (
      <Paper sx={{ p: 2, width, height, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
        <Typography color="text.secondary">No data to aggregate</Typography>
      </Paper>
    );
  }

  return (
    <Box sx={{ width, height }}>
      <Scatter data={chartData} options={options} plugins={[errorBarPlugin]} />
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
  const theme = useTheme();
  const isDarkMode = theme.palette.mode === 'dark';

  // Theme-aware colors for thumbnail
  const colors = useMemo(() => ({
    dataPoint: isDarkMode ? 'rgba(100, 181, 246, 0.9)' : 'rgba(70, 130, 180, 0.8)',
    dataPointBorder: isDarkMode ? 'rgba(100, 181, 246, 1)' : 'rgba(70, 130, 180, 1)',
    validCurve: isDarkMode ? 'rgba(239, 83, 80, 1)' : 'rgba(220, 53, 69, 1)',
    invalidCurve: isDarkMode ? 'rgba(150, 150, 150, 0.9)' : 'rgba(100, 100, 100, 0.9)',
    tickColor: isDarkMode ? 'rgba(255, 255, 255, 0.5)' : 'rgba(0, 0, 0, 0.5)',
    borderColor: isDarkMode ? 'rgba(255, 255, 255, 0.1)' : 'rgba(0, 0, 0, 0.1)',
    emptyBg: isDarkMode ? 'grey.800' : 'grey.100',
  }), [isDarkMode]);

  const chartData = useMemo<ChartData<'scatter'>>(() => {
    const { concentrations, responses } = data;

    const points = concentrations
      .map((conc, idx) => ({ x: conc, y: responses[idx] }))
      .filter(p => p.x > 0 && p.y !== undefined);

    const datasets: ChartData<'scatter'>['datasets'] = [
      {
        label: 'Data',
        data: points,
        backgroundColor: colors.dataPoint,
        borderColor: colors.dataPointBorder,
        pointRadius: 3,
      },
    ];

    // Add fitted curve - prefer backend curve_points for algorithm-agnostic rendering
    const hasBackendCurve = fit?.curvePoints && Array.isArray(fit.curvePoints) && fit.curvePoints.length > 0;

    // Fallback: generate curve using 4PL equation if no backend curve and we have valid 4PL params
    const hasValid4PLParams =
      !hasBackendCurve &&
      fit?.ec50 != null &&
      fit.ec50 > 0 &&
      fit.hill != null &&
      typeof fit.hill === 'number' &&
      fit.minVal != null &&
      typeof fit.minVal === 'number' &&
      fit.maxVal != null &&
      typeof fit.maxVal === 'number';

    if (hasBackendCurve || hasValid4PLParams) {
      let curvePoints: { x: number; y: number }[];

      if (hasBackendCurve) {
        // Use backend-generated curve points (works for all algorithms)
        curvePoints = fit!.curvePoints!.map(([x, y]) => ({ x, y }));
      } else {
        // Fallback: generate using 4PL equation
        let normalizedMin = fit!.minVal!;
        let normalizedMax = fit!.maxVal!;
        if (normalizedMin > normalizedMax) {
          [normalizedMin, normalizedMax] = [normalizedMax, normalizedMin];
        }

        curvePoints = generateFittedCurve(
          concentrations,
          fit!.ec50!,
          fit!.hill!,
          normalizedMin,
          normalizedMax,
          50
        );
      }

      // Use red for valid fits, gray for invalid
      const curveColor = fit?.status === 'valid'
        ? colors.validCurve
        : colors.invalidCurve;

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
  }, [data, fit, colors]);

  // Calculate axis bounds from data, fitted asymptotes, AND control values.
  // Including controls ensures flat/inactive data looks flat in the thumbnail.
  const axisBounds = useMemo(() => {
    const { concentrations, responses } = data;
    const validConcs = concentrations.filter(c => c > 0);
    const validResponses = responses.filter(r => r !== undefined && r !== null) as number[];

    let yMin = validResponses.length ? Math.min(...validResponses) : 0;
    let yMax = validResponses.length ? Math.max(...validResponses) : 100;

    // Include fitted asymptotes if available
    if (fit?.minVal != null && typeof fit.minVal === 'number') {
      yMin = Math.min(yMin, fit.minVal);
    }
    if (fit?.maxVal != null && typeof fit.maxVal === 'number') {
      yMax = Math.max(yMax, fit.maxVal);
    }

    // Include observed control values so the axis reflects the full assay range
    if (data.minControlResponse != null) {
      yMin = Math.min(yMin, data.minControlResponse);
      yMax = Math.max(yMax, data.minControlResponse);
    }
    if (data.maxControlResponse != null) {
      yMin = Math.min(yMin, data.maxControlResponse);
      yMax = Math.max(yMax, data.maxControlResponse);
    }

    // Add 5% padding
    const yRange = yMax - yMin;
    const yPadding = yRange * 0.05;

    return {
      xMin: validConcs.length ? Math.min(...validConcs) : 1,
      xMax: validConcs.length ? Math.max(...validConcs) : 1000,
      yMin: yMin - yPadding,
      yMax: yMax + yPadding,
    };
  }, [data, fit]);

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
        border: { display: true, color: colors.borderColor },
        ticks: {
          display: true,
          font: { size: 8 },
          color: colors.tickColor,
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
        min: axisBounds.yMin,
        max: axisBounds.yMax,
        grid: { display: false },
        border: { display: true, color: colors.borderColor },
        ticks: {
          display: true,
          font: { size: 8 },
          color: colors.tickColor,
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
  }), [colors, axisBounds]);

  if (!data.concentrations.length) {
    return (
      <Box
        sx={{
          width: size,
          height: size,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          bgcolor: colors.emptyBg,
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
