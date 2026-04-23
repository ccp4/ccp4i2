'use client';

import { useMemo } from 'react';
import { Box, Typography } from '@mui/material';
import {
  Chart as ChartJS,
  RadialLinearScale,
  PointElement,
  LineElement,
  Filler,
  Tooltip as ChartTooltip,
  Legend,
  type ChartData,
  type ChartOptions,
  type TooltipItem,
} from 'chart.js';
import { Radar } from 'react-chartjs-2';
import type { CompactRow, ProtocolInfo } from '@/types/compounds/aggregation';
import type { ScorecardConfig } from '@/types/compounds/models';
import { evaluateScorecard } from '@/lib/compounds/scorecard';

ChartJS.register(RadialLinearScale, PointElement, LineElement, Filler, ChartTooltip, Legend);

interface Props {
  config: ScorecardConfig | null | undefined;
  compound: CompactRow;
  /** For richer tooltips — resolve protocol names / units on hover. */
  protocols?: ProtocolInfo[];
  /** `small` for card use (140px, no axis labels); `large` for preview (400px). */
  size?: 'small' | 'large';
}

/**
 * Radar chart of a compound's per-axis scorecard values. Outward = better
 * (axis extent = target_value, centre = poor_value). Axes with missing data
 * render as gaps so a missing measurement doesn't read as a failing one.
 */
export function CompoundSpider({ config, compound, size = 'small' }: Props) {
  const evals = useMemo(() => evaluateScorecard(config, compound), [config, compound]);

  const data: ChartData<'radar'> = useMemo(
    () => ({
      labels: evals.map(({ axis }) => axis.label || '(unnamed)'),
      datasets: [
        {
          label: compound.formatted_id,
          // Coerce missing axes to 0 rather than null — Chart.js 4 radar
          // can fail to construct the polygon path when any datapoint is
          // null, which manifests as a white-screen client error. The
          // visual cost: a missing measurement reads as a pinched-to-centre
          // (i.e. poor) axis. Tooltip still says "no data" for clarity.
          data: evals.map(({ t }) => t ?? 0),
          backgroundColor: 'rgba(25, 118, 210, 0.25)',
          borderColor: 'rgba(25, 118, 210, 1)',
          borderWidth: 1.5,
          pointRadius: size === 'small' ? 2 : 3,
          pointBackgroundColor: evals.map(({ t }) =>
            t == null ? 'rgba(160, 160, 160, 0.8)' : 'rgba(25, 118, 210, 1)',
          ),
          pointBorderColor: '#fff',
          pointBorderWidth: 1,
        },
      ],
    }),
    [evals, compound.formatted_id, size],
  );

  const options: ChartOptions<'radar'> = useMemo(() => {
    const isSmall = size === 'small';
    // Small spiders live in a card next to the structure, so labels have to
    // be compact. Truncate to ~10 chars so typical labels like
    // "TM potency" or "Lipinski" fit without wrapping, but longer ones like
    // "BAF3 mutant selectivity" still show enough to disambiguate.
    const truncate = (label: string, max: number) =>
      label.length > max ? label.slice(0, max - 1) + '…' : label;
    return {
      responsive: true,
      maintainAspectRatio: true,
      // Don't bridge across nulls — missing data stays as a visible gap.
      spanGaps: false,
      // Extra space around the chart so pointLabels don't clip at the edges
      // of the container (especially for the small size).
      layout: { padding: isSmall ? 4 : 8 },
      scales: {
        r: {
          min: 0,
          max: 1,
          ticks: { display: false, stepSize: 0.25 },
          pointLabels: {
            display: true,
            font: { size: isSmall ? 9 : 11 },
            color: 'rgba(0, 0, 0, 0.75)',
            padding: isSmall ? 4 : 8,
            callback: (value: string) =>
              truncate(String(value), isSmall ? 12 : 32),
          },
          grid: { circular: true },
          angleLines: { color: 'rgba(0, 0, 0, 0.1)' },
        },
      },
      plugins: {
        legend: { display: false },
        tooltip: {
          // Enabled for both sizes now — with labels showing, a hover
          // tooltip is still useful for reading the raw value / tier and
          // for the full (untruncated) axis label.
          enabled: true,
          callbacks: {
            title: (items: TooltipItem<'radar'>[]) => {
              const i = items[0]?.dataIndex;
              if (i == null) return '';
              return evals[i]?.axis.label ?? '';
            },
            label: (ctx: TooltipItem<'radar'>) => {
              const ev = evals[ctx.dataIndex];
              if (!ev) return '';
              if (ev.value == null) return 'no data';
              const formatted = formatValue(ev.value);
              const tier = tierLabel(ev.t);
              return `${formatted}  ${tier}`;
            },
          },
        },
      },
    };
  }, [evals, size]);

  // Bump the small size to accommodate axis labels around the perimeter
  // without squashing the plot area.
  const dimension = size === 'small' ? 180 : 360;

  if (!config?.axes?.length) {
    if (size === 'small') return null;
    return (
      <Box
        sx={{
          width: dimension,
          height: dimension,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          border: 1,
          borderColor: 'divider',
          borderRadius: 1,
          p: 2,
        }}
      >
        <Typography variant="caption" color="text.secondary" align="center">
          Add axes to see a preview.
        </Typography>
      </Box>
    );
  }

  return (
    <Box sx={{ width: dimension, height: dimension }}>
      <Radar data={data} options={options} />
    </Box>
  );
}

function formatValue(v: number): string {
  if (!Number.isFinite(v)) return '—';
  if (Math.abs(v) >= 100 || Math.abs(v) < 0.1) return v.toPrecision(3);
  return v.toFixed(2);
}

function tierLabel(t: number | null): string {
  if (t == null) return '';
  if (t >= 1) return '(excellent)';
  if (t >= 2 / 3) return '(good)';
  if (t >= 1 / 3) return '(mid)';
  if (t > 0) return '(poor)';
  return '(failing)';
}
