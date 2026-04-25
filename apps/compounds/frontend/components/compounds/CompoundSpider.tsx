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
  type Chart,
  type ChartData,
  type ChartOptions,
  type Plugin,
  type TooltipItem,
} from 'chart.js';
import { Radar } from 'react-chartjs-2';
import type { CompactRow, ProtocolInfo } from '@/types/compounds/aggregation';
import type { ScorecardConfig } from '@/types/compounds/models';
import {
  evaluateScorecard,
  groupAxesBySector,
  sectorColour,
  type AxisEvaluation,
} from '@/lib/compounds/scorecard';

ChartJS.register(RadialLinearScale, PointElement, LineElement, Filler, ChartTooltip, Legend);

// ---------------------------------------------------------------------------
// Sector wedge plugin — paints translucent coloured pie slices behind the
// polygon for each sector, fading from saturated at the rim to transparent
// near the centre. Sector labels (large size only) sit just outside the
// per-axis label ring at each wedge's angular midpoint.
//
// Chart.js radar puts axis 0 at -π/2 (top) and steps clockwise by 2π/N.
// Wedge boundaries are the midpoint angles between adjacent axes. Once we
// have boundaries, each sector spans from the midpoint-angle BEFORE its
// first axis to the midpoint-angle AFTER its last axis.
// ---------------------------------------------------------------------------

interface WedgeSpec {
  startAngle: number;
  endAngle: number;
  colour: string;
}

const sectorWedgePlugin: Plugin<'radar'> = {
  id: 'sectorWedges',
  beforeDatasetsDraw(chart: Chart<'radar'>) {
    const opts = (chart.options.plugins as Record<string, unknown> | undefined)
      ?.sectorWedges as { wedges?: WedgeSpec[] } | undefined;
    if (!opts?.wedges?.length) return;
    const r = chart.scales.r as unknown as {
      xCenter: number;
      yCenter: number;
      drawingArea: number;
    };
    if (!r || typeof r.drawingArea !== 'number') return;
    const { ctx } = chart;
    const cx = r.xCenter;
    const cy = r.yCenter;
    const outer = r.drawingArea;
    for (const wedge of opts.wedges) {
      const grad = ctx.createRadialGradient(cx, cy, outer * 0.15, cx, cy, outer);
      grad.addColorStop(0, hexToRgba(wedge.colour, 0));
      grad.addColorStop(0.85, hexToRgba(wedge.colour, 0.22));
      grad.addColorStop(1, hexToRgba(wedge.colour, 0.35));
      ctx.save();
      ctx.beginPath();
      ctx.moveTo(cx, cy);
      ctx.arc(cx, cy, outer, wedge.startAngle, wedge.endAngle);
      ctx.closePath();
      ctx.fillStyle = grad;
      ctx.fill();
      ctx.restore();
    }
  },
};

ChartJS.register(sectorWedgePlugin);

function hexToRgba(hex: string, alpha: number): string {
  const m = hex.replace('#', '');
  const bigint = parseInt(m.length === 3 ? m.split('').map((c) => c + c).join('') : m, 16);
  const r = (bigint >> 16) & 255;
  const g = (bigint >> 8) & 255;
  const b = bigint & 255;
  return `rgba(${r}, ${g}, ${b}, ${alpha})`;
}

function buildWedges(
  evals: ReadonlyArray<AxisEvaluation>,
): WedgeSpec[] {
  const N = evals.length;
  if (N === 0) return [];
  // Chart.js radar: axis i is at angle -π/2 + i * 2π/N (top, clockwise).
  const axisAngle = (i: number) => -Math.PI / 2 + (2 * Math.PI * i) / N;
  // Midpoint between axis i and axis i+1 (mod N).
  const boundary = (i: number) => axisAngle(i) + Math.PI / N;

  // Collect contiguous index ranges of the same sector. Because the axes
  // are already sorted by sector when passed in, this is a single pass.
  const groups: Array<{ sector: string | null; first: number; last: number }> = [];
  for (let i = 0; i < N; i++) {
    const sector = evals[i].axis.sector?.trim().toLowerCase() || null;
    const tail = groups[groups.length - 1];
    if (tail && tail.sector === sector) {
      tail.last = i;
    } else {
      groups.push({ sector, first: i, last: i });
    }
  }

  return groups.map((g) => ({
    startAngle: boundary(g.first - 1),
    endAngle: boundary(g.last),
    colour: g.sector ? sectorColour(g.sector) : '#bdbdbd',
  }));
}

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
  // Evaluate axes in their original (config) order, then reorder by sector
  // so adjacent radar axes share a sector. The wedge plugin reads the same
  // sector ordering off the resulting evals array.
  const evals = useMemo(() => {
    const raw = evaluateScorecard(config, compound);
    if (raw.length === 0) return raw;
    const groups = groupAxesBySector(raw.map((e) => e.axis));
    return groups.flatMap((g) => g.items).map(({ index }) => raw[index]);
  }, [config, compound]);

  const wedges = useMemo(() => buildWedges(evals), [evals]);

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
    return {
      responsive: true,
      maintainAspectRatio: true,
      // Don't bridge across nulls — missing data stays as a visible gap.
      spanGaps: false,
      // Generous canvas padding so axis labels don't clip on the
      // left/right (long labels overhanging the polygon's axis end) or
      // top/bottom (now that labels can wrap to two lines, the second
      // line needs vertical room).
      layout: { padding: isSmall ? 32 : 20 },
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
            // Wrap multi-word labels to up to 2 lines; truncate single
            // long words. Chart.js renders array values on separate lines.
            callback: (value: string | number) =>
              wrapOrTruncate(String(value), isSmall ? 12 : 32),
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
        // Sector wedges read by our custom plugin (registered above).
        // Labels were dropped — wedge colour alone communicates grouping
        // and avoids overlap with per-axis labels at the rim.
        sectorWedges: { wedges },
      } as ChartOptions<'radar'>['plugins'],
    };
  }, [evals, size, wedges]);

  // Outer container size. Includes the layout.padding above, so the
  // polygon itself ends up ≈ (dimension - 2 × padding) wide.
  const dimension = size === 'small' ? 220 : 380;

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

/**
 * Format an axis label for display in the radar's perimeter.
 *
 *  - Returns the original string if it already fits (≤ maxLen chars).
 *  - For multi-word labels longer than maxLen, packs words greedily
 *    into a first line that fits, dumps the rest into a second line,
 *    truncates that second line if still too long. Chart.js renders
 *    array-valued pointLabel callbacks across multiple lines.
 *  - For a single long word with no space to break on, falls back to
 *    a one-line ellipsis truncation.
 */
function wrapOrTruncate(label: string, maxLen: number): string | string[] {
  if (label.length <= maxLen) return label;

  const words = label.split(/\s+/).filter(Boolean);
  if (words.length < 2) {
    return label.slice(0, maxLen - 1) + '…';
  }

  let line1 = words[0];
  let i = 1;
  while (i < words.length && (line1 + ' ' + words[i]).length <= maxLen) {
    line1 += ' ' + words[i];
    i++;
  }
  let line2 = words.slice(i).join(' ');
  if (!line2) return line1;
  if (line2.length > maxLen) {
    line2 = line2.slice(0, maxLen - 1) + '…';
  }
  return [line1, line2];
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
