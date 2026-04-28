'use client';

import {
  useState,
  useMemo,
  useRef,
  useEffect,
  useCallback,
  forwardRef,
  useImperativeHandle,
} from 'react';
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
import { ScatterCategorisationStrip } from './aggregation/ScatterCategorisationStrip';
import { SCATTER_CATEGORY_COLOURS } from '@/lib/compounds/scatter-palette';
import { ListSubheader } from '@mui/material';
import {
  CompactRow,
  ProtocolInfo,
  MolecularPropertyName,
} from '@/types/compounds/aggregation';

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

export interface CategorisationSelection {
  id: string;
  name: string;
  /** Set for O(1) membership tests on every point. */
  compoundIds: ReadonlySet<string>;
}

interface ProtocolScatterPlotProps {
  /** Compact aggregation data rows */
  data: CompactRow[];
  /** Available protocols with geomean values (empty if geomean not requested) */
  protocols: ProtocolInfo[];
  /** Molecular properties present in the data (from meta.include_properties) */
  includedProperties?: MolecularPropertyName[];
  /**
   * 'dialog' (default, back-compat) — renders a "Scatter Plot" trigger button
   *   plus a Dialog that opens to show the chart. Used by the row-pencil
   *   "open scatter for this protocol" flow on cards / pivot views.
   * 'inline' — renders the chart UI directly with no Button or Dialog
   *   wrapper. Used when format=scatter is the page's primary view mode.
   */
  mode?: 'dialog' | 'inline';
  /**
   * When non-empty, points are partitioned into colour groups by
   * membership: one dataset per category (in order), plus an "Other"
   * dataset for compounds in none. Each category's `compoundIds` set
   * decides who lives in that colour group. Driven by the chemotype
   * chip strip via `?colour_by=<scaffold-name>,<scaffold-name>`.
   */
  categorisationSelections?: ReadonlyArray<CategorisationSelection>;
  /**
   * Per-scaffold-name → matching compound formatted_ids. The chip strip
   * uses this for relevance-filtering and match-count badges; the chart
   * uses it indirectly via the `categorisationSelections` derived above.
   */
  scaffoldMatches?: ReadonlyMap<string, ReadonlySet<string>>;
}

type AxisKind = 'protocol' | 'property';

interface AxisConfig {
  kind: AxisKind;
  /** Protocol UUID when kind='protocol', MolecularPropertyName when kind='property' */
  key: string | null;
  scale: 'linear' | 'logarithmic';
  min: string;
  max: string;
}

/** Display labels for molecular properties (matches PredicateBuilder) */
const PROPERTY_LABELS: Record<MolecularPropertyName, { label: string; description: string }> = {
  molecular_weight: { label: 'MW', description: 'Molecular Weight' },
  heavy_atom_count: { label: 'HA', description: 'Heavy Atom Count' },
  hbd: { label: 'HBD', description: 'Hydrogen Bond Donors' },
  hba: { label: 'HBA', description: 'Hydrogen Bond Acceptors' },
  clogp: { label: 'cLogP', description: 'Calculated LogP' },
  tpsa: { label: 'TPSA', description: 'Topological Polar Surface Area' },
  rotatable_bonds: { label: 'RotB', description: 'Rotatable Bonds' },
  fraction_sp3: { label: 'Fsp3', description: 'Fraction sp3 Carbons' },
};

/** Encode an axis selection as a single string for the Select value */
function encodeAxisValue(kind: AxisKind, key: string): string {
  return `${kind}:${key}`;
}

function decodeAxisValue(value: string): { kind: AxisKind; key: string } | null {
  const idx = value.indexOf(':');
  if (idx < 0) return null;
  const kind = value.slice(0, idx) as AxisKind;
  const key = value.slice(idx + 1);
  if (kind !== 'protocol' && kind !== 'property') return null;
  return { kind, key };
}

/** Sensible default scale for a given axis kind */
function defaultScale(kind: AxisKind): 'linear' | 'logarithmic' {
  return kind === 'protocol' ? 'logarithmic' : 'linear';
}

/**
 * Snap an axis extremum to a "nice" number that renders cleanly on tick labels.
 *
 * Log scale: snap to the nearest 1, 2, or 5 × 10ⁿ boundary that widens the axis
 * (floor below for min, ceil above for max), giving 1 / 2 / 5 / 10 / 20 / 50 etc.
 *
 * Linear scale: snap to one significant figure at the range's natural step size,
 * giving 10 / 20 / 50 / 100 rather than 8.734 or 468.57.
 *
 * `otherExtreme` is the opposite extreme of the same data series, used only to
 * pick a sensible step size when the value itself is zero or the range is zero.
 */
function niceBound(
  value: number,
  direction: 'down' | 'up',
  scale: 'linear' | 'logarithmic',
  otherExtreme: number,
): number {
  if (!isFinite(value)) return direction === 'down' ? 0 : 1;

  if (scale === 'logarithmic') {
    if (value <= 0) return direction === 'down' ? 0.001 : 1000;
    const exp = Math.floor(Math.log10(value));
    const mantissa = value / Math.pow(10, exp);
    const stops = [1, 2, 5, 10];
    if (direction === 'down') {
      // Largest nice ≤ mantissa
      let best = 1;
      for (const m of stops) if (m <= mantissa) best = m;
      return best === 10 ? Math.pow(10, exp + 1) : best * Math.pow(10, exp);
    } else {
      // Smallest nice ≥ mantissa
      for (const m of stops) if (m >= mantissa) return m * Math.pow(10, exp);
      return 10 * Math.pow(10, exp);
    }
  }

  // Linear: round to 1 sig fig at the natural step of the range.
  const range = Math.max(Math.abs(value - otherExtreme), Math.abs(value));
  if (range === 0) {
    return direction === 'down' ? value - 1 : value + 1;
  }
  const stepExp = Math.floor(Math.log10(range)) - 1;
  const step = Math.pow(10, stepExp);
  return direction === 'down'
    ? Math.floor(value / step) * step
    : Math.ceil(value / step) * step;
}

export interface ProtocolScatterPlotHandle {
  /**
   * Open the scatter dialog with a specific protocol forced on the X axis.
   * Y axis is preserved if it's already a different protocol; otherwise the
   * first protocol that differs from `protocolId` is chosen.
   */
  openForProtocol(protocolId: string): void;
}

/**
 * Scatter plot component for comparing geomean values between two protocols.
 * Each point represents a compound, with coordinates being the geomean values
 * for the selected X and Y protocols.
 */
export const ProtocolScatterPlot = forwardRef<
  ProtocolScatterPlotHandle,
  ProtocolScatterPlotProps
>(function ProtocolScatterPlot({
  data,
  protocols,
  includedProperties = [],
  mode = 'dialog',
  categorisationSelections,
  scaffoldMatches,
}, ref) {
  const router = useRouter();
  const chartRef = useRef<ChartJS<'scatter'>>(null);

  // Build the ordered list of axis candidates: protocols first, then properties.
  const axisOptions = useMemo(() => {
    const options: Array<{ kind: AxisKind; key: string; label: string; description?: string }> = [];
    for (const p of protocols) {
      options.push({ kind: 'protocol', key: p.id, label: p.name });
    }
    for (const propName of includedProperties) {
      const meta = PROPERTY_LABELS[propName];
      options.push({
        kind: 'property',
        key: propName,
        label: meta?.label ?? propName,
        description: meta?.description,
      });
    }
    return options;
  }, [protocols, includedProperties]);

  const [open, setOpen] = useState(false);
  // Pick the first two available axis candidates as defaults
  const defaultX = axisOptions[0];
  const defaultY = axisOptions[1] ?? axisOptions[0];
  const [xAxis, setXAxis] = useState<AxisConfig>({
    kind: defaultX?.kind ?? 'protocol',
    key: defaultX?.key ?? null,
    scale: defaultScale(defaultX?.kind ?? 'protocol'),
    min: '',
    max: '',
  });
  const [yAxis, setYAxis] = useState<AxisConfig>({
    kind: defaultY?.kind ?? 'protocol',
    key: defaultY?.key ?? null,
    scale: defaultScale(defaultY?.kind ?? 'protocol'),
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
  }, [xAxis.key, yAxis.key]);

  // Seed axis selections once options become available (e.g. after data loads)
  useEffect(() => {
    if (axisOptions.length === 0) return;
    if (!xAxis.key) {
      const seed = axisOptions[0];
      setXAxis((prev) => ({
        ...prev,
        kind: seed.kind,
        key: seed.key,
        scale: defaultScale(seed.kind),
      }));
    }
    if (!yAxis.key) {
      const seed = axisOptions[1] ?? axisOptions[0];
      setYAxis((prev) => ({
        ...prev,
        kind: seed.kind,
        key: seed.key,
        scale: defaultScale(seed.kind),
      }));
    }
  }, [axisOptions, xAxis.key, yAxis.key]);

  // Imperative entry point: lets callers (e.g. a plot icon on a protocol
  // header) open the dialog with a forced X axis. Y is preserved when it's
  // already a different protocol; otherwise we pick the first protocol that
  // differs from the forced X so the plot is immediately informative.
  useImperativeHandle(
    ref,
    () => ({
      openForProtocol(protocolId: string) {
        setXAxis((prev) => ({
          ...prev,
          kind: 'protocol',
          key: protocolId,
          scale: 'logarithmic',
        }));
        setYAxis((prev) => {
          const yIsOtherProtocol =
            prev.kind === 'protocol' && prev.key && prev.key !== protocolId;
          if (yIsOtherProtocol) return prev;
          const altProtocol = protocols.find((p) => p.id !== protocolId);
          if (!altProtocol) return prev;
          return {
            ...prev,
            kind: 'protocol',
            key: altProtocol.id,
            scale: 'logarithmic',
          };
        });
        setOpen(true);
      },
    }),
    [protocols],
  );

  // Axis display names (protocol name or property label)
  const xAxisLabel = useMemo(() => {
    const opt = axisOptions.find((o) => o.kind === xAxis.kind && o.key === xAxis.key);
    return opt?.label || 'Select axis';
  }, [axisOptions, xAxis.kind, xAxis.key]);
  const yAxisLabel = useMemo(() => {
    const opt = axisOptions.find((o) => o.kind === yAxis.kind && o.key === yAxis.key);
    return opt?.label || 'Select axis';
  }, [axisOptions, yAxis.kind, yAxis.key]);

  const handleAxisChange = useCallback(
    (setter: React.Dispatch<React.SetStateAction<AxisConfig>>, encoded: string) => {
      const decoded = decodeAxisValue(encoded);
      if (!decoded) return;
      setter((prev) => ({
        ...prev,
        kind: decoded.kind,
        key: decoded.key,
        // When the axis kind changes, reset to the sensible default scale
        // for that kind (log for protocols, linear for properties).
        scale: prev.kind === decoded.kind ? prev.scale : defaultScale(decoded.kind),
        // Also clear any manual min/max since they will be on a different scale
        min: prev.kind === decoded.kind ? prev.min : '',
        max: prev.kind === decoded.kind ? prev.max : '',
      }));
    },
    []
  );

  /** Extract the numeric value for a row on a given axis, or null if missing. */
  const getAxisValue = useCallback((row: CompactRow, axis: AxisConfig): number | null => {
    if (!axis.key) return null;
    if (axis.kind === 'protocol') {
      const v = row.protocols[axis.key]?.geomean;
      return v ?? null;
    }
    // property
    const v = row.properties?.[axis.key as MolecularPropertyName];
    return v ?? null;
  }, []);

  // Build scatter plot data
  const chartData: ChartData<'scatter'> = useMemo(() => {
    if (!xAxis.key || !yAxis.key) {
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

    // Log scales require strictly positive values; linear axes accept any finite number
    const xRequiresPositive = xAxis.scale === 'logarithmic';
    const yRequiresPositive = yAxis.scale === 'logarithmic';

    for (const row of data) {
      const xValue = getAxisValue(row, xAxis);
      const yValue = getAxisValue(row, yAxis);

      if (xValue == null || yValue == null) continue;
      if (!isFinite(xValue) || !isFinite(yValue)) continue;
      if (xRequiresPositive && xValue <= 0) continue;
      if (yRequiresPositive && yValue <= 0) continue;

      points.push({
        x: xValue,
        y: yValue,
        compound_id: row.compound_id,
        formatted_id: row.formatted_id,
        smiles: row.smiles ?? undefined,
        target_name: row.target_name ?? undefined,
      });
    }

    if (!categorisationSelections || categorisationSelections.length === 0) {
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
    }

    // First-match wins: priority order = the array given.
    type Bucket = { x: number; y: number; compound_id: string; formatted_id: string; smiles?: string; target_name?: string };
    const otherBucket: Bucket[] = [];
    const namedBuckets: Bucket[][] = categorisationSelections.map(() => []);
    for (const p of points) {
      let placed = false;
      for (let i = 0; i < categorisationSelections.length; i++) {
        if (categorisationSelections[i].compoundIds.has(p.formatted_id)) {
          namedBuckets[i].push(p);
          placed = true;
          break;
        }
      }
      if (!placed) otherBucket.push(p);
    }
    return {
      datasets: [
        // "Other" goes first so the legend reads it first AND it draws
        // BENEATH the coloured groups (Chart.js draws datasets in order).
        {
          label: 'Other',
          data: otherBucket,
          backgroundColor: 'rgba(150, 150, 150, 0.45)',
          borderColor: 'rgba(120, 120, 120, 0.8)',
          borderWidth: 1,
          pointRadius: 5,
          pointHoverRadius: 7,
        },
        ...categorisationSelections.map((sel, i) => ({
          label: sel.name || `Selection ${i + 1}`,
          data: namedBuckets[i],
          backgroundColor: SCATTER_CATEGORY_COLOURS[i % SCATTER_CATEGORY_COLOURS.length].fill,
          borderColor: SCATTER_CATEGORY_COLOURS[i % SCATTER_CATEGORY_COLOURS.length].border,
          borderWidth: 1,
          pointRadius: 6,
          pointHoverRadius: 8,
        })),
      ],
      _rawPoints: points,
    };
  }, [data, xAxis, yAxis, getAxisValue, categorisationSelections]);

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

  // Calculate axis bounds, snapped to "nice" numbers so the axis-extreme
  // labels read as 1/2/5/10/20/50/100 rather than 7.9432... or 468.572...
  const axisBounds = useMemo(() => {
    // Read from _rawPoints — dataset[0] is no longer the full set when
    // categorisation splits the points into "Other" + per-selection groups.
    const points = (chartData as any)._rawPoints as Array<{ x: number; y: number }> || [];
    if (points.length === 0) {
      return { xMin: 0.01, xMax: 1000, yMin: 0.01, yMax: 1000 };
    }

    const xValues = points.map((p) => p.x);
    const yValues = points.map((p) => p.y);

    const xMinRaw = Math.min(...xValues);
    const xMaxRaw = Math.max(...xValues);
    const yMinRaw = Math.min(...yValues);
    const yMaxRaw = Math.max(...yValues);

    return {
      xMin: niceBound(xMinRaw, 'down', xAxis.scale, xMaxRaw),
      xMax: niceBound(xMaxRaw, 'up', xAxis.scale, xMinRaw),
      yMin: niceBound(yMinRaw, 'down', yAxis.scale, yMaxRaw),
      yMax: niceBound(yMaxRaw, 'up', yAxis.scale, yMinRaw),
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

  // Isoselectivity reference diagonals.
  //
  // When both axes are protocols on a log scale, a scatter of compound
  // potencies against two assays is a classic selectivity plot. Parallel
  // diagonals at y = m*x for m in {0.01, 0.1, 1, 10, 100} make the
  // "how selective is this compound?" question readable at a glance:
  // the y=x diagonal = non-selective; points 10x above = 10-fold preference
  // for the Y assay; points 10x below = 10-fold preference for the X assay.
  //
  // Not drawn when either axis is linear or a molecular property — the
  // concept doesn't apply in those cases.
  const isoselectivityDatasets = useMemo(() => {
    const bothLog = xAxis.scale === 'logarithmic' && yAxis.scale === 'logarithmic';
    const bothProtocol = xAxis.kind === 'protocol' && yAxis.kind === 'protocol';
    if (!bothLog || !bothProtocol) return [];

    const xMin = xAxis.min ? parseFloat(xAxis.min) : axisBounds.xMin;
    const xMax = xAxis.max ? parseFloat(xAxis.max) : axisBounds.xMax;
    if (!isFinite(xMin) || !isFinite(xMax) || xMin <= 0 || xMax <= 0) return [];

    // Multipliers chosen to cover common pharma selectivity windows
    // (±10x = "meaningful", ±100x = "clean"). Solid line at y=x; dashed at ±.
    const rules: Array<{ m: number; label: string; solid: boolean }> = [
      { m: 100, label: '100× y > x', solid: false },
      { m: 10, label: '10× y > x', solid: false },
      { m: 1, label: 'y = x (non-selective)', solid: true },
      { m: 0.1, label: '10× x > y', solid: false },
      { m: 0.01, label: '100× x > y', solid: false },
    ];

    // Two-point lines are enough in log-log space (straight), but pick
    // the min/max of both axes so the line is trimmed to the plotted area.
    const yMin = yAxis.min ? parseFloat(yAxis.min) : axisBounds.yMin;
    const yMax = yAxis.max ? parseFloat(yAxis.max) : axisBounds.yMax;

    return rules
      .map(({ m, label, solid }) => {
        // The visible segment is the x-range intersected with the y-range/m.
        const xStart = Math.max(xMin, yMin / m);
        const xEnd = Math.min(xMax, yMax / m);
        if (!(xStart > 0 && xEnd > xStart)) return null;
        return {
          label: `iso:${label}`,
          data: [
            { x: xStart, y: m * xStart },
            { x: xEnd, y: m * xEnd },
          ],
          backgroundColor: 'transparent',
          borderColor: solid
            ? 'rgba(120, 120, 120, 0.55)'
            : 'rgba(120, 120, 120, 0.3)',
          borderWidth: solid ? 1.25 : 1,
          borderDash: solid ? [] : [4, 4],
          pointRadius: 0,
          pointHoverRadius: 0,
          showLine: true,
          tension: 0,
          order: 2, // Behind regression (order: 1) and behind points (default 0)
        } as any;
      })
      .filter(Boolean);
  }, [xAxis, yAxis, axisBounds]);

  // Final chart data including regression line and isoselectivity diagonals
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

    for (const ds of isoselectivityDatasets) {
      if (ds) datasets.push(ds);
    }

    return { datasets };
  }, [chartData, showRegression, regressionLineData, isoselectivityDatasets]);

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
            // Skip regression and isoselectivity reference lines
            const dsLabel = dataPoint.dataset.label;
            if (dsLabel === 'Regression Line' || dsLabel?.startsWith('iso:')) {
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
            text: xAxisLabel,
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
            text: yAxisLabel,
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
        if (elements.length === 0) return;
        const element = elements[0];
        // Look up the clicked point in its own dataset; the regression line
        // and isoselectivity overlays don't carry compound_id, so they
        // self-filter via the guard below — no need to special-case index.
        const dataset = finalChartData.datasets[element.datasetIndex];
        const point = (dataset?.data as any[] | undefined)?.[element.index];
        if (point?.compound_id) {
          router.push(`/registry/compounds/${point.compound_id}`);
        }
      },
    }),
    [xAxis, yAxis, xAxisLabel, yAxisLabel, axisBounds, finalChartData, router, setHoveredPoint]
  );

  // Reset axis ranges to auto
  const handleResetRanges = () => {
    setXAxis((prev) => ({ ...prev, min: '', max: '' }));
    setYAxis((prev) => ({ ...prev, min: '', max: '' }));
  };

  const pointCount = ((chartData as any)._rawPoints as unknown[] | undefined)?.length || 0;

  // Only show if we have at least 2 axis candidates to plot against each other
  if (axisOptions.length < 2) {
    return null;
  }

  // The chart UI is identical in both modes; only the wrapper differs.
  const panel = (
    <Box sx={{
      display: 'flex',
      flexDirection: 'column',
      // Inline mode fills the page region we're handed; dialog mode lets
      // DialogContent's own padding/overflow rules drive sizing.
      ...(mode === 'inline' && { flex: 1, minHeight: 0 }),
      px: mode === 'inline' ? 0 : 3,
      py: mode === 'inline' ? 0 : 2.5,
      overflow: 'auto',
    }}>
          {/* Categorisation chip strip — only the inline page-level
              scatter view can usefully colour-by saved selections;
              the row-pencil dialog opens from a single-protocol
              context where colouring isn't part of the workflow. */}
          {mode === 'inline' && (
            <ScatterCategorisationStrip
              matchesByScaffold={scaffoldMatches ?? new Map()}
            />
          )}
          {/* Axis controls */}
          <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap', mb: 2 }}>
            {/* X-Axis controls */}
            <Box sx={{ flex: '1 1 300px', minWidth: 280 }}>
              <Typography variant="caption" color="text.secondary" gutterBottom>
                X-Axis
              </Typography>
              <Box sx={{ display: 'flex', gap: 1, alignItems: 'flex-start' }}>
                <FormControl size="small" sx={{ minWidth: 180 }}>
                  <InputLabel>Axis</InputLabel>
                  <Select
                    value={xAxis.key ? encodeAxisValue(xAxis.kind, xAxis.key) : ''}
                    label="Axis"
                    onChange={(e) => handleAxisChange(setXAxis, e.target.value as string)}
                  >
                    {protocols.length > 0 && (
                      <ListSubheader>Protocols</ListSubheader>
                    )}
                    {protocols.map((p) => (
                      <MenuItem key={`protocol:${p.id}`} value={encodeAxisValue('protocol', p.id)}>
                        {p.name}
                      </MenuItem>
                    ))}
                    {includedProperties.length > 0 && (
                      <ListSubheader>Molecular Properties</ListSubheader>
                    )}
                    {includedProperties.map((propName) => {
                      const meta = PROPERTY_LABELS[propName];
                      return (
                        <MenuItem
                          key={`property:${propName}`}
                          value={encodeAxisValue('property', propName)}
                        >
                          {meta?.label ?? propName}
                          {meta?.description ? ` — ${meta.description}` : ''}
                        </MenuItem>
                      );
                    })}
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
                <FormControl size="small" sx={{ minWidth: 180 }}>
                  <InputLabel>Axis</InputLabel>
                  <Select
                    value={yAxis.key ? encodeAxisValue(yAxis.kind, yAxis.key) : ''}
                    label="Axis"
                    onChange={(e) => handleAxisChange(setYAxis, e.target.value as string)}
                  >
                    {protocols.length > 0 && (
                      <ListSubheader>Protocols</ListSubheader>
                    )}
                    {protocols.map((p) => (
                      <MenuItem key={`protocol:${p.id}`} value={encodeAxisValue('protocol', p.id)}>
                        {p.name}
                      </MenuItem>
                    ))}
                    {includedProperties.length > 0 && (
                      <ListSubheader>Molecular Properties</ListSubheader>
                    )}
                    {includedProperties.map((propName) => {
                      const meta = PROPERTY_LABELS[propName];
                      return (
                        <MenuItem
                          key={`property:${propName}`}
                          value={encodeAxisValue('property', propName)}
                        >
                          {meta?.label ?? propName}
                          {meta?.description ? ` — ${meta.description}` : ''}
                        </MenuItem>
                      );
                    })}
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
            {xAxis.key && yAxis.key ? (
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
                    No compounds have values for both selected axes
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
                  Select an axis for both X and Y
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
                          <strong>{xAxisLabel}:</strong> {hoveredPoint.x.toFixed(2)}
                        </Typography>
                        <Typography variant="body2">
                          <strong>{yAxisLabel}:</strong> {hoveredPoint.y.toFixed(2)}
                        </Typography>
                      </Box>
                    </Box>
                  )}
                </Paper>
              </Fade>
            )}
          </Popper>
    </Box>
  );

  if (mode === 'inline') return panel;

  return (
    <>
      <Button
        variant="outlined"
        startIcon={<BubbleChart />}
        onClick={() => setOpen(true)}
        sx={{ mb: 2 }}
      >
        Scatter Plot
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
            <Typography variant="h6">Scatter Plot</Typography>
          </Box>
          <IconButton onClick={() => setOpen(false)} size="small">
            <Close />
          </IconButton>
        </DialogTitle>
        <DialogContent sx={{ display: 'flex', flexDirection: 'column', p: 0 }}>
          {panel}
        </DialogContent>
      </Dialog>
    </>
  );
});
