'use client';

import {
  Box,
  Typography,
  Chip,
  Tooltip,
  FormControl,
  Select,
  MenuItem,
  InputLabel,
} from '@mui/material';
import {
  AggregationType,
  ConcentrationDisplayMode,
  MolecularPropertyName,
  MolecularPropertyThreshold,
  CompoundIdentifiers,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import {
  formatConcentrationValue,
} from '@/lib/compounds/aggregation-api';
import { protocolColour } from '@/lib/compounds/protocol-colour';

/**
 * Module-level empty array for compact/pivot/cards views to fall back to
 * when a response has no `property_thresholds`. Inline `[]` literals
 * allocate a fresh array per render, invalidating downstream `useMemo`s
 * that key on the threshold list.
 */
export const EMPTY_THRESHOLDS: MolecularPropertyThreshold[] = [];

/** Options for concentration display mode selector */
export const CONCENTRATION_DISPLAY_OPTIONS: { value: ConcentrationDisplayMode; label: string }[] = [
  { value: 'natural', label: 'Natural' },
  { value: 'nM', label: 'nM' },
  { value: 'uM', label: 'μM' },
  { value: 'mM', label: 'mM' },
  { value: 'pConc', label: 'pConc' },
];

/** Short labels for molecular properties */
export const PROPERTY_LABELS: Record<MolecularPropertyName, string> = {
  molecular_weight: 'MW',
  heavy_atom_count: 'HA',
  hbd: 'HBD',
  hba: 'HBA',
  clogp: 'cLogP',
  tpsa: 'TPSA',
  rotatable_bonds: 'RotB',
  fraction_sp3: 'Fsp3',
};

/** RAG status colors */
export const RAG_COLORS: Record<string, string> = {
  green: 'inherit',
  amber: '#ff9800',
  red: '#f44336',
};

/**
 * Compact display of a compound's barcode, supplier reference, and aliases.
 * Renders "-" when the identifiers record is present but has no content.
 */
export function IdentifiersCell({ identifiers }: { identifiers?: CompoundIdentifiers | null }) {
  if (!identifiers) {
    return <Typography variant="caption" color="text.secondary">-</Typography>;
  }
  const { barcode, supplier_ref, aliases } = identifiers;
  const hasAny = Boolean(barcode) || Boolean(supplier_ref) || (aliases && aliases.length > 0);
  if (!hasAny) {
    return <Typography variant="caption" color="text.secondary">-</Typography>;
  }
  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', gap: 0.25, minWidth: 0 }}>
      {barcode && (
        <Typography variant="caption" fontFamily="monospace" sx={{ lineHeight: 1.2, whiteSpace: 'nowrap' }}>
          {barcode}
        </Typography>
      )}
      {supplier_ref && (
        <Typography
          variant="caption"
          color="text.secondary"
          sx={{ lineHeight: 1.2, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}
          title={supplier_ref}
        >
          {supplier_ref}
        </Typography>
      )}
      {aliases && aliases.length > 0 && (
        <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.25, mt: 0.25 }}>
          {aliases.map((alias, i) => (
            <Chip
              key={`${alias}-${i}`}
              label={alias}
              size="small"
              variant="outlined"
              sx={{ height: 16, fontSize: '0.65rem', '& .MuiChip-label': { px: 0.5 } }}
            />
          ))}
        </Box>
      )}
    </Box>
  );
}

/** Format property value for display */
export function formatPropertyValue(value: number | string | null | undefined): string {
  if (value == null || value === '') return '-';
  // Convert string values to numbers (handles legacy data stored as strings)
  const numValue = typeof value === 'string' ? parseFloat(value) : value;
  if (isNaN(numValue)) return String(value);
  if (Number.isInteger(numValue)) return String(numValue);
  return numValue.toFixed(2);
}

/**
 * Format aggregation label for column header.
 */
export function formatAggLabel(agg: AggregationType): string {
  const labels: Record<AggregationType, string> = {
    geomean: 'GeoMean',
    count: 'N',
    stdev: 'StDev',
    list: 'Values',
  };
  return labels[agg] || agg;
}

/** Sort direction type */
export type Order = 'asc' | 'desc';

/** Comparator function for sorting - compares non-null values only */
function compareValues(aVal: unknown, bVal: unknown, order: Order): number {
  // Compare values based on type
  if (typeof aVal === 'string' && typeof bVal === 'string') {
    const result = aVal.localeCompare(bVal);
    return order === 'asc' ? result : -result;
  }
  if (typeof aVal === 'number' && typeof bVal === 'number') {
    const result = aVal - bVal;
    return order === 'asc' ? result : -result;
  }
  return 0;
}

export function getComparator<T>(
  order: Order,
  getValue: (item: T) => unknown
): (a: T, b: T) => number {
  return (a, b) => {
    const aVal = getValue(a);
    const bVal = getValue(b);

    // Handle null/undefined - ALWAYS push to end regardless of sort direction
    if (aVal == null && bVal == null) return 0;
    if (aVal == null) return 1;  // a (null) goes after b
    if (bVal == null) return -1; // b (null) goes after a

    // Compare non-null values with order direction
    return compareValues(aVal, bVal, order);
  };
}

/**
 * Represents measurement status for a compound-protocol pair.
 */
type MeasurementStatus = 'not-measured' | 'tested-no-valid' | 'has-data';

/**
 * Determine the measurement status for a protocol cell.
 * Distinguishes between not measured, tested but no valid data, and has data.
 *
 * Logic:
 * - tested === 0 (or undefined/null): compound was never tested → 'not-measured'
 * - tested > 0 && count === 0: compound was tested but no valid results → 'tested-no-valid'
 * - count > 0: has valid data → 'has-data'
 */
export function getMeasurementStatus(
  protocolData: { count?: number; tested?: number } | undefined,
  hasCount: boolean
): MeasurementStatus {
  // If no protocol data at all, compound was not measured
  if (!protocolData || Object.keys(protocolData).length === 0) {
    return 'not-measured';
  }
  // Check tested count first - if 0 or undefined, compound was never tested
  const tested = protocolData.tested ?? 0;
  if (tested === 0) {
    return 'not-measured';
  }
  // If tested > 0 but count is 0, compound was tested but no valid results
  if (hasCount && protocolData.count === 0) {
    return 'tested-no-valid';
  }
  // Has actual data
  return 'has-data';
}

/**
 * Format measurement with statistics for cards view.
 * Shows value ± stdev (n=count) when all data is available.
 *
 * In pConc mode, uses stdev_log (standard deviation calculated in log-space)
 * for meaningful error display.
 */
export function formatMeasurementWithStats(
  geomean: number | null | undefined,
  stdev: number | null | undefined,
  stdevLog: number | null | undefined,
  count: number | undefined,
  unit: string | null | undefined,
  concentrationDisplay: ConcentrationDisplayMode,
  hasGeomean: boolean,
  hasStdev: boolean,
  hasCount: boolean
): string {
  if (geomean == null) return '-';

  const formatted = formatConcentrationValue(geomean, unit, concentrationDisplay);
  let result = formatted.displayValue;

  // Add stdev if available and requested
  if (hasStdev) {
    if (concentrationDisplay === 'pConc') {
      // In pConc mode, use log-space stdev directly (it's already in log units)
      if (stdevLog != null) {
        result += ` ± ${stdevLog.toFixed(2)}`;
      }
    } else if (stdev != null) {
      // In linear modes, use regular stdev with unit conversion
      const stdevFormatted = formatConcentrationValue(stdev, unit, concentrationDisplay);
      result += ` ± ${stdevFormatted.displayValue}`;
    }
  }

  // Add unit (only in natural mode)
  if (formatted.displayUnit && concentrationDisplay === 'natural') {
    result += ` ${formatted.displayUnit}`;
  }

  // Add count if available and requested
  if (hasCount && count != null) {
    result += ` (n=${count})`;
  }

  return result;
}

/**
 * Display component for "tested but no valid data" status.
 * Shows a subtle indicator that compound was tested but all results were invalid.
 * Optionally shows breakdown of tested vs no_analysis counts.
 * Click to view the underlying data series.
 */
export function TestedNoValidBadge({
  tested,
  noAnalysis,
  invalid,
  unassigned,
  onClick,
}: {
  tested?: number;
  noAnalysis?: number;
  invalid?: number;
  unassigned?: number;
  onClick?: (e: React.MouseEvent) => void;
}) {
  // Build informative tooltip showing breakdown of non-valid results
  const parts: string[] = [];
  if (noAnalysis != null && noAnalysis > 0) {
    parts.push(`${noAnalysis} no analysis`);
  }
  if (invalid != null && invalid > 0) {
    parts.push(`${invalid} invalid`);
  }
  if (unassigned != null && unassigned > 0) {
    parts.push(`${unassigned} unassigned`);
  }

  let tooltipText = 'Tested but no valid results';
  if (tested != null && tested > 0) {
    if (parts.length > 0) {
      tooltipText = `${tested} test${tested > 1 ? 's' : ''}: ${parts.join(', ')}`;
    } else {
      tooltipText = `${tested} test${tested > 1 ? 's' : ''}, 0 valid results`;
    }
  }
  if (onClick) {
    tooltipText += ' — click to view';
  }

  // Determine label and color based on primary issue
  const hasNoAnalysis = noAnalysis != null && noAnalysis > 0;
  const hasInvalid = invalid != null && invalid > 0;
  const totalNonValid = (noAnalysis || 0) + (invalid || 0) + (unassigned || 0);

  // Pick a representative label
  let label: string;
  if (hasNoAnalysis && noAnalysis === tested) {
    label = 'no analysis';
  } else if (hasInvalid && invalid === tested) {
    label = 'invalid';
  } else if (totalNonValid > 0) {
    label = `0/${tested}`;
  } else {
    label = '0 valid';
  }

  return (
    <Tooltip title={tooltipText}>
      <Chip
        label={label}
        size="small"
        variant="outlined"
        onClick={onClick}
        sx={{
          height: 20,
          fontSize: '0.7rem',
          color: hasNoAnalysis ? 'warning.main' : hasInvalid ? 'error.main' : 'text.secondary',
          borderColor: hasNoAnalysis ? 'warning.light' : hasInvalid ? 'error.light' : 'grey.400',
          bgcolor: hasNoAnalysis ? 'warning.50' : hasInvalid ? 'error.50' : 'grey.50',
          cursor: onClick ? 'pointer' : 'default',
          '&:hover': onClick ? {
            borderColor: hasNoAnalysis ? 'warning.main' : hasInvalid ? 'error.main' : 'grey.600',
            bgcolor: hasNoAnalysis ? 'warning.100' : hasInvalid ? 'error.100' : 'grey.100',
          } : {},
        }}
      />
    </Tooltip>
  );
}

/**
 * Display component for "not tested" status.
 * Shows a subtle grayed indicator that compound was never tested for this protocol.
 * Visually distinct from TestedNoValidBadge (tested but invalid).
 */
export function NotTestedIndicator() {
  return (
    <Tooltip title="Not tested for this protocol">
      <Typography
        variant="body2"
        component="span"
        sx={{
          color: 'grey.400',
          fontFamily: 'monospace',
          fontStyle: 'italic',
          fontSize: '0.75rem',
        }}
      >
        n/t
      </Typography>
    </Tooltip>
  );
}

/**
 * Concentration display mode selector component.
 */
export function ConcentrationDisplaySelector({
  value,
  onChange,
}: {
  value: ConcentrationDisplayMode;
  onChange: (mode: ConcentrationDisplayMode) => void;
}) {
  return (
    <FormControl size="small" sx={{ minWidth: 110 }}>
      <InputLabel>Concentration</InputLabel>
      <Select
        value={value}
        label="Concentration"
        onChange={(e) => onChange(e.target.value as ConcentrationDisplayMode)}
      >
        {CONCENTRATION_DISPLAY_OPTIONS.map((opt) => (
          <MenuItem key={opt.value} value={opt.value}>
            {opt.label}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  );
}

/**
 * Inline legend for the protocol-threshold colour scale.
 *
 * Renders a compact gradient sample (poor → excellent) when at least one
 * protocol in the current response has thresholds configured. The component
 * is self-gating: if no protocol has thresholds, it renders nothing, so
 * chemists only see the legend once it's actually explaining a colour in
 * the table.
 */
export function ThresholdLegend({ protocols }: { protocols: ProtocolInfo[] }) {
  const configured = protocols.filter(
    (p) => p.target_value != null && p.poor_value != null,
  );
  if (configured.length === 0) return null;

  const gradient = `linear-gradient(to right, ${protocolColour(0, { target_value: 1, poor_value: 0, threshold_scale: 'linear' }).background}, ${protocolColour(0.5, { target_value: 1, poor_value: 0, threshold_scale: 'linear' }).background}, ${protocolColour(1, { target_value: 1, poor_value: 0, threshold_scale: 'linear' }).background})`;

  const tooltipBody = (
    <Box sx={{ p: 0.5 }}>
      <Typography variant="caption" sx={{ display: 'block', mb: 0.5 }}>
        Cells are tinted against each protocol&apos;s curated thresholds.
      </Typography>
      <Typography variant="caption" sx={{ display: 'block', color: 'text.secondary' }}>
        {configured.length} of {protocols.length} protocol{protocols.length !== 1 ? 's' : ''} configured.
        Uncoloured cells mean the protocol has no thresholds set.
      </Typography>
    </Box>
  );

  return (
    <Tooltip title={tooltipBody} arrow>
      <Box
        sx={{
          display: 'flex',
          alignItems: 'center',
          gap: 0.75,
          px: 1,
          py: 0.25,
          border: 1,
          borderColor: 'divider',
          borderRadius: 1,
          bgcolor: 'background.paper',
        }}
      >
        <Typography variant="caption" color="text.secondary">
          poor
        </Typography>
        <Box
          sx={{
            width: 80,
            height: 12,
            borderRadius: 0.5,
            background: gradient,
            border: 1,
            borderColor: 'divider',
          }}
        />
        <Typography variant="caption" color="text.secondary">
          excellent
        </Typography>
      </Box>
    </Tooltip>
  );
}
