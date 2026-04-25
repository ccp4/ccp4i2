'use client';

import { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import {
  Box,
  Paper,
  Typography,
  LinearProgress,
  IconButton,
  TextField,
  InputAdornment,
} from '@mui/material';
import { Search, Clear } from '@mui/icons-material';
import {
  AggregationResponse,
  AggregationType,
  CardContent,
  CompactRow,
  isCompactResponse,
  isMediumResponse,
  ProtocolInfo,
  ConcentrationDisplayMode,
  OutputFormat,
} from '@/types/compounds/aggregation';
import {
  ProtocolScatterPlot,
  type ProtocolScatterPlotHandle,
} from './ProtocolScatterPlot';
import {
  ProtocolThresholdQuickEdit,
  type ProtocolThresholdUpdate,
} from './ProtocolThresholdQuickEdit';
import { useAuth } from '@/lib/compounds/auth-context';
import type { ScorecardConfig } from '@/types/compounds/models';
import {
  ConcentrationDisplaySelector,
  ThresholdLegend,
} from './aggregation/shared';
import { CompactTable } from './aggregation/CompactTable';
import { MediumTable } from './aggregation/MediumTable';
import { LongTable } from './aggregation/LongTable';
import { PivotTable } from './aggregation/PivotTable';
import { CardsView } from './aggregation/CardsView';
import { BulletsView } from './aggregation/BulletsView';

// Module-level empty array for stable identity when the scatter has no
// protocol axes to offer. Allocating `[]` inline on each render defeats
// ProtocolScatterPlot's prop-identity checks.
const NO_PROTOCOLS: ProtocolInfo[] = [];

interface AggregationTableProps {
  data: AggregationResponse | null | undefined;
  loading?: boolean;
  aggregations: AggregationType[];
  /** Output format to determine which view to render (default: auto-detect from data) */
  outputFormat?: OutputFormat;
  /** Concentration display mode (default: 'natural') */
  concentrationDisplay?: ConcentrationDisplayMode;
  /** Callback when concentration display mode changes */
  onConcentrationDisplayChange?: (mode: ConcentrationDisplayMode) => void;
  /** Target scorecard config, used by Cards view to render a per-compound spider.
   *  Passed down from pages that filter by a single target (so the config is unambiguous). */
  scorecardConfig?: ScorecardConfig | null;
  /** Cards-view body selector. URL-encoded on the aggregation page so deep
   *  links preserve the chosen layout; CardsView falls back to local state
   *  when this is omitted. */
  cardContent?: CardContent;
  onCardContentChange?: (next: CardContent) => void;
  /** If true, table fills available parent height instead of using fixed maxHeight */
  fillHeight?: boolean;
}

/**
 * Main aggregation table component.
 * Renders either compact or long format based on the response data.
 */
export function AggregationTable({
  data,
  loading,
  aggregations,
  outputFormat,
  concentrationDisplay = 'natural',
  onConcentrationDisplayChange,
  scorecardConfig,
  cardContent,
  onCardContentChange,
  fillHeight = false,
}: AggregationTableProps) {
  // Search state for filtering compounds in the table
  const [searchTerm, setSearchTerm] = useState('');

  // Clear search when data changes (new query)
  useEffect(() => {
    setSearchTerm('');
  }, [data]);

  // Threshold quick-edit state. `overrides` merges into data.protocols
  // so colours update live after save without needing a fresh query.
  // Cleared when data changes (next fetch reflects server state).
  const { canContribute } = useAuth();
  const [overrides, setOverrides] = useState<Record<string, ProtocolThresholdUpdate>>({});
  const [editingProtocol, setEditingProtocol] = useState<ProtocolInfo | null>(null);
  useEffect(() => {
    setOverrides({});
  }, [data]);

  const handleProtocolSaved = useCallback(
    (id: string, update: ProtocolThresholdUpdate) => {
      setOverrides((prev) => ({ ...prev, [id]: update }));
    },
    [],
  );

  // Ref to the scatter plot so sub-views can trigger a pre-filled open.
  const scatterRef = useRef<ProtocolScatterPlotHandle>(null);
  const handleScatterProtocol = useCallback((protocol: ProtocolInfo) => {
    scatterRef.current?.openForProtocol(protocol.id);
  }, []);

  // Merge overrides into data.protocols when rendering compact-style responses.
  // Non-compact responses (MediumRow/LongRow) don't carry protocols[] so pass data through.
  const effectiveData = useMemo(() => {
    if (!data || !isCompactResponse(data)) return data;
    if (Object.keys(overrides).length === 0) return data;
    return {
      ...data,
      protocols: data.protocols.map((p) =>
        overrides[p.id] ? { ...p, ...overrides[p.id] } : p,
      ),
    };
  }, [data, overrides]);

  // Use internal state if no external control is provided
  const [internalDisplay, setInternalDisplay] = useState<ConcentrationDisplayMode>(concentrationDisplay);

  // Use controlled or uncontrolled mode
  const displayMode = onConcentrationDisplayChange ? concentrationDisplay : internalDisplay;
  const handleDisplayChange = (mode: ConcentrationDisplayMode) => {
    if (onConcentrationDisplayChange) {
      onConcentrationDisplayChange(mode);
    } else {
      setInternalDisplay(mode);
    }
  };

  // NB: this callback MUST stay above the early-return guards below. Declaring
  // it after them makes the component's hook count vary with `loading`/`data`
  // state (React error #310 — "Rendered more hooks than during the previous
  // render") and white-screens every table view.
  const handleEditProtocol = useCallback(
    (protocol: ProtocolInfo) => setEditingProtocol(protocol),
    [],
  );
  const onEditProtocol = canContribute ? handleEditProtocol : undefined;

  if (loading) {
    return (
      <Paper sx={{ p: 3, ...(fillHeight && { height: '100%' }) }}>
        <LinearProgress />
        <Typography sx={{ mt: 2 }} color="text.secondary" align="center">
          Running aggregation query...
        </Typography>
      </Paper>
    );
  }

  if (!data) {
    return (
      <Paper sx={{ p: 3, ...(fillHeight && { height: '100%' }) }}>
        <Typography color="text.secondary" align="center">
          Select a target, protocol, or compound above to see results.
        </Typography>
      </Paper>
    );
  }

  if (data.meta.total_measurements === 0) {
    return (
      <Paper sx={{ p: 3, ...(fillHeight && { height: '100%' }) }}>
        <Typography color="text.secondary" align="center">
          No data found matching your criteria.
        </Typography>
      </Paper>
    );
  }

  // Determine which table to render based on response type and outputFormat
  const renderTable = () => {
    // Narrow effectiveData for the type-checker — the outer function has
    // already early-returned when data is null, so this can't actually fire.
    if (!effectiveData) return null;
    // For pivot and cards formats, we need compact response data
    if (outputFormat === 'pivot' && isCompactResponse(effectiveData)) {
      return (
        <PivotTable
          data={effectiveData}
          aggregations={aggregations}
          concentrationDisplay={displayMode}
          onEditProtocol={onEditProtocol}
          onScatterProtocol={handleScatterProtocol}
        />
      );
    }
    if (outputFormat === 'cards' && isCompactResponse(effectiveData)) {
      return (
        <CardsView
          data={effectiveData}
          aggregations={aggregations}
          concentrationDisplay={displayMode}
          onEditProtocol={onEditProtocol}
          onScatterProtocol={handleScatterProtocol}
          scorecardConfig={scorecardConfig}
          cardContent={cardContent}
          onCardContentChange={onCardContentChange}
        />
      );
    }
    if (outputFormat === 'bullets' && isCompactResponse(effectiveData)) {
      return (
        <BulletsView
          data={effectiveData}
          scorecardConfig={scorecardConfig}
          searchTerm={searchTerm}
          concentrationDisplay={displayMode}
          fillHeight={fillHeight}
        />
      );
    }
    // Default: auto-detect from response type
    if (isCompactResponse(effectiveData)) {
      return (
        <CompactTable
          data={effectiveData}
          aggregations={aggregations}
          concentrationDisplay={displayMode}
          searchTerm={searchTerm}
          fillHeight={fillHeight}
          onEditProtocol={onEditProtocol}
          onScatterProtocol={handleScatterProtocol}
        />
      );
    } else if (isMediumResponse(data)) {
      return (
        <MediumTable
          data={data}
          aggregations={aggregations}
          concentrationDisplay={displayMode}
          searchTerm={searchTerm}
          fillHeight={fillHeight}
        />
      );
    } else {
      return (
        <LongTable
          data={data}
          concentrationDisplay={displayMode}
          searchTerm={searchTerm}
          fillHeight={fillHeight}
        />
      );
    }
  };

  return (
    <Paper sx={{ width: '100%', overflow: 'clip', ...(fillHeight && { height: '100%', display: 'flex', flexDirection: 'column' }) }}>
      <Box sx={{ p: 2, ...(fillHeight && { flex: 1, display: 'flex', flexDirection: 'column', minHeight: 0, overflow: 'clip' }) }}>
        <Box sx={{ display: 'flex', justifyContent: 'flex-end', alignItems: 'center', mb: 1, flexShrink: 0, gap: 1 }}>
          {effectiveData && isCompactResponse(effectiveData) && (
            <ThresholdLegend protocols={effectiveData.protocols} />
          )}
          <TextField
            size="small"
            placeholder="Search compounds..."
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
            slotProps={{
              input: {
                startAdornment: (
                  <InputAdornment position="start">
                    <Search fontSize="small" color="action" />
                  </InputAdornment>
                ),
                endAdornment: searchTerm ? (
                  <InputAdornment position="end">
                    <IconButton size="small" onClick={() => setSearchTerm('')} edge="end">
                      <Clear fontSize="small" />
                    </IconButton>
                  </InputAdornment>
                ) : null,
              },
            }}
            sx={{ width: 220 }}
          />
          <ConcentrationDisplaySelector
            value={displayMode}
            onChange={handleDisplayChange}
          />
        </Box>
        <Box sx={{ ...(fillHeight && { flex: 1, minHeight: 0, display: 'flex', flexDirection: 'column' }) }}>
          {renderTable()}
        </Box>
      </Box>

      {/* Scatter plot is mounted at the top level so it's reachable from
          Pivot, Cards, and Compact views via the scatterRef imperative
          handle. Rendered for compact-style responses that have at least
          two axis candidates (protocols with geomean or molecular properties). */}
      {effectiveData && isCompactResponse(effectiveData) && (() => {
        const hasGeomean = aggregations.includes('geomean');
        const protocolAxes = hasGeomean ? effectiveData.protocols.length : 0;
        const propertyAxes = effectiveData.meta.include_properties?.length || 0;
        if (protocolAxes + propertyAxes < 2) return null;
        return (
          <ProtocolScatterPlot
            ref={scatterRef}
            data={effectiveData.data as CompactRow[]}
            protocols={hasGeomean ? effectiveData.protocols : NO_PROTOCOLS}
            includedProperties={effectiveData.meta.include_properties || []}
          />
        );
      })()}

      {/* Threshold quick-edit dialog — opened from per-protocol pencil icons. */}
      <ProtocolThresholdQuickEdit
        open={editingProtocol !== null}
        protocol={editingProtocol}
        onClose={() => setEditingProtocol(null)}
        onSaved={handleProtocolSaved}
      />
    </Paper>
  );
}
