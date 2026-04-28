'use client';

import { useState, useRef, useEffect, useCallback } from 'react';
import { Box, Paper } from '@mui/material';
import {
  AggregationResponse,
  AggregationType,
  CardContent,
  ConcentrationDisplayMode,
  OutputFormat,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import {
  ProtocolThresholdQuickEdit,
} from './ProtocolThresholdQuickEdit';
import type { ProtocolScatterPlotHandle, CategorisationSelection } from './ProtocolScatterPlot';
import type { ScorecardConfig } from '@/types/compounds/models';
import { AggregationStatus } from './aggregation/AggregationStatus';
import { AggregationToolbar } from './aggregation/AggregationToolbar';
import { AggregationView } from './aggregation/AggregationView';
import { AggregationScatter } from './aggregation/AggregationScatter';
import { useThresholdOverrides } from './aggregation/useThresholdOverrides';

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
  /** Selections used to colour scatter points by membership when format=scatter. */
  categorisationSelections?: ReadonlyArray<CategorisationSelection>;
}

/**
 * Top-level wrapper for the aggregation table. Owns search and
 * concentration-display state; threshold-override state lives in
 * useThresholdOverrides; the dispatch into per-format views is in
 * AggregationView; the imperatively-driven scatter plot is in
 * AggregationScatter.
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
  categorisationSelections,
}: AggregationTableProps) {
  const [searchTerm, setSearchTerm] = useState('');
  useEffect(() => {
    setSearchTerm('');
  }, [data]);

  const {
    effectiveData,
    editingProtocol,
    setEditingProtocol,
    handleProtocolSaved,
    onEditProtocol,
  } = useThresholdOverrides(data);

  const scatterRef = useRef<ProtocolScatterPlotHandle>(null);
  const handleScatterProtocol = useCallback((protocol: ProtocolInfo) => {
    scatterRef.current?.openForProtocol(protocol.id);
  }, []);

  const [internalDisplay, setInternalDisplay] = useState<ConcentrationDisplayMode>(concentrationDisplay);
  const displayMode = onConcentrationDisplayChange ? concentrationDisplay : internalDisplay;
  const handleDisplayChange = (mode: ConcentrationDisplayMode) => {
    if (onConcentrationDisplayChange) onConcentrationDisplayChange(mode);
    else setInternalDisplay(mode);
  };

  if (loading) return <AggregationStatus state="loading" fillHeight={fillHeight} />;
  if (!data) return <AggregationStatus state="no-selection" fillHeight={fillHeight} />;
  if (data.meta.total_measurements === 0) {
    return <AggregationStatus state="no-data" fillHeight={fillHeight} />;
  }

  return (
    <Paper sx={{ width: '100%', overflow: 'clip', ...(fillHeight && { height: '100%', display: 'flex', flexDirection: 'column' }) }}>
      <Box sx={{ p: 2, ...(fillHeight && { flex: 1, display: 'flex', flexDirection: 'column', minHeight: 0, overflow: 'clip' }) }}>
        <AggregationToolbar
          effectiveData={effectiveData}
          searchTerm={searchTerm}
          onSearchTermChange={setSearchTerm}
          displayMode={displayMode}
          onDisplayChange={handleDisplayChange}
        />
        <Box sx={{ ...(fillHeight && { flex: 1, minHeight: 0, display: 'flex', flexDirection: 'column' }) }}>
          {effectiveData && (
            <AggregationView
              effectiveData={effectiveData}
              data={data}
              aggregations={aggregations}
              outputFormat={outputFormat}
              displayMode={displayMode}
              searchTerm={searchTerm}
              fillHeight={fillHeight}
              scorecardConfig={scorecardConfig}
              cardContent={cardContent}
              onCardContentChange={onCardContentChange}
              onEditProtocol={onEditProtocol}
              onScatterProtocol={handleScatterProtocol}
              categorisationSelections={categorisationSelections}
            />
          )}
        </Box>
      </Box>

      <AggregationScatter
        ref={scatterRef}
        effectiveData={effectiveData}
        aggregations={aggregations}
      />

      <ProtocolThresholdQuickEdit
        open={editingProtocol !== null}
        protocol={editingProtocol}
        onClose={() => setEditingProtocol(null)}
        onSaved={handleProtocolSaved}
      />
    </Paper>
  );
}
