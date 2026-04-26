'use client';

import { forwardRef } from 'react';
import {
  AggregationResponse,
  AggregationType,
  CompactRow,
  isCompactResponse,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import {
  ProtocolScatterPlot,
  type ProtocolScatterPlotHandle,
} from '../ProtocolScatterPlot';

// Stable empty-array reference — inline `[]` re-allocates on every render
// and defeats ProtocolScatterPlot's prop-identity checks.
const NO_PROTOCOLS: ProtocolInfo[] = [];

interface Props {
  effectiveData: AggregationResponse | null | undefined;
  aggregations: AggregationType[];
}

/** Scatter plot rendered for compact-style responses with at least two
 *  axis candidates (geomean protocols + molecular property columns).
 *  Returns null when there's nothing useful to plot. The forwarded ref
 *  exposes the imperative `openForProtocol(id)` handle so sub-views can
 *  pre-fill the X axis from a per-row pencil click. */
export const AggregationScatter = forwardRef<ProtocolScatterPlotHandle, Props>(
  function AggregationScatter({ effectiveData, aggregations }, ref) {
    if (!effectiveData || !isCompactResponse(effectiveData)) return null;
    const hasGeomean = aggregations.includes('geomean');
    const protocolAxes = hasGeomean ? effectiveData.protocols.length : 0;
    const propertyAxes = effectiveData.meta.include_properties?.length || 0;
    if (protocolAxes + propertyAxes < 2) return null;
    return (
      <ProtocolScatterPlot
        ref={ref}
        data={effectiveData.data as CompactRow[]}
        protocols={hasGeomean ? effectiveData.protocols : NO_PROTOCOLS}
        includedProperties={effectiveData.meta.include_properties || []}
      />
    );
  },
);
