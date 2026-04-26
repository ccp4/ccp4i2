'use client';

import {
  AggregationResponse,
  AggregationType,
  CardContent,
  CompactRow,
  ConcentrationDisplayMode,
  isCompactResponse,
  isMediumResponse,
  OutputFormat,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import type { ScorecardConfig } from '@/types/compounds/models';
import { CompactTable } from './CompactTable';
import { MediumTable } from './MediumTable';
import { LongTable } from './LongTable';
import { PivotTable } from './PivotTable';
import { CardsView } from './CardsView';
import { BulletsView } from './BulletsView';

interface Props {
  /** Display data — incorporates threshold overrides for compact responses. */
  effectiveData: AggregationResponse;
  /** Original data — passed straight through to MediumTable / LongTable, since
   *  threshold overrides only patch compact responses' protocols[]. */
  data: AggregationResponse;
  aggregations: AggregationType[];
  outputFormat?: OutputFormat;
  displayMode: ConcentrationDisplayMode;
  searchTerm: string;
  fillHeight: boolean;
  scorecardConfig?: ScorecardConfig | null;
  cardContent?: CardContent;
  onCardContentChange?: (next: CardContent) => void;
  onEditProtocol?: (protocol: ProtocolInfo) => void;
  onScatterProtocol?: (protocol: ProtocolInfo) => void;
}

/** Dispatch into the right per-format sub-view component. */
export function AggregationView({
  effectiveData,
  data,
  aggregations,
  outputFormat,
  displayMode,
  searchTerm,
  fillHeight,
  scorecardConfig,
  cardContent,
  onCardContentChange,
  onEditProtocol,
  onScatterProtocol,
}: Props) {
  if (outputFormat === 'pivot' && isCompactResponse(effectiveData)) {
    return (
      <PivotTable
        data={effectiveData}
        aggregations={aggregations}
        concentrationDisplay={displayMode}
        onEditProtocol={onEditProtocol}
        onScatterProtocol={onScatterProtocol}
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
        onScatterProtocol={onScatterProtocol}
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
  if (isCompactResponse(effectiveData)) {
    return (
      <CompactTable
        data={effectiveData}
        aggregations={aggregations}
        concentrationDisplay={displayMode}
        searchTerm={searchTerm}
        fillHeight={fillHeight}
        onEditProtocol={onEditProtocol}
        onScatterProtocol={onScatterProtocol}
      />
    );
  }
  if (isMediumResponse(data)) {
    return (
      <MediumTable
        data={data}
        aggregations={aggregations}
        concentrationDisplay={displayMode}
        searchTerm={searchTerm}
        fillHeight={fillHeight}
      />
    );
  }
  return (
    <LongTable
      data={data}
      concentrationDisplay={displayMode}
      searchTerm={searchTerm}
      fillHeight={fillHeight}
    />
  );
}

// Re-export for callers that want the row type alongside the dispatcher.
export type { CompactRow };
