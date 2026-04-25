'use client';

import { useState, useRef, useCallback, useMemo } from 'react';
import {
  Box,
  Paper,
  Typography,
  Button,
  Chip,
  Tooltip,
  FormControl,
  Select,
  MenuItem,
  IconButton,
} from '@mui/material';
import { ZoomIn, ContentCopy, Check, EditOutlined, BubbleChart } from '@mui/icons-material';
import html2canvas from 'html2canvas';
import { useRouter } from 'next/navigation';
import {
  AggregationType,
  CardContent,
  CompactRow,
  ProtocolInfo,
  ConcentrationDisplayMode,
  MolecularPropertyName,
  MolecularPropertyThreshold,
  MolecularPropertyValues,
  getRagStatus,
  CompactAggregationResponse,
} from '@/types/compounds/aggregation';
import { MoleculeChip, CompoundNameChip } from '../MoleculeView';
import { DataSeriesDetailModal } from '../DataSeriesDetailModal';
import { CompoundSpider } from '../CompoundSpider';
import { protocolColour } from '@/lib/compounds/protocol-colour';
import type { ScorecardConfig } from '@/types/compounds/models';
import {
  Order,
  PROPERTY_LABELS,
  RAG_COLORS,
  IdentifiersCell,
  EMPTY_THRESHOLDS,
  formatPropertyValue,
  getComparator,
  getMeasurementStatus,
  formatMeasurementWithStats,
  TestedNoValidBadge,
  NotTestedIndicator,
} from './shared';

/** Sort key type for cards view */
type CardsSortKey = 'compound' | `property_${MolecularPropertyName}` | `protocol_${string}`;

/**
 * Cards view component (grid of compound cards).
 * Each card shows structure, compound info, and all protocols with formatted measurements.
 */
export function CardsView({
  data,
  aggregations,
  concentrationDisplay = 'natural',
  onEditProtocol,
  onScatterProtocol,
  scorecardConfig,
  cardContent: cardContentProp,
  onCardContentChange,
}: {
  data: CompactAggregationResponse;
  aggregations: AggregationType[];
  concentrationDisplay: ConcentrationDisplayMode;
  onEditProtocol?: (protocol: ProtocolInfo) => void;
  onScatterProtocol?: (protocol: ProtocolInfo) => void;
  scorecardConfig?: ScorecardConfig | null;
  /** Optional controlled value for the body selector — when omitted,
   *  CardsView falls back to internal local state. The aggregation page
   *  threads this through so the toggle is URL-encoded. */
  cardContent?: CardContent;
  onCardContentChange?: (next: CardContent) => void;
}) {
  const router = useRouter();

  // Modal state for data series detail
  const [modalOpen, setModalOpen] = useState(false);
  const [selectedCompound, setSelectedCompound] = useState<{
    id: string;
    name: string;
  } | null>(null);
  const [selectedProtocol, setSelectedProtocol] = useState<{
    id: string;
    name: string;
  } | null>(null);

  // Sorting state
  const [orderBy, setOrderBy] = useState<CardsSortKey>('compound');
  const [order, setOrder] = useState<Order>('asc');

  // Card body content selector — chemists may want spider-only for SAR
  // glance, protocols-only for raw numbers, or both. When the parent
  // controls this (URL-encoded on the aggregation page), we use the
  // controlled value; otherwise we fall back to local state defaulted
  // to whichever makes sense given the scorecard configuration.
  const hasScorecard = Boolean(scorecardConfig?.axes?.length);
  const [cardContentLocal, setCardContentLocal] = useState<CardContent>(
    hasScorecard ? 'both' : 'protocols',
  );
  const cardContent = cardContentProp ?? cardContentLocal;
  const setCardContent = (next: CardContent) => {
    if (onCardContentChange) onCardContentChange(next);
    else setCardContentLocal(next);
  };

  // Copy to clipboard state
  const [copiedCardId, setCopiedCardId] = useState<string | null>(null);
  const cardRefs = useRef<Map<string, HTMLDivElement>>(new Map());

  const rows = data.data as CompactRow[];
  const protocols = data.protocols;
  const includeProperties = data.meta.include_properties || [];
  const propertyThresholds = data.property_thresholds ?? EMPTY_THRESHOLDS;
  const showBatch = data.meta.group_by_batch;
  const showIdentifiers = Boolean(data.meta.include_identifiers);

  // Check which aggregations are selected
  const hasGeomean = aggregations.includes('geomean');
  const hasStdev = aggregations.includes('stdev');
  const hasCount = aggregations.includes('count');
  const hasList = aggregations.includes('list');

  // Create a map for quick threshold lookup
  const thresholdMap = useMemo(() => {
    const map: Record<string, MolecularPropertyThreshold> = {};
    for (const t of propertyThresholds) {
      map[t.property_name] = t;
    }
    return map;
  }, [propertyThresholds]);

  // Build sort options for dropdown
  const sortOptions = useMemo(() => {
    const options: { value: CardsSortKey; label: string }[] = [
      { value: 'compound', label: 'Compound ID' },
    ];
    // Add property options
    for (const propName of includeProperties) {
      options.push({
        value: `property_${propName}` as CardsSortKey,
        label: PROPERTY_LABELS[propName] || propName,
      });
    }
    // Add protocol options
    for (const protocol of protocols) {
      options.push({
        value: `protocol_${protocol.id}` as CardsSortKey,
        label: protocol.name,
      });
    }
    return options;
  }, [includeProperties, protocols]);

  // Get sort value for a compound
  const getSortValue = useCallback((row: CompactRow, key: CardsSortKey): unknown => {
    if (key === 'compound') return row.formatted_id;
    if (key.startsWith('property_')) {
      const propName = key.replace('property_', '') as MolecularPropertyName;
      return row.properties?.[propName] ?? null;
    }
    if (key.startsWith('protocol_')) {
      const protocolId = key.replace('protocol_', '');
      const protocolData = row.protocols[protocolId];
      if (hasGeomean) return protocolData?.geomean ?? null;
      if (hasCount) return protocolData?.count ?? null;
      return null;
    }
    return null;
  }, [hasGeomean, hasCount]);

  // Sorted cards
  const sortedRows = useMemo(() => {
    const comparator = getComparator<CompactRow>(order, (row) => getSortValue(row, orderBy));
    return [...rows].sort(comparator);
  }, [rows, order, orderBy, getSortValue]);

  // Handle click on protocol to show data series detail
  const handleProtocolClick = (
    event: React.MouseEvent,
    row: CompactRow,
    protocol: ProtocolInfo
  ) => {
    event.stopPropagation(); // Prevent navigation to compound
    setSelectedCompound({ id: row.compound_id, name: row.formatted_id });
    setSelectedProtocol({ id: protocol.id, name: protocol.name });
    setModalOpen(true);
  };

  const handleCloseModal = () => {
    setModalOpen(false);
    setSelectedCompound(null);
    setSelectedProtocol(null);
  };

  // Copy card to clipboard as image
  const handleCopyCard = useCallback(async (cardId: string, event: React.MouseEvent) => {
    event.stopPropagation(); // Prevent navigation to compound
    const cardElement = cardRefs.current.get(cardId);
    if (!cardElement) return;

    try {
      const canvas = await html2canvas(cardElement, {
        backgroundColor: '#ffffff',
        scale: 2, // Higher resolution for presentations
        logging: false,
      });

      canvas.toBlob(async (blob: Blob | null) => {
        if (!blob) return;
        try {
          await navigator.clipboard.write([
            new ClipboardItem({ 'image/png': blob }),
          ]);
          setCopiedCardId(cardId);
          setTimeout(() => setCopiedCardId(null), 2000);
        } catch (err) {
          console.error('Failed to copy to clipboard:', err);
        }
      }, 'image/png');
    } catch (err) {
      console.error('Failed to capture card:', err);
    }
  }, []);

  return (
    <>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2, flexWrap: 'wrap', gap: 1 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
          <Typography variant="body2" color="text.secondary">
            {data.meta.compound_count} compounds
            {data.meta.group_by_batch && data.meta.row_count && data.meta.row_count !== data.meta.compound_count && (
              <> ({data.meta.row_count} cards with batch split)</>
            )}
            , {data.meta.protocol_count} protocols,{' '}
            {data.meta.total_measurements} measurements
          </Typography>
          <Tooltip title="Click any protocol to view dose-response charts">
            <Chip
              icon={<ZoomIn fontSize="small" />}
              label="Click protocol for details"
              size="small"
              variant="outlined"
              color="info"
            />
          </Tooltip>
        </Box>
        {/* Sort controls + content selector */}
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          {hasScorecard && (
            <>
              <Typography variant="body2" color="text.secondary">
                Show:
              </Typography>
              <FormControl size="small" sx={{ minWidth: 110 }}>
                <Select
                  value={cardContent}
                  onChange={(e) => setCardContent(e.target.value as CardContent)}
                  size="small"
                  sx={{ fontSize: '0.875rem' }}
                >
                  <MenuItem value="both">Both</MenuItem>
                  <MenuItem value="spider">Spider only</MenuItem>
                  <MenuItem value="protocols">Protocols only</MenuItem>
                </Select>
              </FormControl>
            </>
          )}
          <Typography variant="body2" color="text.secondary">
            Sort by:
          </Typography>
          <FormControl size="small" sx={{ minWidth: 140 }}>
            <Select
              value={orderBy}
              onChange={(e) => setOrderBy(e.target.value as CardsSortKey)}
              size="small"
              sx={{ fontSize: '0.875rem' }}
            >
              {sortOptions.map((opt) => (
                <MenuItem key={opt.value} value={opt.value}>
                  {opt.label}
                </MenuItem>
              ))}
            </Select>
          </FormControl>
          <Tooltip title={order === 'asc' ? 'Ascending' : 'Descending'}>
            <Button
              size="small"
              variant="outlined"
              onClick={() => setOrder(order === 'asc' ? 'desc' : 'asc')}
              sx={{ minWidth: 36, px: 1 }}
            >
              {order === 'asc' ? '↑' : '↓'}
            </Button>
          </Tooltip>
        </Box>
      </Box>

      <Box
        sx={{
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fill, minmax(400px, 1fr))',
          gap: 2,
          maxHeight: 600,
          overflow: 'auto',
          p: 1,
        }}
      >
        {sortedRows.map((row) => {
          const cardId = showBatch ? `${row.compound_id}-${row.batch_id}` : row.compound_id;
          const isCopied = copiedCardId === cardId;
          return (
          <Paper
            key={cardId}
            ref={(el) => {
              if (el) cardRefs.current.set(cardId, el);
              else cardRefs.current.delete(cardId);
            }}
            variant="outlined"
            sx={{
              p: 2,
              cursor: 'pointer',
              position: 'relative',
              '&:hover': { boxShadow: 2, borderColor: 'primary.main' },
              '&:hover .copy-button': { opacity: 1 },
              display: 'flex',
              flexDirection: 'column',
            }}
            onClick={() => router.push(`/registry/compounds/${row.compound_id}`)}
          >
            {/* Copy to clipboard button */}
            <Tooltip title={isCopied ? 'Copied!' : 'Copy card as image'}>
              <IconButton
                className="copy-button"
                size="small"
                onClick={(e) => handleCopyCard(cardId, e)}
                sx={{
                  position: 'absolute',
                  top: 8,
                  right: 8,
                  opacity: 0,
                  transition: 'opacity 0.2s',
                  bgcolor: 'background.paper',
                  border: 1,
                  borderColor: 'divider',
                  '&:hover': { bgcolor: 'action.hover' },
                  zIndex: 1,
                }}
              >
                {isCopied ? <Check fontSize="small" color="success" /> : <ContentCopy fontSize="small" />}
              </IconButton>
            </Tooltip>
            {/* Header: Structure + Compound ID (+ optional per-compound spider) */}
            <Box sx={{ display: 'flex', gap: 2, mb: 1.5, alignItems: 'flex-start' }}>
              {row.smiles ? (
                <MoleculeChip smiles={row.smiles} size={180} />
              ) : (
                <Box
                  sx={{
                    width: 180,
                    height: 180,
                    bgcolor: 'grey.100',
                    borderRadius: 1,
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    flexShrink: 0,
                  }}
                >
                  <Typography variant="caption" color="text.secondary">-</Typography>
                </Box>
              )}
              <Box sx={{ flex: 1, minWidth: 0 }}>
                <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                  <CompoundNameChip formattedId={row.formatted_id} smiles={row.smiles} chipColor="primary" />
                  {showBatch && row.batch_number != null && (
                    <Typography variant="caption" color="text.secondary">
                      /{row.batch_number}
                    </Typography>
                  )}
                </Box>
                {row.target_name && (
                  <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5 }}>
                    Target: {row.target_name}
                  </Typography>
                )}
              </Box>
              {hasScorecard && cardContent !== 'protocols' && (
                <Box sx={{ flexShrink: 0 }}>
                  <CompoundSpider config={scorecardConfig} compound={row} size="small" />
                </Box>
              )}
            </Box>

            {/* Identifiers (barcode / supplier ref / aliases) */}
            {showIdentifiers && (
              <Box
                sx={{
                  py: 1,
                  borderTop: 1,
                  borderColor: 'divider',
                }}
              >
                <IdentifiersCell identifiers={row.identifiers} />
              </Box>
            )}

            {/* Molecular Properties (if any) */}
            {includeProperties.length > 0 && (
              <Box
                sx={{
                  display: 'flex',
                  flexWrap: 'wrap',
                  gap: 1,
                  py: 1,
                  borderTop: 1,
                  borderColor: 'divider',
                }}
              >
                {includeProperties.map((propName) => {
                  const value = row.properties?.[propName as keyof MolecularPropertyValues] as number | null | undefined;
                  const threshold = thresholdMap[propName];
                  const ragStatus = getRagStatus(value, threshold);
                  return (
                    <Tooltip key={propName} title={threshold?.property_display || propName}>
                      <Typography
                        variant="caption"
                        sx={{
                          bgcolor: 'grey.100',
                          px: 1,
                          py: 0.25,
                          borderRadius: 0.5,
                          color: RAG_COLORS[ragStatus],
                          fontFamily: 'monospace',
                        }}
                      >
                        {PROPERTY_LABELS[propName]}: {formatPropertyValue(value)}
                      </Typography>
                    </Tooltip>
                  );
                })}
              </Box>
            )}

            {/* Protocols — hidden when the user picked spider-only */}
            {cardContent !== 'spider' && (
            <Box sx={{ flex: 1, pt: 1, borderTop: 1, borderColor: 'divider' }}>
              {protocols.map((protocol) => {
                const protocolData = row.protocols[protocol.id];
                const measurementStatus = getMeasurementStatus(protocolData, hasCount);

                // Not measured - compound was never tested with this protocol
                if (measurementStatus === 'not-measured') {
                  return (
                    <Box key={protocol.id} sx={{ mb: 1 }}>
                      <Typography variant="body2" fontWeight={500}>
                        {protocol.name}
                      </Typography>
                      <NotTestedIndicator />
                    </Box>
                  );
                }

                // Tested but no valid data - still clickable to see what data exists
                if (measurementStatus === 'tested-no-valid') {
                  return (
                    <Box
                      key={protocol.id}
                      sx={{
                        mb: 1,
                        cursor: 'pointer',
                        p: 0.5,
                        mx: -0.5,
                        borderRadius: 0.5,
                        '&:hover': { bgcolor: 'action.hover' },
                      }}
                      onClick={(e) => handleProtocolClick(e, row, protocol)}
                    >
                      <Typography variant="body2" fontWeight={500}>
                        {protocol.name}
                      </Typography>
                      <TestedNoValidBadge
                        tested={protocolData?.tested}
                        noAnalysis={protocolData?.no_analysis}
                        invalid={protocolData?.invalid}
                        unassigned={protocolData?.unassigned}
                        onClick={(e) => handleProtocolClick(e, row, protocol)}
                      />
                    </Box>
                  );
                }

                // Format the measurement based on available aggregations
                let measurementDisplay: string;
                if (hasGeomean) {
                  measurementDisplay = formatMeasurementWithStats(
                    protocolData!.geomean,
                    protocolData!.stdev,
                    protocolData!.stdev_log,
                    protocolData!.count,
                    protocol.kpi_unit,
                    concentrationDisplay,
                    hasGeomean,
                    hasStdev,
                    hasCount
                  );
                } else if (hasCount) {
                  measurementDisplay = `n=${protocolData!.count ?? '-'}`;
                } else if (hasList) {
                  measurementDisplay = protocolData!.list || '-';
                } else {
                  measurementDisplay = '-';
                }

                const cardColour = hasGeomean
                  ? protocolColour(protocolData?.geomean, protocol)
                  : { background: null, t: null };

                return (
                  <Box
                    key={protocol.id}
                    sx={{
                      mb: 1,
                      cursor: 'pointer',
                      p: 0.5,
                      mx: -0.5,
                      borderRadius: 0.5,
                      bgcolor: cardColour.background ?? undefined,
                      '&:hover': cardColour.background
                        ? { filter: 'brightness(0.95)' }
                        : { bgcolor: 'action.hover' },
                    }}
                    onClick={(e) => handleProtocolClick(e, row, protocol)}
                  >
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.25 }}>
                      <Typography variant="body2" fontWeight={500}>
                        {protocol.name}
                      </Typography>
                      {onEditProtocol && (
                        <Tooltip title="Edit thresholds for this protocol" arrow>
                          <IconButton
                            size="small"
                            onClick={(e) => {
                              e.stopPropagation();
                              onEditProtocol(protocol);
                            }}
                            sx={{ p: 0.25, opacity: 0.5, '&:hover': { opacity: 1 } }}
                          >
                            <EditOutlined fontSize="inherit" />
                          </IconButton>
                        </Tooltip>
                      )}
                      {onScatterProtocol && (
                        <Tooltip title="Scatter this protocol vs another" arrow>
                          <IconButton
                            size="small"
                            onClick={(e) => {
                              e.stopPropagation();
                              onScatterProtocol(protocol);
                            }}
                            sx={{ p: 0.25, opacity: 0.5, '&:hover': { opacity: 1 } }}
                          >
                            <BubbleChart fontSize="inherit" />
                          </IconButton>
                        </Tooltip>
                      )}
                    </Box>
                    <Typography
                      variant="body2"
                      fontFamily="monospace"
                      sx={{
                        overflow: 'hidden',
                        textOverflow: 'ellipsis',
                        whiteSpace: 'nowrap',
                      }}
                    >
                      {measurementDisplay}
                    </Typography>
                  </Box>
                );
              })}
            </Box>
            )}
          </Paper>
        );
        })}
      </Box>

      {/* Data series detail modal */}
      {selectedCompound && selectedProtocol && (
        <DataSeriesDetailModal
          open={modalOpen}
          onClose={handleCloseModal}
          compoundId={selectedCompound.id}
          protocolId={selectedProtocol.id}
          compoundName={selectedCompound.name}
          protocolName={selectedProtocol.name}
        />
      )}
    </>
  );
}
