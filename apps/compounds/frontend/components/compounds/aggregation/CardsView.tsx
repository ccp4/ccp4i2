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
import { ZoomIn, ContentCopy, Check, EditOutlined, BubbleChart, Medication, Slideshow } from '@mui/icons-material';
import html2canvas from 'html2canvas';
import { exportCardsToPptx } from '@/lib/compounds/pptx-export';
import { pinCaptureFonts } from '@/lib/compounds/html2canvas-fonts';
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
import { MoleculeChip } from '../MoleculeView';
import { DataSeriesDetailModal } from '../DataSeriesDetailModal';
import { CompoundSpider } from '../CompoundSpider';
import { ScorecardValuesTable } from '../ScorecardValuesTable';
import { protocolColour } from '@/lib/compounds/protocol-colour';
import {
  evaluateScorecard,
  groupAxesBySector,
  sectorColour,
  type AxisEvaluation,
} from '@/lib/compounds/scorecard';
import type { ScorecardConfig } from '@/types/compounds/models';
import {
  formatAxisValueForBullet,
  tierLabel,
} from './BulletsView';
import {
  Order,
  MONOSPACE_FONT_STACK,
  SANS_SERIF_FONT_STACK,
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
type CardsSortKey =
  | 'compound'
  | `property_${MolecularPropertyName}`
  | `protocol_${string}`
  | `axis_${number}`;

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

  // PowerPoint export state
  const [exportCardsPerSlide, setExportCardsPerSlide] = useState<number>(6);
  const [exportProgress, setExportProgress] = useState<{ done: number; total: number } | null>(null);
  const [exportError, setExportError] = useState<string | null>(null);

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
  const scorecardAxes = scorecardConfig?.axes ?? [];
  const sortOptions = useMemo(() => {
    const options: { value: CardsSortKey; label: string }[] = [
      { value: 'compound', label: 'Compound ID' },
    ];
    // Scorecard axes first (most relevant when a scorecard is in play)
    for (let i = 0; i < scorecardAxes.length; i++) {
      const axis = scorecardAxes[i];
      options.push({
        value: `axis_${i}` as CardsSortKey,
        label: `★ ${axis.label || `axis ${i + 1}`}`,
      });
    }
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
  }, [includeProperties, protocols, scorecardAxes]);

  // Get sort value for a compound. For scorecard axes we sort by the
  // normalised score `t` (0–1, direction-aware) so 'best first' makes
  // sense regardless of whether the underlying axis is lower-better
  // (IC50) or higher-better (solubility).
  const getSortValue = useCallback((row: CompactRow, key: CardsSortKey): unknown => {
    if (key === 'compound') return row.formatted_id;
    if (key.startsWith('axis_')) {
      const idx = Number(key.slice(5));
      const axis = scorecardAxes[idx];
      if (!axis || !scorecardConfig) return null;
      const evals = evaluateScorecard(scorecardConfig, row);
      return evals[idx]?.t ?? null;
    }
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
  }, [hasGeomean, hasCount, scorecardAxes, scorecardConfig]);

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

  // Copy card to clipboard as image. Safari and recent Chrome both
  // invalidate the user-gesture chain if the page awaits before calling
  // navigator.clipboard.write — html2canvas easily takes >100ms, which
  // breaks the gesture window and triggers NotAllowedError. Workaround:
  // hand a Promise<Blob> directly to ClipboardItem so the spec sees the
  // write call as a synchronous continuation of the original click.
  const handleCopyCard = useCallback((cardId: string, event: React.MouseEvent) => {
    event.stopPropagation(); // Prevent navigation to compound
    const cardElement = cardRefs.current.get(cardId);
    if (!cardElement) return;

    const blobPromise = (async () => {
      // Wait for web fonts to finish loading — html2canvas otherwise
      // captures with system-font fallback metrics, producing mis-spaced
      // text (e.g. adjacent words overlapping in copied PNGs).
      if (typeof document !== 'undefined' && document.fonts?.ready) {
        await document.fonts.ready;
      }
      const canvas = await html2canvas(cardElement, {
        backgroundColor: '#ffffff',
        scale: 2, // Higher resolution for presentations
        logging: false,
        onclone: pinCaptureFonts,
      });
      return new Promise<Blob>((resolve, reject) => {
        canvas.toBlob((blob: Blob | null) => {
          if (blob) resolve(blob);
          else reject(new Error('Could not convert canvas to PNG'));
        }, 'image/png');
      });
    })();

    navigator.clipboard
      .write([new ClipboardItem({ 'image/png': blobPromise })])
      .then(() => {
        setCopiedCardId(cardId);
        setTimeout(() => setCopiedCardId(null), 2000);
      })
      .catch((err) => {
        console.error('Failed to copy to clipboard:', err);
      });
  }, []);

  // Render-order-stable list of card DOM nodes for the PowerPoint
  // export. The Map is populated by ref callbacks as cards mount, but
  // its iteration order isn't necessarily the current sort order, so
  // walk sortedRows to drive the ordering.
  const handleExportPptx = useCallback(async () => {
    setExportError(null);
    const nodes: HTMLElement[] = [];
    for (const row of sortedRows) {
      const cardId = showBatch ? `${row.compound_id}-${row.batch_id}` : row.compound_id;
      const node = cardRefs.current.get(cardId);
      if (node) nodes.push(node);
    }
    if (nodes.length === 0) {
      setExportError('No cards to export.');
      return;
    }
    try {
      await exportCardsToPptx(nodes, {
        cardsPerSlide: exportCardsPerSlide,
        onProgress: (done, total) => setExportProgress({ done, total }),
      });
    } catch (err) {
      setExportError(err instanceof Error ? err.message : 'Export failed');
    } finally {
      setExportProgress(null);
    }
  }, [sortedRows, showBatch, exportCardsPerSlide]);

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
                  <MenuItem value="compact">Compact</MenuItem>
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
          {/* PowerPoint export — captures every card in the current
              sort order, then bundles into a .pptx download. The
              cards-per-slide selector matches pptx-export's gridFor
              breakpoints (4/6/9/12). */}
          <FormControl size="small" sx={{ minWidth: 90 }}>
            <Select
              value={exportCardsPerSlide}
              onChange={(e) => setExportCardsPerSlide(Number(e.target.value))}
              size="small"
              sx={{ fontSize: '0.875rem' }}
            >
              <MenuItem value={4}>4 / slide</MenuItem>
              <MenuItem value={6}>6 / slide</MenuItem>
              <MenuItem value={9}>9 / slide</MenuItem>
              <MenuItem value={12}>12 / slide</MenuItem>
            </Select>
          </FormControl>
          <Tooltip
            title={
              exportProgress
                ? `Capturing ${exportProgress.done}/${exportProgress.total}...`
                : 'Download as PowerPoint (.pptx)'
            }
          >
            <span>
              <Button
                size="small"
                variant="outlined"
                startIcon={<Slideshow fontSize="small" />}
                onClick={handleExportPptx}
                disabled={exportProgress !== null}
              >
                {exportProgress
                  ? `${exportProgress.done}/${exportProgress.total}`
                  : 'PPTX'}
              </Button>
            </span>
          </Tooltip>
        </Box>
      </Box>
      {exportError && (
        <Typography variant="caption" color="error" sx={{ mb: 1, display: 'block' }}>
          PowerPoint export failed: {exportError}
        </Typography>
      )}

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
              // Force letterSpacing to a unit-less zero on every text node
              // inside the card. MUI's caption / body variants ship em-based
              // letterSpacing in the theme (~0.03em), which html2canvas
              // mis-resolves into overlapping glyphs in the captured PNG —
              // visible as a strike-through on the compound-ID badge and as
              // squashed words ("TargetmEGFR") in the supporting text.
              // Browser-side change is imperceptible at these font sizes.
              '& .MuiTypography-root': { letterSpacing: 0 },
              '& .MuiChip-label': { letterSpacing: 0 },
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
            {cardContent === 'compact' && hasScorecard && scorecardConfig ? (
              <CompactCardBody
                row={row}
                showBatch={showBatch}
                config={scorecardConfig}
                protocols={protocols}
                concentrationDisplay={concentrationDisplay}
              />
            ) : hasScorecard && cardContent !== 'protocols' ? (
              // Two-column header when spider is shown:
              //   left  ≈ 50% : ID chip (top) → target (line below) → structure
              //   right ≈ 50% : spider, vertically centred
              // Each side gets ~half the card width, so neither element is
              // squeezed and there's near-zero white space between them.
              <Box
                sx={{
                  display: 'flex',
                  gap: 2,
                  mb: 1.5,
                  alignItems: 'center',
                }}
              >
                <Box sx={{ flex: 1, minWidth: 0, display: 'flex', flexDirection: 'column', gap: 0.5 }}>
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, flexWrap: 'wrap' }}>
                    <CompoundIdBadge formattedId={row.formatted_id} />
                    {showBatch && row.batch_number != null && (
                      <Typography variant="caption" color="text.secondary">
                        /{row.batch_number}
                      </Typography>
                    )}
                  </Box>
                  {row.target_name && (
                    <Typography variant="caption" color="text.secondary" sx={{ display: 'block' }}>
                      Target: {row.target_name}
                    </Typography>
                  )}
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
                      }}
                    >
                      <Typography variant="caption" color="text.secondary">-</Typography>
                    </Box>
                  )}
                </Box>
                <Box sx={{ flex: 1, minWidth: 0, display: 'flex', justifyContent: 'center' }}>
                  <CompoundSpider config={scorecardConfig} compound={row} size="small" />
                </Box>
              </Box>
            ) : (
              // No spider: keep the original single-row header (structure
              // + name + target) so the card stays compact for chemists
              // who picked Protocols-only or who haven't configured a scorecard.
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
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, flexWrap: 'wrap' }}>
                    <CompoundIdBadge formattedId={row.formatted_id} />
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
              </Box>
            )}

            {/* Identifiers (barcode / supplier ref / aliases) */}
            {showIdentifiers && cardContent !== 'compact' && (
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
            {includeProperties.length > 0 && cardContent !== 'compact' && (
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
                          fontFamily: MONOSPACE_FONT_STACK,
                        }}
                      >
                        {PROPERTY_LABELS[propName]}: {formatPropertyValue(value)}
                      </Typography>
                    </Tooltip>
                  );
                })}
              </Box>
            )}

            {/* Composite scorecard values — replaces raw protocol rows
                when scorecard is configured AND user picked "Both", so the
                key+value section reflects the same composited data the
                spider is plotting (ratios, worst-of, lipinski, etc.). */}
            {cardContent === 'both' && hasScorecard && scorecardConfig && (
              <Box sx={{ pt: 1, borderTop: 1, borderColor: 'divider' }}>
                <ScorecardValuesTable
                  config={scorecardConfig}
                  compound={row}
                  protocols={protocols}
                  concentrationDisplay={concentrationDisplay}
                />
              </Box>
            )}

            {/* Raw per-protocol rows — shown when user picked Protocols
                only, or when no scorecard is configured at all. */}
            {(cardContent === 'protocols' || (!hasScorecard && cardContent !== 'spider')) && (
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
                      sx={{
                        fontFamily: MONOSPACE_FONT_STACK,
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

// ---------------------------------------------------------------------------
// Card-internal compound-ID badge.
//
// MUI's <Chip> renders fine in the browser but trips up html2canvas: its
// flexbox internals + text-overflow:ellipsis cause the captured PNG to
// truncate the label even when there's plenty of room. The hover-popup
// CompoundNameChip (used in tables) isn't useful here anyway because the
// structure is already shown right next to the badge at full size. Inline
// the styling so the captured DOM is dead simple.
// ---------------------------------------------------------------------------
function CompoundIdBadge({ formattedId }: { formattedId: string }) {
  return (
    <Box
      component="div"
      sx={{
        display: 'inline-flex',
        alignItems: 'center',
        gap: 0.5,
        px: 1,
        py: 0.25,
        border: '1px solid',
        borderColor: 'primary.main',
        borderRadius: 999,
        color: 'primary.main',
        fontFamily: MONOSPACE_FONT_STACK,
        fontSize: '0.875rem',
        letterSpacing: 0,
        whiteSpace: 'nowrap',
        lineHeight: 1.4,
        // Defensive: force NBSP so any whitespace inside the formatted
        // ID survives html2canvas, and disable kerning so the dash in
        // NCL-00030851 doesn't visually merge with the adjacent digits.
        fontKerning: 'none',
      }}
    >
      <Medication sx={{ fontSize: '1rem' }} />
      <span>{formattedId}</span>
    </Box>
  );
}

// ---------------------------------------------------------------------------
// Compact card body: chemistry-first layout. ID chip → large structure →
// vertically stacked two-line bullets (axis label + value inside each
// bar's body track). Skips the identifiers / properties rows so the card
// stays focused on the structure and the scorecard signals.
// ---------------------------------------------------------------------------
function CompactCardBody({
  row,
  showBatch,
  config,
  protocols,
  concentrationDisplay,
}: {
  row: CompactRow;
  showBatch?: boolean;
  config: ScorecardConfig;
  protocols: ProtocolInfo[];
  concentrationDisplay: ConcentrationDisplayMode;
}) {
  // Reorder axes by sector so adjacent bullets share a sector — same
  // ordering the spider, bullets view, and key panel use, so a chemist
  // only learns one mapping.
  const evals = useMemo(() => {
    const raw = evaluateScorecard(config, row);
    if (raw.length === 0) return raw;
    const groups = groupAxesBySector(raw.map((e) => e.axis));
    return groups.flatMap((g) => g.items).map(({ index }) => raw[index]);
  }, [config, row]);
  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: 1 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, flexWrap: 'wrap', justifyContent: 'center' }}>
        <CompoundIdBadge formattedId={row.formatted_id} />
        {showBatch && row.batch_number != null && (
          <Typography variant="caption" color="text.secondary">
            /{row.batch_number}
          </Typography>
        )}
      </Box>
      {row.smiles ? (
        // fillWidth: the molecule expands to the card's available inner
        // width. 4:3 aspect ratio matches the natural shape of most
        // drug-like molecules — square boxes leave conspicuous dead
        // space above and below the chemistry. Generous `size` keeps
        // the SVG's intrinsic label / bond rendering crisp at the
        // larger display dimension; vector content scales for free.
        <Box sx={{ width: '100%' }}>
          <MoleculeChip smiles={row.smiles} size={400} fillWidth aspectRatio={4 / 3} />
        </Box>
      ) : (
        <Box
          sx={{
            width: '100%',
            aspectRatio: '4 / 3',
            bgcolor: 'grey.100',
            borderRadius: 1,
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
          }}
        >
          <Typography variant="caption" color="text.secondary">-</Typography>
        </Box>
      )}
      {evals.length > 0 && (
        <Box sx={{ width: '100%', display: 'flex', flexDirection: 'column', gap: 0.5, mt: 0.5 }}>
          {evals.map((evaluation, i) => (
            <StackedBulletCell
              key={i}
              evaluation={evaluation}
              protocols={protocols}
              concentrationDisplay={concentrationDisplay}
            />
          ))}
        </Box>
      )}
    </Box>
  );
}

// ---------------------------------------------------------------------------
// Two-line bullet for the compact card mode.
//
// Bar **length** encodes compliance (t, the normalised score), bar
// **colour** encodes sector — same palette as the spider wedges and the
// scorecard key panel, so a chemist only learns one colour mapping.
// Track is faint grey; fill is the sector colour at 60% alpha, light
// enough to keep the dark text inside the bar readable across all
// sectors (potency-blue and pk-pink are darkest).
//
// A single-row native <table> places label and value side by side with
// label right-aligned and value left-aligned, meeting at the bar's
// midline. Native tables capture cleanly under html2canvas (CSS Grid
// and nested flex are both prior pain points); since the layout is
// horizontal, font sizes can grow without making the bar taller.
// ---------------------------------------------------------------------------

const BULLET_TRACK_BG = 'rgba(0, 0, 0, 0.05)';
const BULLET_FILL_ALPHA = 0.6;
const BULLET_HEIGHT = 32;

function hexToRgba(hex: string, alpha: number): string {
  const m = hex.replace('#', '');
  const expanded = m.length === 3 ? m.split('').map((c) => c + c).join('') : m;
  const bigint = parseInt(expanded, 16);
  const r = (bigint >> 16) & 255;
  const g = (bigint >> 8) & 255;
  const b = bigint & 255;
  return `rgba(${r}, ${g}, ${b}, ${alpha})`;
}

const BULLET_LABEL_CELL_STYLE: React.CSSProperties = {
  width: '50%',
  textAlign: 'right',
  verticalAlign: 'middle',
  paddingRight: 6,
  fontFamily: SANS_SERIF_FONT_STACK,
  fontSize: '0.85rem',
  fontWeight: 600,
  lineHeight: 1.1,
  letterSpacing: 0,
  fontKerning: 'none',
  color: 'rgba(0,0,0,0.75)',
  whiteSpace: 'nowrap',
};

const BULLET_VALUE_CELL_STYLE: React.CSSProperties = {
  width: '50%',
  textAlign: 'left',
  verticalAlign: 'middle',
  paddingLeft: 6,
  fontFamily: MONOSPACE_FONT_STACK,
  fontSize: '0.9rem',
  fontWeight: 700,
  lineHeight: 1.1,
  letterSpacing: 0,
  fontKerning: 'none',
  color: 'rgba(0,0,0,0.9)',
  whiteSpace: 'nowrap',
};

function StackedBulletCell({
  evaluation,
  protocols,
  concentrationDisplay,
}: {
  evaluation: AxisEvaluation;
  protocols: ProtocolInfo[];
  concentrationDisplay: ConcentrationDisplayMode;
}) {
  const { axis, value, t } = evaluation;
  const label = axis.label || '(unnamed)';

  if (t == null) {
    return (
      <Tooltip title="No data" arrow>
        <Box
          sx={{
            height: BULLET_HEIGHT,
            border: '1px dashed',
            borderColor: 'divider',
            borderRadius: 0.5,
          }}
        >
          <table style={{ width: '100%', height: '100%', borderCollapse: 'collapse', tableLayout: 'fixed' }}>
            <tbody>
              <tr>
                <td style={BULLET_LABEL_CELL_STYLE}>{label}</td>
                <td style={{ ...BULLET_VALUE_CELL_STYLE, color: 'rgba(0,0,0,0.4)', fontWeight: 400 }}>no data</td>
              </tr>
            </tbody>
          </table>
        </Box>
      </Tooltip>
    );
  }

  const fillColour = hexToRgba(sectorColour(axis.sector), BULLET_FILL_ALPHA);
  const display = formatAxisValueForBullet(axis, value, protocols, concentrationDisplay);
  return (
    <Tooltip title={`${display}  ${tierLabel(t)}`} arrow>
      <Box
        sx={{
          position: 'relative',
          height: BULLET_HEIGHT,
          backgroundColor: BULLET_TRACK_BG,
          border: 1,
          borderColor: 'divider',
          borderRadius: 0.5,
          overflow: 'hidden',
        }}
      >
        <Box
          sx={{
            position: 'absolute',
            inset: 0,
            width: `${(t * 100).toFixed(1)}%`,
            backgroundColor: fillColour,
          }}
        />
        <table
          style={{
            position: 'relative',
            width: '100%',
            height: '100%',
            borderCollapse: 'collapse',
            tableLayout: 'fixed',
          }}
        >
          <tbody>
            <tr>
              <td style={BULLET_LABEL_CELL_STYLE}>{label}</td>
              <td style={BULLET_VALUE_CELL_STYLE}>{display}</td>
            </tr>
          </tbody>
        </table>
      </Box>
    </Tooltip>
  );
}
