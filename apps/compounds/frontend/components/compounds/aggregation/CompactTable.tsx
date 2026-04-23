'use client';

import { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import {
  Box,
  Typography,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TableSortLabel,
  Button,
  Chip,
  Tooltip,
  IconButton,
} from '@mui/material';
import { useVirtualizer } from '@tanstack/react-virtual';
import { Download, ZoomIn, EditOutlined, BubbleChart } from '@mui/icons-material';
import { useRouter } from 'next/navigation';
import {
  AggregationType,
  CompactAggregationResponse,
  CompactRow,
  ProtocolInfo,
  ConcentrationDisplayMode,
  MolecularPropertyThreshold,
  MolecularPropertyValues,
  getRagStatus,
} from '@/types/compounds/aggregation';
import { MoleculeChip, CompoundNameChip } from '../MoleculeView';
import { DataSeriesDetailModal } from '../DataSeriesDetailModal';
import {
  generateCompactCsv,
  downloadCsv,
  formatConcentrationValue,
  getConcentrationHeaderUnit,
  isConcentrationUnit,
} from '@/lib/compounds/aggregation-api';
import { protocolColour } from '@/lib/compounds/protocol-colour';
import {
  Order,
  PROPERTY_LABELS,
  RAG_COLORS,
  IdentifiersCell,
  EMPTY_THRESHOLDS,
  formatPropertyValue,
  formatAggLabel,
  getComparator,
  getMeasurementStatus,
  TestedNoValidBadge,
  NotTestedIndicator,
} from './shared';

/** Sort key for compact table - either a fixed column or protocol-based */
type CompactSortKey = 'compound' | 'batch' | 'target' | `protocol_${string}_${AggregationType}`;

export function CompactTable({
  data,
  aggregations,
  concentrationDisplay = 'natural',
  searchTerm = '',
  fillHeight = false,
  onEditProtocol,
  onScatterProtocol,
}: {
  data: CompactAggregationResponse;
  aggregations: AggregationType[];
  concentrationDisplay: ConcentrationDisplayMode;
  searchTerm?: string;
  fillHeight?: boolean;
  onEditProtocol?: (protocol: ProtocolInfo) => void;
  onScatterProtocol?: (protocol: ProtocolInfo) => void;
}) {
  const router = useRouter();
  const parentRef = useRef<HTMLDivElement>(null);
  const topScrollRef = useRef<HTMLDivElement>(null);

  // Modal state
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
  const [orderBy, setOrderBy] = useState<CompactSortKey>('compound');
  const [order, setOrder] = useState<Order>('asc');

  const rows = data.data;
  const protocols = data.protocols;

  // Filter rows by search term
  const filteredRows = useMemo(() => {
    if (!searchTerm) return rows;
    const term = searchTerm.toLowerCase();
    return rows.filter(row => row.formatted_id?.toLowerCase().includes(term));
  }, [rows, searchTerm]);

  // Extract property info from response
  const includeProperties = data.meta.include_properties || [];
  const showIdentifiersColumn = Boolean(data.meta.include_identifiers);
  const propertyThresholds = data.property_thresholds ?? EMPTY_THRESHOLDS;

  // Create a map for quick threshold lookup
  const thresholdMap = useMemo(() => {
    const map: Record<string, MolecularPropertyThreshold> = {};
    for (const t of propertyThresholds) {
      map[t.property_name] = t;
    }
    return map;
  }, [propertyThresholds]);

  // Handle sort request
  const handleRequestSort = (property: CompactSortKey) => {
    const isAsc = orderBy === property && order === 'asc';
    setOrder(isAsc ? 'desc' : 'asc');
    setOrderBy(property);
  };

  // Get value for sorting based on orderBy key
  const getSortValue = useCallback((row: CompactRow, key: CompactSortKey): unknown => {
    if (key === 'compound') return row.formatted_id;
    if (key === 'batch') return row.batch_number ?? -1;
    if (key === 'target') return row.target_name ?? '';

    // Protocol column: protocol_{id}_{aggType}
    const match = key.match(/^protocol_(.+)_(\w+)$/);
    if (match) {
      const [, protocolId, aggType] = match;
      const protocolData = row.protocols[protocolId];
      if (!protocolData) return null;
      const value = protocolData[aggType as AggregationType];
      // For 'list' type, return null to keep stable sort
      if (aggType === 'list') return null;
      return value ?? null;
    }
    return null;
  }, []);

  // Sorted rows (from filtered set)
  const sortedRows = useMemo(() => {
    const comparator = getComparator<CompactRow>(order, (row) => getSortValue(row, orderBy));
    return [...filteredRows].sort(comparator);
  }, [filteredRows, order, orderBy, getSortValue]);
  const showBatchColumn = data.meta.group_by_batch;

  // Track horizontal scroll position to show/hide scroll shadow
  const [canScrollRight, setCanScrollRight] = useState(true);
  // Track table width for top scrollbar
  const [tableScrollWidth, setTableScrollWidth] = useState(0);

  const updateScrollShadow = useCallback(() => {
    const el = parentRef.current;
    if (el) {
      const hasMoreRight = el.scrollLeft < el.scrollWidth - el.clientWidth - 1;
      setCanScrollRight(hasMoreRight);
      setTableScrollWidth(el.scrollWidth);
    }
  }, []);

  // Sync horizontal scroll between top scrollbar and table
  const isSyncingScroll = useRef(false);

  const handleTopScroll = useCallback(() => {
    if (isSyncingScroll.current) return;
    const topEl = topScrollRef.current;
    const tableEl = parentRef.current;
    if (topEl && tableEl) {
      isSyncingScroll.current = true;
      tableEl.scrollLeft = topEl.scrollLeft;
      requestAnimationFrame(() => {
        isSyncingScroll.current = false;
      });
    }
  }, []);

  const handleTableScroll = useCallback(() => {
    if (isSyncingScroll.current) return;
    const topEl = topScrollRef.current;
    const tableEl = parentRef.current;
    if (topEl && tableEl) {
      isSyncingScroll.current = true;
      topEl.scrollLeft = tableEl.scrollLeft;
      requestAnimationFrame(() => {
        isSyncingScroll.current = false;
      });
    }
    updateScrollShadow();
  }, [updateScrollShadow]);

  useEffect(() => {
    const el = parentRef.current;
    if (el) {
      // Wait for table to be fully laid out before measuring scroll width
      // setTimeout ensures we run after React's commit phase and browser paint
      const timeoutId = setTimeout(() => {
        updateScrollShadow();
      }, 50);
      el.addEventListener('scroll', handleTableScroll);
      window.addEventListener('resize', updateScrollShadow);
      return () => {
        clearTimeout(timeoutId);
        el.removeEventListener('scroll', handleTableScroll);
        window.removeEventListener('resize', updateScrollShadow);
      };
    }
  }, [handleTableScroll, updateScrollShadow, protocols.length, aggregations.length]);

  // Virtualization for smooth scrolling with large datasets
  const rowVirtualizer = useVirtualizer({
    count: sortedRows.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 76, // Approximate row height with molecule chip
    overscan: 5,
  });

  const handleExport = () => {
    const csv = generateCompactCsv(data, aggregations);
    downloadCsv(csv, `aggregation_compact_${new Date().toISOString().slice(0, 10)}.csv`);
  };

  const handleProtocolCellClick = (
    event: React.MouseEvent,
    row: CompactRow,
    protocol: ProtocolInfo
  ) => {
    event.stopPropagation(); // Prevent row click navigation
    setSelectedCompound({ id: row.compound_id, name: row.formatted_id });
    setSelectedProtocol({ id: protocol.id, name: protocol.name });
    setModalOpen(true);
  };

  const handleCloseModal = () => {
    setModalOpen(false);
    setSelectedCompound(null);
    setSelectedProtocol(null);
  };

  // Calculate expected table width for top scrollbar
  // Use measured scrollWidth if available, otherwise calculate from column count
  const expectedTableWidth = 320 + (showIdentifiersColumn ? 160 : 0) + protocols.length * aggregations.length * 100;
  const topScrollbarWidth = tableScrollWidth > 0 ? tableScrollWidth : expectedTableWidth;

  return (
    <>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
          <Typography variant="body2" color="text.secondary">
            {searchTerm ? `${filteredRows.length} of ` : ''}{data.meta.compound_count} compounds
            {data.meta.group_by_batch && data.meta.row_count && data.meta.row_count !== data.meta.compound_count && (
              <> ({data.meta.row_count} rows with batch split)</>
            )}
            , {data.meta.protocol_count} protocols,{' '}
            {data.meta.total_measurements} measurements
          </Typography>
          <Tooltip title="Click any protocol cell to view dose-response charts">
            <Chip
              icon={<ZoomIn fontSize="small" />}
              label="Click cell for details"
              size="small"
              variant="outlined"
              color="info"
            />
          </Tooltip>
        </Box>
        <Button startIcon={<Download />} onClick={handleExport} size="small">
          Export CSV
        </Button>
      </Box>

      <Box sx={{ position: 'relative', width: '100%', ...(fillHeight && { flex: 1, display: 'flex', flexDirection: 'column', minHeight: 0 }) }}>
        {/* Top horizontal scrollbar - synced with table */}
        <Box
          ref={topScrollRef}
          onScroll={handleTopScroll}
          sx={{
            overflowX: 'auto',
            overflowY: 'hidden',
            width: '100%',
            mb: 1,
            // Style the scrollbar to be more visible
            '&::-webkit-scrollbar': {
              height: 10,
            },
            '&::-webkit-scrollbar-track': {
              bgcolor: 'grey.100',
              borderRadius: 1,
            },
            '&::-webkit-scrollbar-thumb': {
              bgcolor: 'grey.400',
              borderRadius: 1,
              '&:hover': {
                bgcolor: 'grey.500',
              },
            },
          }}
        >
          {/* Invisible content to create scroll width */}
          <Box sx={{ width: topScrollbarWidth, height: 1 }} />
        </Box>
        {/* Right scroll shadow indicator */}
        <Box
          sx={{
            position: 'absolute',
            top: 44, // Account for top scrollbar
            right: 0,
            bottom: 0,
            width: 40,
            background: 'linear-gradient(to right, transparent, rgba(0,0,0,0.08))',
            pointerEvents: 'none',
            zIndex: 2,
            opacity: canScrollRight ? 1 : 0,
            transition: 'opacity 0.2s ease',
          }}
        />
        <TableContainer
          ref={parentRef}
          sx={fillHeight ? { flex: 1, minHeight: 0, overflow: 'auto' } : { maxHeight: 600, overflow: 'auto' }}
        >
          {/* Calculate min width: 3 fixed columns (320px) + property columns (70px each) + protocol columns (100px each) */}
          <Table
          stickyHeader
          size="small"
          sx={{
            tableLayout: 'fixed',
            minWidth: 320 + includeProperties.length * 70 + protocols.length * aggregations.length * 100,
          }}
        >
          <TableHead>
            <TableRow>
              <TableCell sx={{ fontWeight: 600, width: 100 }}>Structure</TableCell>
              <TableCell sx={{ fontWeight: 600, width: 120 }} sortDirection={orderBy === 'compound' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'compound'}
                  direction={orderBy === 'compound' ? order : 'asc'}
                  onClick={() => handleRequestSort('compound')}
                >
                  Compound
                </TableSortLabel>
              </TableCell>
              {showBatchColumn && (
                <TableCell sx={{ fontWeight: 600, width: 60 }} sortDirection={orderBy === 'batch' ? order : false}>
                  <TableSortLabel
                    active={orderBy === 'batch'}
                    direction={orderBy === 'batch' ? order : 'asc'}
                    onClick={() => handleRequestSort('batch')}
                  >
                    Batch
                  </TableSortLabel>
                </TableCell>
              )}
              {showIdentifiersColumn && (
                <TableCell sx={{ fontWeight: 600, width: 160 }}>Identifiers</TableCell>
              )}
              <TableCell sx={{ fontWeight: 600, width: 100 }} sortDirection={orderBy === 'target' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'target'}
                  direction={orderBy === 'target' ? order : 'asc'}
                  onClick={() => handleRequestSort('target')}
                >
                  Target
                </TableSortLabel>
              </TableCell>
              {/* Molecular property columns */}
              {includeProperties.map((propName) => (
                <TableCell
                  key={propName}
                  sx={{ fontWeight: 600, width: 70 }}
                  align="right"
                >
                  <Tooltip title={thresholdMap[propName]?.property_display || propName}>
                    <span>{PROPERTY_LABELS[propName] || propName}</span>
                  </Tooltip>
                </TableCell>
              ))}
              {protocols.map((protocol) => (
                aggregations.map((agg, aggIdx) => {
                  const sortKey: CompactSortKey = `protocol_${protocol.id}_${agg}`;
                  const isSortable = agg !== 'list'; // List columns aren't sortable
                  const isFirstAggForProtocol = aggIdx === 0;
                  return (
                    <TableCell
                      key={`${protocol.id}-${agg}`}
                      sx={{ fontWeight: 600, width: 100 }}
                      align="right"
                      sortDirection={orderBy === sortKey ? order : false}
                    >
                      <Tooltip title={protocol.name}>
                        <span>
                          {isSortable ? (
                            <TableSortLabel
                              active={orderBy === sortKey}
                              direction={orderBy === sortKey ? order : 'asc'}
                              onClick={() => handleRequestSort(sortKey)}
                            >
                              {protocol.name.length > 15
                                ? `${protocol.name.slice(0, 15)}...`
                                : protocol.name}
                            </TableSortLabel>
                          ) : (
                            protocol.name.length > 15
                              ? `${protocol.name.slice(0, 15)}...`
                              : protocol.name
                          )}
                          {onEditProtocol && isFirstAggForProtocol && (
                            <Tooltip title="Edit thresholds for this protocol" arrow>
                              <IconButton
                                size="small"
                                onClick={(e) => {
                                  e.stopPropagation();
                                  onEditProtocol(protocol);
                                }}
                                sx={{ ml: 0.25, p: 0.25, opacity: 0.5, '&:hover': { opacity: 1 } }}
                              >
                                <EditOutlined fontSize="inherit" />
                              </IconButton>
                            </Tooltip>
                          )}
                          {onScatterProtocol && isFirstAggForProtocol && (
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
                          <br />
                          <Typography variant="caption" color="text.secondary">
                            {concentrationDisplay === 'pConc' && agg !== 'count' && agg !== 'list' && isConcentrationUnit(protocol.kpi_unit)
                              ? `p${formatAggLabel(agg)}`
                              : formatAggLabel(agg)}
                            {/* Show unit for value-based aggregations */}
                            {agg !== 'count' && agg !== 'list' && protocol.kpi_unit && concentrationDisplay !== 'pConc' && (
                              <> ({getConcentrationHeaderUnit(protocol.kpi_unit, concentrationDisplay)})</>
                            )}
                          </Typography>
                        </span>
                      </Tooltip>
                    </TableCell>
                  );
                })
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {/* Spacer for rows above viewport */}
            {rowVirtualizer.getVirtualItems().length > 0 &&
              rowVirtualizer.getVirtualItems()[0].start > 0 && (
              <TableRow>
                <TableCell
                  colSpan={3 + (showBatchColumn ? 1 : 0) + (showIdentifiersColumn ? 1 : 0) + includeProperties.length + protocols.length * aggregations.length}
                  sx={{
                    height: rowVirtualizer.getVirtualItems()[0].start,
                    padding: 0,
                    border: 'none',
                  }}
                />
              </TableRow>
            )}

            {/* Virtualized rows */}
            {rowVirtualizer.getVirtualItems().map((virtualRow) => {
              const row = sortedRows[virtualRow.index];
              return (
                <TableRow
                  key={showBatchColumn ? `${row.compound_id}-${row.batch_id || 'no-batch'}` : row.compound_id}
                  data-index={virtualRow.index}
                  ref={rowVirtualizer.measureElement}
                  hover
                  sx={{ cursor: 'pointer' }}
                  onClick={() => router.push(`/registry/compounds/${row.compound_id}`)}
                >
                  <TableCell>
                    {row.smiles ? (
                      <MoleculeChip smiles={row.smiles} size={80} />
                    ) : (
                      <Box
                        sx={{
                          width: 80,
                          height: 80,
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
                  </TableCell>
                  <TableCell>
                    <CompoundNameChip formattedId={row.formatted_id} smiles={row.smiles} label={String(Number(row.formatted_id.split('-').pop()))} />
                  </TableCell>
                  {showBatchColumn && (
                    <TableCell>
                      <Typography variant="body2" fontFamily="monospace">
                        {row.batch_number != null ? `/${row.batch_number}` : '-'}
                      </Typography>
                    </TableCell>
                  )}
                  {showIdentifiersColumn && (
                    <TableCell onClick={(e) => e.stopPropagation()}>
                      <IdentifiersCell identifiers={row.identifiers} />
                    </TableCell>
                  )}
                  <TableCell>
                    <Typography variant="body2">{row.target_name || '-'}</Typography>
                  </TableCell>
                  {/* Molecular property cells */}
                  {includeProperties.map((propName) => {
                    const propValue = row.properties?.[propName as keyof MolecularPropertyValues] as number | null | undefined;
                    const threshold = thresholdMap[propName];
                    const ragStatus = getRagStatus(propValue, threshold);
                    return (
                      <TableCell key={propName} align="right">
                        <Typography
                          variant="body2"
                          fontFamily="monospace"
                          sx={{ color: RAG_COLORS[ragStatus] }}
                        >
                          {formatPropertyValue(propValue)}
                        </Typography>
                      </TableCell>
                    );
                  })}
                  {protocols.map((protocol) => (
                    aggregations.map((agg) => {
                      const protocolData = row.protocols[protocol.id];
                      const value = protocolData?.[agg];
                      const measurementStatus = getMeasurementStatus(protocolData, aggregations.includes('count'));
                      const hasData = measurementStatus === 'has-data';
                      const testedNoValid = measurementStatus === 'tested-no-valid';
                      const notTested = measurementStatus === 'not-measured';

                      // Threshold-based colour only applies to the geomean KPI column
                      const cellColour = agg === 'geomean' && hasData
                        ? protocolColour(protocolData?.geomean, protocol)
                        : { background: null, t: null };

                      return (
                        <TableCell
                          key={`${protocol.id}-${agg}`}
                          align="right"
                          onClick={(e) => hasData && handleProtocolCellClick(e, row, protocol)}
                          sx={{
                            cursor: hasData ? 'zoom-in' : 'default',
                            bgcolor: cellColour.background ?? undefined,
                            '&:hover': hasData
                              ? cellColour.background
                                ? { filter: 'brightness(0.95)' }
                                : { bgcolor: 'action.hover' }
                              : {},
                          }}
                        >
                          {/* Show "not tested" indicator on first column for not-measured */}
                          {notTested && agg === aggregations[0] ? (
                            <NotTestedIndicator />
                          ) : notTested ? (
                            <Typography variant="body2" color="text.disabled">-</Typography>
                          ) : /* Show "tested but no valid" badge only on first aggregation column for this protocol */
                          testedNoValid && agg === aggregations[0] ? (
                            <TestedNoValidBadge
                              tested={protocolData?.tested}
                              noAnalysis={protocolData?.no_analysis}
                              invalid={protocolData?.invalid}
                              unassigned={protocolData?.unassigned}
                              onClick={(e) => handleProtocolCellClick(e, row, protocol)}
                            />
                          ) : testedNoValid ? (
                            <Typography variant="body2" color="text.disabled">-</Typography>
                          ) : agg === 'list' ? (
                            <Tooltip title={value || '-'}>
                              <Typography
                                variant="body2"
                                fontFamily="monospace"
                                sx={{
                                  maxWidth: 100,
                                  overflow: 'hidden',
                                  textOverflow: 'ellipsis',
                                  whiteSpace: 'nowrap',
                                }}
                              >
                                {value || '-'}
                              </Typography>
                            </Tooltip>
                          ) : agg === 'count' ? (
                            <Typography variant="body2" fontFamily="monospace">
                              {value ?? '-'}
                            </Typography>
                          ) : (
                            <Typography variant="body2" fontFamily="monospace">
                              {formatConcentrationValue(
                                value as number | null,
                                protocol.kpi_unit,
                                concentrationDisplay
                              ).displayValue}
                            </Typography>
                          )}
                        </TableCell>
                      );
                    })
                  ))}
                </TableRow>
              );
            })}

            {/* Spacer for rows below viewport */}
            {rowVirtualizer.getVirtualItems().length > 0 && (
              <TableRow>
                <TableCell
                  colSpan={3 + (showBatchColumn ? 1 : 0) + includeProperties.length + protocols.length * aggregations.length}
                  sx={{
                    height:
                      rowVirtualizer.getTotalSize() -
                      (rowVirtualizer.getVirtualItems()[
                        rowVirtualizer.getVirtualItems().length - 1
                      ]?.end ?? 0),
                    padding: 0,
                    border: 'none',
                  }}
                />
              </TableRow>
            )}
          </TableBody>
        </Table>
      </TableContainer>
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
