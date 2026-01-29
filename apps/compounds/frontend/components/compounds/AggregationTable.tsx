'use client';

import { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import {
  Box,
  Paper,
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
  LinearProgress,
  Tooltip,
  FormControl,
  Select,
  MenuItem,
  InputLabel,
} from '@mui/material';
import { useVirtualizer } from '@tanstack/react-virtual';
import { Download, Medication, ZoomIn } from '@mui/icons-material';
import { useRouter } from 'next/navigation';
import {
  AggregationResponse,
  AggregationType,
  CompactRow,
  MediumRow,
  LongRow,
  isCompactResponse,
  isMediumResponse,
  ProtocolInfo,
  ConcentrationDisplayMode,
} from '@/types/compounds/aggregation';
import { MoleculeChip } from './MoleculeView';
import { DataSeriesDetailModal } from './DataSeriesDetailModal';
import { ProtocolScatterPlot } from './ProtocolScatterPlot';
import {
  formatKpiUnit,
  generateCompactCsv,
  generateLongCsv,
  downloadCsv,
  formatConcentrationValue,
  getConcentrationHeaderUnit,
  isConcentrationUnit,
} from '@/lib/compounds/aggregation-api';

/**
 * Generate CSV content from medium aggregation data.
 */
function generateMediumCsv(
  data: AggregationResponse,
  aggregations: AggregationType[]
): string {
  const { data: rows } = data;

  // Build header row (include Unit column)
  const headers = ['Compound ID', 'SMILES', 'Target', 'Protocol', 'Unit'];
  for (const agg of aggregations) {
    headers.push(agg.charAt(0).toUpperCase() + agg.slice(1));
  }

  const csvRows = [headers.join(',')];

  // Build data rows
  for (const row of rows as MediumRow[]) {
    const values = [
      `"${row.formatted_id}"`,
      `"${row.smiles || ''}"`,
      `"${row.target_name || ''}"`,
      `"${row.protocol_name}"`,
      `"${row.kpi_unit ? formatKpiUnit(row.kpi_unit) : ''}"`,
    ];

    for (const agg of aggregations) {
      const value = row[agg as keyof MediumRow];
      if (agg === 'list') {
        values.push(`"${value || ''}"`);
      } else {
        values.push(value !== null && value !== undefined ? String(value) : '');
      }
    }

    csvRows.push(values.join(','));
  }

  return csvRows.join('\n');
}

/** Options for concentration display mode selector */
const CONCENTRATION_DISPLAY_OPTIONS: { value: ConcentrationDisplayMode; label: string }[] = [
  { value: 'natural', label: 'Natural' },
  { value: 'nM', label: 'nM' },
  { value: 'uM', label: 'μM' },
  { value: 'mM', label: 'mM' },
  { value: 'pConc', label: 'pConc' },
];

interface AggregationTableProps {
  data: AggregationResponse | null | undefined;
  loading?: boolean;
  aggregations: AggregationType[];
  /** Concentration display mode (default: 'natural') */
  concentrationDisplay?: ConcentrationDisplayMode;
  /** Callback when concentration display mode changes */
  onConcentrationDisplayChange?: (mode: ConcentrationDisplayMode) => void;
}

/**
 * Format aggregation label for column header.
 */
function formatAggLabel(agg: AggregationType): string {
  const labels: Record<AggregationType, string> = {
    geomean: 'GeoMean',
    count: 'N',
    stdev: 'StDev',
    list: 'Values',
  };
  return labels[agg] || agg;
}

/**
 * Compact table component (one row per compound).
 * Clicking on protocol cells opens a modal with data series details.
 */
/** Sort direction type */
type Order = 'asc' | 'desc';

/** Sort key for compact table - either a fixed column or protocol-based */
type CompactSortKey = 'compound' | 'batch' | 'target' | `protocol_${string}_${AggregationType}`;

/** Comparator function for sorting */
function descendingComparator<T>(a: T, b: T, getValue: (item: T) => unknown): number {
  const aVal = getValue(a);
  const bVal = getValue(b);

  // Handle null/undefined - push to end
  if (aVal == null && bVal == null) return 0;
  if (aVal == null) return 1;
  if (bVal == null) return -1;

  // Compare values
  if (typeof aVal === 'string' && typeof bVal === 'string') {
    return bVal.localeCompare(aVal);
  }
  if (typeof aVal === 'number' && typeof bVal === 'number') {
    return bVal - aVal;
  }
  return 0;
}

function getComparator<T>(
  order: Order,
  getValue: (item: T) => unknown
): (a: T, b: T) => number {
  return order === 'desc'
    ? (a, b) => descendingComparator(a, b, getValue)
    : (a, b) => -descendingComparator(a, b, getValue);
}

function CompactTable({
  data,
  aggregations,
  concentrationDisplay = 'natural',
}: {
  data: AggregationResponse & { protocols: ProtocolInfo[] };
  aggregations: AggregationType[];
  concentrationDisplay: ConcentrationDisplayMode;
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

  const rows = data.data as CompactRow[];
  const protocols = data.protocols;

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

  // Sorted rows
  const sortedRows = useMemo(() => {
    const comparator = getComparator<CompactRow>(order, (row) => getSortValue(row, orderBy));
    return [...rows].sort(comparator);
  }, [rows, order, orderBy, getSortValue]);
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
  const expectedTableWidth = 320 + protocols.length * aggregations.length * 100;
  const topScrollbarWidth = tableScrollWidth > 0 ? tableScrollWidth : expectedTableWidth;

  return (
    <>
      {/* Protocol comparison scatter plot */}
      {aggregations.includes('geomean') && protocols.length >= 2 && (
        <ProtocolScatterPlot data={rows} protocols={protocols} />
      )}

      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
          <Typography variant="body2" color="text.secondary">
            {data.meta.compound_count} compounds
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

      <Box sx={{ position: 'relative', width: '100%' }}>
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
          sx={{ maxHeight: 600, overflow: 'auto' }}
        >
          {/* Calculate min width: 3 fixed columns (320px) + protocol columns (100px each) */}
          <Table
          stickyHeader
          size="small"
          sx={{
            tableLayout: 'fixed',
            minWidth: 320 + protocols.length * aggregations.length * 100,
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
              <TableCell sx={{ fontWeight: 600, width: 100 }} sortDirection={orderBy === 'target' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'target'}
                  direction={orderBy === 'target' ? order : 'asc'}
                  onClick={() => handleRequestSort('target')}
                >
                  Target
                </TableSortLabel>
              </TableCell>
              {protocols.map((protocol) => (
                aggregations.map((agg) => {
                  const sortKey: CompactSortKey = `protocol_${protocol.id}_${agg}`;
                  const isSortable = agg !== 'list'; // List columns aren't sortable
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
                  colSpan={3 + (showBatchColumn ? 1 : 0) + protocols.length * aggregations.length}
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
                    <Chip
                      icon={<Medication fontSize="small" />}
                      label={row.formatted_id}
                      size="small"
                      variant="outlined"
                    />
                  </TableCell>
                  {showBatchColumn && (
                    <TableCell>
                      <Typography variant="body2" fontFamily="monospace">
                        {row.batch_number != null ? `/${row.batch_number}` : '-'}
                      </Typography>
                    </TableCell>
                  )}
                  <TableCell>
                    <Typography variant="body2">{row.target_name || '-'}</Typography>
                  </TableCell>
                  {protocols.map((protocol) => (
                    aggregations.map((agg) => {
                      const protocolData = row.protocols[protocol.id];
                      const value = protocolData?.[agg];
                      const hasData = protocolData && Object.keys(protocolData).length > 0;

                      return (
                        <TableCell
                          key={`${protocol.id}-${agg}`}
                          align="right"
                          onClick={(e) => hasData && handleProtocolCellClick(e, row, protocol)}
                          sx={{
                            cursor: hasData ? 'zoom-in' : 'default',
                            '&:hover': hasData ? {
                              bgcolor: 'action.hover',
                            } : {},
                          }}
                        >
                          {agg === 'list' ? (
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
                  colSpan={3 + (showBatchColumn ? 1 : 0) + protocols.length * aggregations.length}
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

/** Sort key for medium table */
type MediumSortKey = 'compound' | 'batch' | 'target' | 'protocol' | AggregationType;

/**
 * Medium table component (one row per compound-protocol pair).
 * Clicking a row opens a modal with data series details and charts.
 */
function MediumTable({
  data,
  aggregations,
  concentrationDisplay = 'natural',
}: {
  data: AggregationResponse;
  aggregations: AggregationType[];
  concentrationDisplay: ConcentrationDisplayMode;
}) {
  const parentRef = useRef<HTMLDivElement>(null);

  // Modal state
  const [modalOpen, setModalOpen] = useState(false);
  const [selectedRow, setSelectedRow] = useState<MediumRow | null>(null);

  // Sorting state
  const [orderBy, setOrderBy] = useState<MediumSortKey>('compound');
  const [order, setOrder] = useState<Order>('asc');

  const rows = data.data as MediumRow[];
  const showBatchColumn = data.meta.group_by_batch;

  // Handle sort request
  const handleRequestSort = (property: MediumSortKey) => {
    const isAsc = orderBy === property && order === 'asc';
    setOrder(isAsc ? 'desc' : 'asc');
    setOrderBy(property);
  };

  // Get value for sorting based on orderBy key
  const getSortValue = useCallback((row: MediumRow, key: MediumSortKey): unknown => {
    if (key === 'compound') return row.formatted_id;
    if (key === 'batch') return row.batch_number ?? -1;
    if (key === 'target') return row.target_name ?? '';
    if (key === 'protocol') return row.protocol_name;
    // Aggregation columns
    if (key === 'list') return null; // Not sortable
    const value = row[key as keyof MediumRow];
    return value ?? null;
  }, []);

  // Sorted rows
  const sortedRows = useMemo(() => {
    const comparator = getComparator<MediumRow>(order, (row) => getSortValue(row, orderBy));
    return [...rows].sort(comparator);
  }, [rows, order, orderBy, getSortValue]);

  // Virtualization for smooth scrolling with large datasets
  const rowVirtualizer = useVirtualizer({
    count: sortedRows.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 76, // Approximate row height with molecule chip
    overscan: 5,
  });

  const handleExport = () => {
    const csv = generateMediumCsv(data, aggregations);
    downloadCsv(csv, `aggregation_medium_${new Date().toISOString().slice(0, 10)}.csv`);
  };

  const handleRowClick = (row: MediumRow) => {
    setSelectedRow(row);
    setModalOpen(true);
  };

  const handleCloseModal = () => {
    setModalOpen(false);
    setSelectedRow(null);
  };

  return (
    <>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
          <Typography variant="body2" color="text.secondary">
            {data.meta.compound_count} compounds
            {data.meta.group_by_batch && data.meta.row_count && data.meta.row_count !== data.meta.compound_count && (
              <> ({data.meta.row_count} rows with batch split)</>
            )}
            , {data.meta.protocol_count} protocols,{' '}
            {data.meta.total_measurements} measurements
          </Typography>
          <Tooltip title="Click any row to view dose-response charts">
            <Chip
              icon={<ZoomIn fontSize="small" />}
              label="Click row for details"
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

      <TableContainer
        ref={parentRef}
        sx={{ maxHeight: 600, overflow: 'auto' }}
      >
        <Table stickyHeader size="small" sx={{ tableLayout: 'fixed' }}>
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
              <TableCell sx={{ fontWeight: 600, width: 100 }} sortDirection={orderBy === 'target' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'target'}
                  direction={orderBy === 'target' ? order : 'asc'}
                  onClick={() => handleRequestSort('target')}
                >
                  Target
                </TableSortLabel>
              </TableCell>
              <TableCell sx={{ fontWeight: 600, width: 120 }} sortDirection={orderBy === 'protocol' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'protocol'}
                  direction={orderBy === 'protocol' ? order : 'asc'}
                  onClick={() => handleRequestSort('protocol')}
                >
                  Protocol
                </TableSortLabel>
              </TableCell>
              {aggregations.map((agg) => {
                const isSortable = agg !== 'list';
                return (
                  <TableCell
                    key={agg}
                    sx={{ fontWeight: 600, width: 80 }}
                    align="right"
                    sortDirection={orderBy === agg ? order : false}
                  >
                    {isSortable ? (
                      <TableSortLabel
                        active={orderBy === agg}
                        direction={orderBy === agg ? order : 'asc'}
                        onClick={() => handleRequestSort(agg)}
                      >
                        {formatAggLabel(agg)}
                      </TableSortLabel>
                    ) : (
                      formatAggLabel(agg)
                    )}
                  </TableCell>
                );
              })}
            </TableRow>
          </TableHead>
          <TableBody>
            {/* Spacer for rows above viewport */}
            {rowVirtualizer.getVirtualItems().length > 0 &&
              rowVirtualizer.getVirtualItems()[0].start > 0 && (
              <TableRow>
                <TableCell
                  colSpan={4 + (showBatchColumn ? 1 : 0) + aggregations.length}
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
                  key={`${row.compound_id}-${row.batch_id || 'no-batch'}-${row.protocol_id}-${virtualRow.index}`}
                  data-index={virtualRow.index}
                  ref={rowVirtualizer.measureElement}
                  hover
                  sx={{ cursor: 'pointer' }}
                  onClick={() => handleRowClick(row)}
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
                    <Chip
                      icon={<Medication fontSize="small" />}
                      label={row.formatted_id}
                      size="small"
                      variant="outlined"
                    />
                  </TableCell>
                  {showBatchColumn && (
                    <TableCell>
                      <Typography variant="body2" fontFamily="monospace">
                        {row.batch_number != null ? `/${row.batch_number}` : '-'}
                      </Typography>
                    </TableCell>
                  )}
                  <TableCell>
                    <Typography variant="body2">{row.target_name || '-'}</Typography>
                  </TableCell>
                  <TableCell>
                    <Typography variant="body2">{row.protocol_name}</Typography>
                  </TableCell>
                  {aggregations.map((agg) => {
                    const value = row[agg as keyof MediumRow];
                    const formatted = formatConcentrationValue(
                      value as number | null,
                      row.kpi_unit,
                      concentrationDisplay
                    );
                    return (
                      <TableCell key={agg} align="right">
                        {agg === 'list' ? (
                          <Tooltip title={String(value || '-')}>
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
                              {String(value || '-')}
                            </Typography>
                          </Tooltip>
                        ) : agg === 'count' ? (
                          <Typography variant="body2" fontFamily="monospace">
                            {value ?? '-'}
                          </Typography>
                        ) : (
                          <Tooltip title={formatted.displayUnit || 'No unit'}>
                            <Typography variant="body2" fontFamily="monospace">
                              {formatted.displayValue}
                              {formatted.displayUnit && value != null && concentrationDisplay === 'natural' && (
                                <Typography
                                  component="span"
                                  variant="caption"
                                  color="text.secondary"
                                  sx={{ ml: 0.5 }}
                                >
                                  {formatted.displayUnit}
                                </Typography>
                              )}
                            </Typography>
                          </Tooltip>
                        )}
                      </TableCell>
                    );
                  })}
                </TableRow>
              );
            })}

            {/* Spacer for rows below viewport */}
            {rowVirtualizer.getVirtualItems().length > 0 && (
              <TableRow>
                <TableCell
                  colSpan={4 + (showBatchColumn ? 1 : 0) + aggregations.length}
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

      {/* Data series detail modal */}
      {selectedRow && (
        <DataSeriesDetailModal
          open={modalOpen}
          onClose={handleCloseModal}
          compoundId={selectedRow.compound_id}
          protocolId={selectedRow.protocol_id}
          compoundName={selectedRow.formatted_id}
          protocolName={selectedRow.protocol_name}
        />
      )}
    </>
  );
}

/** Sort key for long table */
type LongSortKey = 'compound' | 'batch' | 'target' | 'protocol' | 'date' | 'kpi' | 'status';

/**
 * Long table component (one row per measurement).
 */
function LongTable({
  data,
  concentrationDisplay = 'natural',
}: {
  data: AggregationResponse;
  concentrationDisplay: ConcentrationDisplayMode;
}) {
  const router = useRouter();
  const parentRef = useRef<HTMLDivElement>(null);

  // Sorting state
  const [orderBy, setOrderBy] = useState<LongSortKey>('compound');
  const [order, setOrder] = useState<Order>('asc');

  const rows = data.data as LongRow[];
  // Long format always includes batch info, but only show column if any row has batch data
  const showBatchColumn = rows.some(row => row.batch_number != null);

  // Handle sort request
  const handleRequestSort = (property: LongSortKey) => {
    const isAsc = orderBy === property && order === 'asc';
    setOrder(isAsc ? 'desc' : 'asc');
    setOrderBy(property);
  };

  // Get value for sorting based on orderBy key
  const getSortValue = useCallback((row: LongRow, key: LongSortKey): unknown => {
    if (key === 'compound') return row.formatted_id ?? row.compound_name ?? '';
    if (key === 'batch') return row.batch_number ?? -1;
    if (key === 'target') return row.target_name ?? '';
    if (key === 'protocol') return row.protocol_name;
    if (key === 'date') return row.assay_date ?? '';
    if (key === 'kpi') return row.kpi_value ?? null;
    if (key === 'status') return row.status ?? '';
    return null;
  }, []);

  // Sorted rows
  const sortedRows = useMemo(() => {
    const comparator = getComparator<LongRow>(order, (row) => getSortValue(row, orderBy));
    return [...rows].sort(comparator);
  }, [rows, order, orderBy, getSortValue]);

  // Virtualization for smooth scrolling with large datasets
  const rowVirtualizer = useVirtualizer({
    count: sortedRows.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 76, // Approximate row height with molecule chip
    overscan: 5,
  });

  const handleExport = () => {
    const csv = generateLongCsv(data);
    downloadCsv(csv, `aggregation_long_${new Date().toISOString().slice(0, 10)}.csv`);
  };

  return (
    <>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
        <Typography variant="body2" color="text.secondary">
          {data.meta.compound_count} compounds, {data.meta.protocol_count} protocols,{' '}
          {data.meta.total_measurements} measurements
        </Typography>
        <Button startIcon={<Download />} onClick={handleExport} size="small">
          Export CSV
        </Button>
      </Box>

      <TableContainer
        ref={parentRef}
        sx={{ maxHeight: 600, overflow: 'auto' }}
      >
        <Table stickyHeader size="small" sx={{ tableLayout: 'fixed' }}>
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
              <TableCell sx={{ fontWeight: 600, width: 100 }} sortDirection={orderBy === 'target' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'target'}
                  direction={orderBy === 'target' ? order : 'asc'}
                  onClick={() => handleRequestSort('target')}
                >
                  Target
                </TableSortLabel>
              </TableCell>
              <TableCell sx={{ fontWeight: 600, width: 120 }} sortDirection={orderBy === 'protocol' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'protocol'}
                  direction={orderBy === 'protocol' ? order : 'asc'}
                  onClick={() => handleRequestSort('protocol')}
                >
                  Protocol
                </TableSortLabel>
              </TableCell>
              <TableCell sx={{ fontWeight: 600, width: 100 }} sortDirection={orderBy === 'date' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'date'}
                  direction={orderBy === 'date' ? order : 'asc'}
                  onClick={() => handleRequestSort('date')}
                >
                  Date
                </TableSortLabel>
              </TableCell>
              <TableCell sx={{ fontWeight: 600, width: 80 }} align="right" sortDirection={orderBy === 'kpi' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'kpi'}
                  direction={orderBy === 'kpi' ? order : 'asc'}
                  onClick={() => handleRequestSort('kpi')}
                >
                  KPI
                </TableSortLabel>
              </TableCell>
              <TableCell sx={{ fontWeight: 600, width: 80 }} sortDirection={orderBy === 'status' ? order : false}>
                <TableSortLabel
                  active={orderBy === 'status'}
                  direction={orderBy === 'status' ? order : 'asc'}
                  onClick={() => handleRequestSort('status')}
                >
                  Status
                </TableSortLabel>
              </TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {/* Spacer for rows above viewport */}
            {rowVirtualizer.getVirtualItems().length > 0 &&
              rowVirtualizer.getVirtualItems()[0].start > 0 && (
              <TableRow>
                <TableCell
                  colSpan={7 + (showBatchColumn ? 1 : 0)}
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
                  key={row.data_series_id}
                  data-index={virtualRow.index}
                  ref={rowVirtualizer.measureElement}
                  hover
                  sx={{ cursor: 'pointer' }}
                  onClick={() => row.compound_id && router.push(`/registry/compounds/${row.compound_id}`)}
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
                    {row.formatted_id ? (
                      <Chip
                        icon={<Medication fontSize="small" />}
                        label={row.formatted_id}
                        size="small"
                        variant="outlined"
                      />
                    ) : (
                      <Typography variant="body2" color="text.secondary">
                        {row.compound_name || 'Unknown'}
                      </Typography>
                    )}
                  </TableCell>
                  {showBatchColumn && (
                    <TableCell>
                      <Typography variant="body2" fontFamily="monospace">
                        {row.batch_number != null ? `/${row.batch_number}` : '-'}
                      </Typography>
                    </TableCell>
                  )}
                  <TableCell>
                    <Typography variant="body2">{row.target_name || '-'}</Typography>
                  </TableCell>
                  <TableCell>
                    <Typography variant="body2">{row.protocol_name}</Typography>
                  </TableCell>
                  <TableCell>
                    <Typography variant="body2">
                      {row.assay_date ? new Date(row.assay_date).toLocaleDateString() : '-'}
                    </Typography>
                  </TableCell>
                  <TableCell align="right">
                    {(() => {
                      const formatted = formatConcentrationValue(
                        row.kpi_value,
                        row.kpi_unit,
                        concentrationDisplay
                      );
                      return (
                        <Tooltip title={formatted.displayUnit || 'No unit'}>
                          <Typography variant="body2" fontFamily="monospace" fontWeight={500}>
                            {formatted.displayValue}
                            {formatted.displayUnit && row.kpi_value != null && concentrationDisplay === 'natural' && (
                              <Typography
                                component="span"
                                variant="caption"
                                color="text.secondary"
                                sx={{ ml: 0.5 }}
                              >
                                {formatted.displayUnit}
                              </Typography>
                            )}
                          </Typography>
                        </Tooltip>
                      );
                    })()}
                  </TableCell>
                  <TableCell>
                    <Chip
                      label={row.status || 'unknown'}
                      size="small"
                      color={row.status === 'valid' ? 'success' : row.status === 'invalid' ? 'error' : 'default'}
                    />
                  </TableCell>
                </TableRow>
              );
            })}

            {/* Spacer for rows below viewport */}
            {rowVirtualizer.getVirtualItems().length > 0 && (
              <TableRow>
                <TableCell
                  colSpan={7 + (showBatchColumn ? 1 : 0)}
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
    </>
  );
}

/**
 * Concentration display mode selector component.
 */
function ConcentrationDisplaySelector({
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
 * Main aggregation table component.
 * Renders either compact or long format based on the response data.
 */
export function AggregationTable({
  data,
  loading,
  aggregations,
  concentrationDisplay = 'natural',
  onConcentrationDisplayChange,
}: AggregationTableProps) {
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

  if (loading) {
    return (
      <Paper sx={{ p: 3 }}>
        <LinearProgress />
        <Typography sx={{ mt: 2 }} color="text.secondary" align="center">
          Running aggregation query...
        </Typography>
      </Paper>
    );
  }

  if (!data) {
    return (
      <Paper sx={{ p: 3 }}>
        <Typography color="text.secondary" align="center">
          Select a target, protocol, or compound above to see results.
        </Typography>
      </Paper>
    );
  }

  if (data.meta.total_measurements === 0) {
    return (
      <Paper sx={{ p: 3 }}>
        <Typography color="text.secondary" align="center">
          No data found matching your criteria.
        </Typography>
      </Paper>
    );
  }

  // Determine which table to render based on response type
  const renderTable = () => {
    if (isCompactResponse(data)) {
      return (
        <CompactTable
          data={data}
          aggregations={aggregations}
          concentrationDisplay={displayMode}
        />
      );
    } else if (isMediumResponse(data)) {
      return (
        <MediumTable
          data={data}
          aggregations={aggregations}
          concentrationDisplay={displayMode}
        />
      );
    } else {
      return <LongTable data={data} concentrationDisplay={displayMode} />;
    }
  };

  return (
    <Paper sx={{ width: '100%', overflow: 'hidden' }}>
      <Box sx={{ p: 2 }}>
        <Box sx={{ display: 'flex', justifyContent: 'flex-end', mb: 1 }}>
          <ConcentrationDisplaySelector
            value={displayMode}
            onChange={handleDisplayChange}
          />
        </Box>
        {renderTable()}
      </Box>
    </Paper>
  );
}
