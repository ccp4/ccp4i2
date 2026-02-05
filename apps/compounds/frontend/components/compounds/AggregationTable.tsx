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
  MolecularPropertyName,
  MolecularPropertyThreshold,
  MolecularPropertyValues,
  getRagStatus,
  OutputFormat,
  CompactAggregationResponse,
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

/** Short labels for molecular properties */
const PROPERTY_LABELS: Record<MolecularPropertyName, string> = {
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
const RAG_COLORS: Record<string, string> = {
  green: 'inherit',
  amber: '#ff9800',
  red: '#f44336',
};

/** Format property value for display */
function formatPropertyValue(value: number | string | null | undefined): string {
  if (value == null || value === '') return '-';
  // Convert string values to numbers (handles legacy data stored as strings)
  const numValue = typeof value === 'string' ? parseFloat(value) : value;
  if (isNaN(numValue)) return String(value);
  if (Number.isInteger(numValue)) return String(numValue);
  return numValue.toFixed(2);
}

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
  /** If true, table fills available parent height instead of using fixed maxHeight */
  fillHeight?: boolean;
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

function getComparator<T>(
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

  // Extract property info from response
  const includeProperties = data.meta.include_properties || [];
  const propertyThresholds = (data as { property_thresholds?: MolecularPropertyThreshold[] }).property_thresholds || [];

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
                  colSpan={3 + (showBatchColumn ? 1 : 0) + includeProperties.length + protocols.length * aggregations.length}
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
                  {aggregations.map((agg, aggIndex) => {
                    const value = row[agg as keyof MediumRow];
                    const formatted = formatConcentrationValue(
                      value as number | null,
                      row.kpi_unit,
                      concentrationDisplay
                    );
                    // Check if this is a "tested but no valid" row (count === 0)
                    const isTestedNoValid = aggregations.includes('count') && row.count === 0;

                    return (
                      <TableCell key={agg} align="right">
                        {/* Show badge on first column for tested-no-valid rows */}
                        {isTestedNoValid && aggIndex === 0 ? (
                          <TestedNoValidBadge
                            tested={row.tested}
                            noAnalysis={row.no_analysis}
                            invalid={row.invalid}
                            unassigned={row.unassigned}
                            onClick={(e) => {
                              e.stopPropagation();
                              handleRowClick(row);
                            }}
                          />
                        ) : isTestedNoValid ? (
                          <Typography variant="body2" color="text.disabled">-</Typography>
                        ) : agg === 'list' ? (
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
                            {String(value ?? '-')}
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
function getMeasurementStatus(
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
function formatMeasurementWithStats(
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
function TestedNoValidBadge({
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
function NotTestedIndicator() {
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

/** Sort key for pivot table - sort compounds (columns) by property or protocol */
type PivotSortKey = 'compound' | `property_${MolecularPropertyName}` | `protocol_${string}`;

/**
 * Pivot table component (compounds as columns, properties/protocols as rows).
 * Horizontally scrollable for many compounds. Click row labels to sort columns.
 */
function PivotTable({
  data,
  aggregations,
  concentrationDisplay = 'natural',
}: {
  data: CompactAggregationResponse;
  aggregations: AggregationType[];
  concentrationDisplay: ConcentrationDisplayMode;
}) {
  const router = useRouter();
  const scrollContainerRef = useRef<HTMLDivElement>(null);

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

  // Sorting state - sort compounds (columns) by a property or protocol value
  const [orderBy, setOrderBy] = useState<PivotSortKey>('compound');
  const [order, setOrder] = useState<Order>('asc');

  const rows = data.data as CompactRow[];
  const protocols = data.protocols;
  const includeProperties = data.meta.include_properties || [];
  const propertyThresholds = data.property_thresholds || [];
  const showBatchColumn = data.meta.group_by_batch;

  // Create a map for quick threshold lookup
  const thresholdMap = useMemo(() => {
    const map: Record<string, MolecularPropertyThreshold> = {};
    for (const t of propertyThresholds) {
      map[t.property_name] = t;
    }
    return map;
  }, [propertyThresholds]);

  // Check which aggregations are selected for combined formatting
  const hasGeomean = aggregations.includes('geomean');
  const hasStdev = aggregations.includes('stdev');
  const hasCount = aggregations.includes('count');
  const hasList = aggregations.includes('list');

  // Handle sort request
  const handleRequestSort = (key: PivotSortKey) => {
    const isAsc = orderBy === key && order === 'asc';
    setOrder(isAsc ? 'desc' : 'asc');
    setOrderBy(key);
  };

  // Get sort value for a compound row
  const getSortValue = useCallback((row: CompactRow, key: PivotSortKey): unknown => {
    if (key === 'compound') return row.formatted_id;
    if (key.startsWith('property_')) {
      const propName = key.replace('property_', '') as MolecularPropertyName;
      return row.properties?.[propName] ?? null;
    }
    if (key.startsWith('protocol_')) {
      const protocolId = key.replace('protocol_', '');
      const protocolData = row.protocols[protocolId];
      // Sort by geomean if available, otherwise by count
      if (hasGeomean) return protocolData?.geomean ?? null;
      if (hasCount) return protocolData?.count ?? null;
      return null;
    }
    return null;
  }, [hasGeomean, hasCount]);

  // Sorted compounds (columns)
  const sortedRows = useMemo(() => {
    const comparator = getComparator<CompactRow>(order, (row) => getSortValue(row, orderBy));
    return [...rows].sort(comparator);
  }, [rows, order, orderBy, getSortValue]);

  // Handle click on protocol cell to show data series detail
  const handleProtocolCellClick = (
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

  // Build row definitions: properties first, then protocols
  // When geomean is selected, combine value ± stdev (n=count) in a single row
  type PivotRowDef =
    | { type: 'property'; name: MolecularPropertyName; label: string }
    | { type: 'protocol'; protocol: ProtocolInfo }
    | { type: 'protocol-list'; protocol: ProtocolInfo };

  const pivotRows = useMemo<PivotRowDef[]>(() => {
    const result: PivotRowDef[] = [];

    // Add property rows
    for (const propName of includeProperties) {
      result.push({
        type: 'property',
        name: propName,
        label: PROPERTY_LABELS[propName] || propName,
      });
    }

    // Add protocol rows - one row per protocol (combined format)
    for (const protocol of protocols) {
      // Main measurement row (geomean ± stdev (n=count) or just count)
      if (hasGeomean || hasCount) {
        result.push({ type: 'protocol', protocol });
      }
      // Separate list row if list is selected
      if (hasList) {
        result.push({ type: 'protocol-list', protocol });
      }
    }

    return result;
  }, [includeProperties, protocols, hasGeomean, hasCount, hasList]);

  return (
    <>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
          <Typography variant="body2" color="text.secondary">
            {data.meta.compound_count} compounds
            {data.meta.group_by_batch && data.meta.row_count && data.meta.row_count !== data.meta.compound_count && (
              <> ({data.meta.row_count} columns with batch split)</>
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
        <Tooltip title="Click row labels to sort compounds">
          <Chip label="Click to sort" size="small" variant="outlined" />
        </Tooltip>
      </Box>

      <Box
        ref={scrollContainerRef}
        sx={{
          overflow: 'auto',
          maxHeight: 600,
          border: 1,
          borderColor: 'divider',
          borderRadius: 1,
        }}
      >
        <Table size="small" sx={{ minWidth: 200 + sortedRows.length * 140 }}>
          <TableHead>
            <TableRow>
              {/* Fixed property column header */}
              <TableCell
                sx={{
                  fontWeight: 600,
                  position: 'sticky',
                  left: 0,
                  bgcolor: 'background.paper',
                  zIndex: 3,
                  width: 150,
                  borderRight: 1,
                  borderRightColor: 'divider',
                }}
              >
                <TableSortLabel
                  active={orderBy === 'compound'}
                  direction={orderBy === 'compound' ? order : 'asc'}
                  onClick={() => handleRequestSort('compound')}
                >
                  Property
                </TableSortLabel>
              </TableCell>
              {/* Compound column headers with structure */}
              {sortedRows.map((row) => (
                <TableCell
                  key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                  align="center"
                  sx={{
                    fontWeight: 600,
                    width: 130,
                    cursor: 'pointer',
                    '&:hover': { bgcolor: 'action.hover' },
                  }}
                  onClick={() => router.push(`/registry/compounds/${row.compound_id}`)}
                >
                  <Box sx={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: 0.5 }}>
                    {row.smiles ? (
                      <MoleculeChip smiles={row.smiles} size={90} />
                    ) : (
                      <Box
                        sx={{
                          width: 90,
                          height: 90,
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
                    <Typography variant="caption" fontWeight={600}>
                      {row.formatted_id}
                      {showBatchColumn && row.batch_number != null && `/${row.batch_number}`}
                    </Typography>
                  </Box>
                </TableCell>
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {pivotRows.map((rowDef) => {
              if (rowDef.type === 'property') {
                const threshold = thresholdMap[rowDef.name];
                const sortKey: PivotSortKey = `property_${rowDef.name}`;
                return (
                  <TableRow key={`prop-${rowDef.name}`}>
                    <TableCell
                      sx={{
                        position: 'sticky',
                        left: 0,
                        bgcolor: 'background.paper',
                        zIndex: 2,
                        fontWeight: 500,
                        borderRight: 1,
                        borderRightColor: 'divider',
                        cursor: 'pointer',
                      }}
                      onClick={() => handleRequestSort(sortKey)}
                    >
                      <Tooltip title={`Sort by ${threshold?.property_display || rowDef.name}`}>
                        <TableSortLabel
                          active={orderBy === sortKey}
                          direction={orderBy === sortKey ? order : 'asc'}
                        >
                          {rowDef.label}
                        </TableSortLabel>
                      </Tooltip>
                    </TableCell>
                    {sortedRows.map((row) => {
                      const value = row.properties?.[rowDef.name as keyof MolecularPropertyValues] as number | null | undefined;
                      const ragStatus = getRagStatus(value, threshold);
                      return (
                        <TableCell
                          key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                          align="center"
                        >
                          <Typography
                            variant="body2"
                            fontFamily="monospace"
                            sx={{ color: RAG_COLORS[ragStatus] }}
                          >
                            {formatPropertyValue(value)}
                          </Typography>
                        </TableCell>
                      );
                    })}
                  </TableRow>
                );
              }

              if (rowDef.type === 'protocol-list') {
                // Separate row for list values - not sortable
                const { protocol } = rowDef;
                return (
                  <TableRow key={`proto-list-${protocol.id}`}>
                    <TableCell
                      sx={{
                        position: 'sticky',
                        left: 0,
                        bgcolor: 'background.paper',
                        zIndex: 2,
                        borderRight: 1,
                        borderRightColor: 'divider',
                      }}
                    >
                      <Typography variant="body2">
                        {protocol.name}
                        <Typography component="span" variant="caption" color="text.secondary" sx={{ ml: 0.5 }}>
                          (values)
                        </Typography>
                      </Typography>
                    </TableCell>
                    {sortedRows.map((row) => {
                      const protocolData = row.protocols[protocol.id];
                      const measurementStatus = getMeasurementStatus(protocolData, hasCount);
                      const value = protocolData?.list;

                      // For not-measured, show not-tested indicator
                      if (measurementStatus === 'not-measured') {
                        return (
                          <TableCell
                            key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                            align="center"
                          >
                            <NotTestedIndicator />
                          </TableCell>
                        );
                      }

                      // For tested-no-valid in list row, show dash (badge shown in main row)
                      if (measurementStatus === 'tested-no-valid') {
                        return (
                          <TableCell
                            key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                            align="center"
                          >
                            <Typography variant="body2" color="text.disabled">-</Typography>
                          </TableCell>
                        );
                      }

                      return (
                        <TableCell
                          key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                          align="center"
                          sx={{ cursor: 'pointer', '&:hover': { bgcolor: 'action.hover' } }}
                          onClick={(e) => handleProtocolCellClick(e, row, protocol)}
                        >
                          <Tooltip title={value || '-'}>
                            <Typography
                              variant="body2"
                              fontFamily="monospace"
                              sx={{
                                maxWidth: 80,
                                overflow: 'hidden',
                                textOverflow: 'ellipsis',
                                whiteSpace: 'nowrap',
                              }}
                            >
                              {value || '-'}
                            </Typography>
                          </Tooltip>
                        </TableCell>
                      );
                    })}
                  </TableRow>
                );
              }

              // Protocol measurement row - combined format: value ± stdev (n=count)
              const { protocol } = rowDef;
              const sortKey: PivotSortKey = `protocol_${protocol.id}`;
              return (
                <TableRow key={`proto-${protocol.id}`}>
                  <TableCell
                    sx={{
                      position: 'sticky',
                      left: 0,
                      bgcolor: 'background.paper',
                      zIndex: 2,
                      fontWeight: 500,
                      borderRight: 1,
                      borderRightColor: 'divider',
                      cursor: 'pointer',
                    }}
                    onClick={() => handleRequestSort(sortKey)}
                  >
                    <Tooltip title={`Sort by ${protocol.name}`}>
                      <TableSortLabel
                        active={orderBy === sortKey}
                        direction={orderBy === sortKey ? order : 'asc'}
                      >
                        {protocol.name}
                      </TableSortLabel>
                    </Tooltip>
                    {protocol.kpi_unit && hasGeomean && (
                      <Typography component="span" variant="caption" color="text.secondary" sx={{ ml: 0.5 }}>
                        ({getConcentrationHeaderUnit(protocol.kpi_unit, concentrationDisplay)})
                      </Typography>
                    )}
                  </TableCell>
                  {sortedRows.map((row) => {
                    const protocolData = row.protocols[protocol.id];
                    const measurementStatus = getMeasurementStatus(protocolData, hasCount);

                    // Handle not measured - compound was never tested with this protocol
                    if (measurementStatus === 'not-measured') {
                      return (
                        <TableCell
                          key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                          align="center"
                        >
                          <NotTestedIndicator />
                        </TableCell>
                      );
                    }

                    // Handle tested but no valid data
                    if (measurementStatus === 'tested-no-valid') {
                      return (
                        <TableCell
                          key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                          align="center"
                        >
                          <TestedNoValidBadge
                            tested={protocolData?.tested}
                            noAnalysis={protocolData?.no_analysis}
                            invalid={protocolData?.invalid}
                            unassigned={protocolData?.unassigned}
                            onClick={(e) => handleProtocolCellClick(e, row, protocol)}
                          />
                        </TableCell>
                      );
                    }

                    // Use combined format when geomean is selected
                    let displayValue: string;
                    if (hasGeomean) {
                      displayValue = formatMeasurementWithStats(
                        protocolData?.geomean,
                        protocolData?.stdev,
                        protocolData?.stdev_log,
                        protocolData?.count,
                        protocol.kpi_unit,
                        concentrationDisplay,
                        hasGeomean,
                        hasStdev,
                        hasCount
                      );
                    } else if (hasCount) {
                      // Only count selected
                      displayValue = protocolData?.count != null ? `n=${protocolData.count}` : '-';
                    } else {
                      displayValue = '-';
                    }

                    return (
                      <TableCell
                        key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                        align="center"
                        sx={{ cursor: 'pointer', '&:hover': { bgcolor: 'action.hover' } }}
                        onClick={(e) => handleProtocolCellClick(e, row, protocol)}
                      >
                        <Typography variant="body2" fontFamily="monospace">
                          {displayValue}
                        </Typography>
                      </TableCell>
                    );
                  })}
                </TableRow>
              );
            })}
          </TableBody>
        </Table>
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

/** Sort key type for cards view */
type CardsSortKey = 'compound' | `property_${MolecularPropertyName}` | `protocol_${string}`;

/**
 * Cards view component (grid of compound cards).
 * Each card shows structure, compound info, and all protocols with formatted measurements.
 */
function CardsView({
  data,
  aggregations,
  concentrationDisplay = 'natural',
}: {
  data: CompactAggregationResponse;
  aggregations: AggregationType[];
  concentrationDisplay: ConcentrationDisplayMode;
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

  const rows = data.data as CompactRow[];
  const protocols = data.protocols;
  const includeProperties = data.meta.include_properties || [];
  const propertyThresholds = data.property_thresholds || [];
  const showBatch = data.meta.group_by_batch;

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
        {/* Sort controls */}
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
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
          gridTemplateColumns: 'repeat(auto-fill, minmax(300px, 1fr))',
          gap: 2,
          maxHeight: 600,
          overflow: 'auto',
          p: 1,
        }}
      >
        {sortedRows.map((row) => (
          <Paper
            key={showBatch ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
            variant="outlined"
            sx={{
              p: 2,
              cursor: 'pointer',
              '&:hover': { boxShadow: 2, borderColor: 'primary.main' },
              display: 'flex',
              flexDirection: 'column',
            }}
            onClick={() => router.push(`/registry/compounds/${row.compound_id}`)}
          >
            {/* Header: Structure + Compound ID */}
            <Box sx={{ display: 'flex', gap: 2, mb: 1.5 }}>
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
                    flexShrink: 0,
                  }}
                >
                  <Typography variant="caption" color="text.secondary">-</Typography>
                </Box>
              )}
              <Box sx={{ flex: 1, minWidth: 0 }}>
                <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                  <Chip
                    icon={<Medication fontSize="small" />}
                    label={row.formatted_id}
                    size="small"
                    variant="outlined"
                    color="primary"
                  />
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

            {/* Protocols */}
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
          </Paper>
        ))}
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
  outputFormat,
  concentrationDisplay = 'natural',
  onConcentrationDisplayChange,
  fillHeight = false,
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
    // For pivot and cards formats, we need compact response data
    if (outputFormat === 'pivot' && isCompactResponse(data)) {
      return (
        <PivotTable
          data={data}
          aggregations={aggregations}
          concentrationDisplay={displayMode}
        />
      );
    }
    if (outputFormat === 'cards' && isCompactResponse(data)) {
      return (
        <CardsView
          data={data}
          aggregations={aggregations}
          concentrationDisplay={displayMode}
        />
      );
    }
    // Default: auto-detect from response type
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
    <Paper sx={{ width: '100%', overflow: 'hidden', ...(fillHeight && { height: '100%', display: 'flex', flexDirection: 'column' }) }}>
      <Box sx={{ p: 2, ...(fillHeight && { flex: 1, display: 'flex', flexDirection: 'column', minHeight: 0, overflow: 'hidden' }) }}>
        <Box sx={{ display: 'flex', justifyContent: 'flex-end', mb: 1, flexShrink: 0 }}>
          <ConcentrationDisplaySelector
            value={displayMode}
            onChange={handleDisplayChange}
          />
        </Box>
        <Box sx={{ ...(fillHeight && { flex: 1, minHeight: 0, overflow: 'auto' }) }}>
          {renderTable()}
        </Box>
      </Box>
    </Paper>
  );
}
