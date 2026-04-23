'use client';

import { useState, useRef, useCallback, useMemo } from 'react';
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
} from '@mui/material';
import { useVirtualizer } from '@tanstack/react-virtual';
import { Download, ZoomIn } from '@mui/icons-material';
import {
  AggregationResponse,
  AggregationType,
  MediumRow,
  ConcentrationDisplayMode,
} from '@/types/compounds/aggregation';
import { MoleculeChip, CompoundNameChip } from '../MoleculeView';
import { DataSeriesDetailModal } from '../DataSeriesDetailModal';
import {
  formatKpiUnit,
  downloadCsv,
  formatConcentrationValue,
} from '@/lib/compounds/aggregation-api';
import {
  Order,
  IdentifiersCell,
  formatAggLabel,
  getComparator,
  TestedNoValidBadge,
} from './shared';

/**
 * Generate CSV content from medium aggregation data.
 */
function generateMediumCsv(
  data: AggregationResponse,
  aggregations: AggregationType[]
): string {
  const { data: rows } = data;
  const includeIdentifiers = Boolean(data.meta.include_identifiers);

  // Build header row (include Unit column)
  const headers = [
    'Compound ID',
    'SMILES',
    'Target',
    ...(includeIdentifiers ? ['Barcode', 'Supplier Ref', 'Aliases'] : []),
    'Protocol',
    'Unit',
  ];
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
    ];

    if (includeIdentifiers) {
      const ids = row.identifiers;
      const aliases = ids && Array.isArray(ids.aliases) ? ids.aliases.join('; ') : '';
      values.push(
        `"${ids?.barcode || ''}"`,
        `"${ids?.supplier_ref || ''}"`,
        `"${aliases}"`,
      );
    }

    values.push(
      `"${row.protocol_name}"`,
      `"${row.kpi_unit ? formatKpiUnit(row.kpi_unit) : ''}"`,
    );

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

/** Sort key for medium table */
type MediumSortKey = 'compound' | 'batch' | 'target' | 'protocol' | AggregationType;

/**
 * Medium table component (one row per compound-protocol pair).
 * Clicking a row opens a modal with data series details and charts.
 */
export function MediumTable({
  data,
  aggregations,
  concentrationDisplay = 'natural',
  searchTerm = '',
  fillHeight = false,
}: {
  data: AggregationResponse;
  aggregations: AggregationType[];
  concentrationDisplay: ConcentrationDisplayMode;
  searchTerm?: string;
  fillHeight?: boolean;
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
  const showIdentifiersColumn = Boolean(data.meta.include_identifiers);

  // Filter rows by search term
  const filteredRows = useMemo(() => {
    if (!searchTerm) return rows;
    const term = searchTerm.toLowerCase();
    return rows.filter(row => row.formatted_id?.toLowerCase().includes(term));
  }, [rows, searchTerm]);

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

  // Sorted rows (from filtered set)
  const sortedRows = useMemo(() => {
    const comparator = getComparator<MediumRow>(order, (row) => getSortValue(row, orderBy));
    return [...filteredRows].sort(comparator);
  }, [filteredRows, order, orderBy, getSortValue]);

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
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2, flexShrink: 0 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
          <Typography variant="body2" color="text.secondary">
            {searchTerm ? `${filteredRows.length} of ` : ''}{data.meta.compound_count} compounds
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
        sx={fillHeight ? { flex: 1, minHeight: 0, overflow: 'auto' } : { maxHeight: 600, overflow: 'auto' }}
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
                  colSpan={4 + (showBatchColumn ? 1 : 0) + (showIdentifiersColumn ? 1 : 0) + aggregations.length}
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
