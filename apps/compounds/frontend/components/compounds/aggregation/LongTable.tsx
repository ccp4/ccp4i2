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
import { Download } from '@mui/icons-material';
import { useRouter } from 'next/navigation';
import {
  AggregationResponse,
  LongRow,
  ConcentrationDisplayMode,
} from '@/types/compounds/aggregation';
import { MoleculeChip, CompoundNameChip } from '../MoleculeView';
import {
  generateLongCsv,
  downloadCsv,
  formatConcentrationValue,
} from '@/lib/compounds/aggregation-api';
import {
  Order,
  IdentifiersCell,
  getComparator,
} from './shared';

/** Sort key for long table */
type LongSortKey = 'compound' | 'batch' | 'target' | 'protocol' | 'date' | 'kpi' | 'status';

/**
 * Long table component (one row per measurement).
 */
export function LongTable({
  data,
  concentrationDisplay = 'natural',
  searchTerm = '',
  fillHeight = false,
}: {
  data: AggregationResponse;
  concentrationDisplay: ConcentrationDisplayMode;
  searchTerm?: string;
  fillHeight?: boolean;
}) {
  const router = useRouter();
  const parentRef = useRef<HTMLDivElement>(null);

  // Sorting state
  const [orderBy, setOrderBy] = useState<LongSortKey>('compound');
  const [order, setOrder] = useState<Order>('asc');

  const rows = data.data as LongRow[];
  // Long format always includes batch info, but only show column if any row has batch data
  const showBatchColumn = rows.some(row => row.batch_number != null);
  const showIdentifiersColumn = Boolean(data.meta.include_identifiers);

  // Filter rows by search term
  const filteredRows = useMemo(() => {
    if (!searchTerm) return rows;
    const term = searchTerm.toLowerCase();
    return rows.filter(row =>
      row.formatted_id?.toLowerCase().includes(term) ||
      row.compound_name?.toLowerCase().includes(term)
    );
  }, [rows, searchTerm]);

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

  // Sorted rows (from filtered set)
  const sortedRows = useMemo(() => {
    const comparator = getComparator<LongRow>(order, (row) => getSortValue(row, orderBy));
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
    const csv = generateLongCsv(data);
    downloadCsv(csv, `aggregation_long_${new Date().toISOString().slice(0, 10)}.csv`);
  };

  return (
    <>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2, flexShrink: 0 }}>
        <Typography variant="body2" color="text.secondary">
          {searchTerm ? `${filteredRows.length} of ` : ''}{data.meta.compound_count} compounds, {data.meta.protocol_count} protocols,{' '}
          {data.meta.total_measurements} measurements
        </Typography>
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
                  colSpan={7 + (showBatchColumn ? 1 : 0) + (showIdentifiersColumn ? 1 : 0)}
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
                      <CompoundNameChip formattedId={row.formatted_id} smiles={row.smiles} label={String(Number(row.formatted_id.split('-').pop()))} />
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
