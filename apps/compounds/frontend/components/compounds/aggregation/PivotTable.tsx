'use client';

import { useState, useRef, useCallback, useMemo } from 'react';
import {
  Box,
  Typography,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  TableSortLabel,
  Chip,
  Tooltip,
  IconButton,
} from '@mui/material';
import { ZoomIn, EditOutlined, BubbleChart } from '@mui/icons-material';
import { useRouter } from 'next/navigation';
import {
  AggregationType,
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
import {
  getConcentrationHeaderUnit,
} from '@/lib/compounds/aggregation-api';
import { protocolColour } from '@/lib/compounds/protocol-colour';
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

/** Sort key for pivot table - sort compounds (columns) by property or protocol */
type PivotSortKey = 'compound' | `property_${MolecularPropertyName}` | `protocol_${string}`;

/**
 * Pivot table component (compounds as columns, properties/protocols as rows).
 * Horizontally scrollable for many compounds. Click row labels to sort columns.
 */
export function PivotTable({
  data,
  aggregations,
  concentrationDisplay = 'natural',
  onEditProtocol,
  onScatterProtocol,
}: {
  data: CompactAggregationResponse;
  aggregations: AggregationType[];
  concentrationDisplay: ConcentrationDisplayMode;
  onEditProtocol?: (protocol: ProtocolInfo) => void;
  onScatterProtocol?: (protocol: ProtocolInfo) => void;
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
  const propertyThresholds = data.property_thresholds ?? EMPTY_THRESHOLDS;
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

  // Build row definitions: identifiers (optional), properties, then protocols
  // When geomean is selected, combine value ± stdev (n=count) in a single row
  type PivotRowDef =
    | { type: 'identifiers' }
    | { type: 'property'; name: MolecularPropertyName; label: string }
    | { type: 'protocol'; protocol: ProtocolInfo }
    | { type: 'protocol-list'; protocol: ProtocolInfo };

  const showIdentifiersRow = Boolean(data.meta.include_identifiers);

  const pivotRows = useMemo<PivotRowDef[]>(() => {
    const result: PivotRowDef[] = [];

    if (showIdentifiersRow) {
      result.push({ type: 'identifiers' });
    }

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
  }, [showIdentifiersRow, includeProperties, protocols, hasGeomean, hasCount, hasList]);

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
                    <Tooltip title={row.formatted_id + (showBatchColumn && row.batch_number != null ? `/${row.batch_number}` : '')}>
                      <Typography variant="caption" fontWeight={600} sx={{ fontFamily: 'monospace' }}>
                        {row.formatted_id}
                        {showBatchColumn && row.batch_number != null && `/${row.batch_number}`}
                      </Typography>
                    </Tooltip>
                  </Box>
                </TableCell>
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {pivotRows.map((rowDef) => {
              if (rowDef.type === 'identifiers') {
                return (
                  <TableRow key="identifiers">
                    <TableCell
                      sx={{
                        position: 'sticky',
                        left: 0,
                        bgcolor: 'background.paper',
                        zIndex: 2,
                        fontWeight: 500,
                        borderRight: 1,
                        borderRightColor: 'divider',
                      }}
                    >
                      <Typography variant="body2">Identifiers</Typography>
                    </TableCell>
                    {sortedRows.map((row) => (
                      <TableCell
                        key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                        align="center"
                        sx={{ verticalAlign: 'top' }}
                      >
                        <IdentifiersCell identifiers={row.identifiers} />
                      </TableCell>
                    ))}
                  </TableRow>
                );
              }
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
                    {onEditProtocol && (
                      <Tooltip title="Edit thresholds for this protocol" arrow>
                        <IconButton
                          size="small"
                          onClick={(e) => {
                            e.stopPropagation();
                            onEditProtocol(protocol);
                          }}
                          sx={{ ml: 0.5, opacity: 0.5, '&:hover': { opacity: 1 } }}
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
                          sx={{ opacity: 0.5, '&:hover': { opacity: 1 } }}
                        >
                          <BubbleChart fontSize="inherit" />
                        </IconButton>
                      </Tooltip>
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

                    // Apply threshold-based background colour when geomean is shown
                    // and the protocol has curated thresholds
                    const colour = hasGeomean
                      ? protocolColour(protocolData?.geomean, protocol)
                      : { background: null, t: null };

                    return (
                      <TableCell
                        key={showBatchColumn ? `${row.compound_id}-${row.batch_id}` : row.compound_id}
                        align="center"
                        sx={{
                          cursor: 'pointer',
                          bgcolor: colour.background ?? undefined,
                          '&:hover': colour.background
                            ? { filter: 'brightness(0.95)' }
                            : { bgcolor: 'action.hover' },
                        }}
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
