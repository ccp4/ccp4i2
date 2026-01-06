'use client';

import { useState, useMemo } from 'react';
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
  TablePagination,
  Button,
  Chip,
  LinearProgress,
  Tooltip,
} from '@mui/material';
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
} from '@/types/aggregation';
import { MoleculeChip } from './MoleculeView';
import { DataSeriesDetailModal } from './DataSeriesDetailModal';
import {
  formatKpiValue,
  generateCompactCsv,
  generateLongCsv,
  downloadCsv,
} from '@/lib/aggregation-api';

/**
 * Generate CSV content from medium aggregation data.
 */
function generateMediumCsv(
  data: AggregationResponse,
  aggregations: AggregationType[]
): string {
  const { data: rows } = data;

  // Build header row
  const headers = ['Compound ID', 'SMILES', 'Target', 'Protocol'];
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

interface AggregationTableProps {
  data: AggregationResponse | null | undefined;
  loading?: boolean;
  aggregations: AggregationType[];
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
 */
function CompactTable({
  data,
  aggregations,
}: {
  data: AggregationResponse & { protocols: ProtocolInfo[] };
  aggregations: AggregationType[];
}) {
  const router = useRouter();
  const [page, setPage] = useState(0);
  const [rowsPerPage, setRowsPerPage] = useState(25);

  const rows = data.data as CompactRow[];
  const protocols = data.protocols;

  const paginatedRows = useMemo(() => {
    const start = page * rowsPerPage;
    return rows.slice(start, start + rowsPerPage);
  }, [rows, page, rowsPerPage]);

  const handleExport = () => {
    const csv = generateCompactCsv(data, aggregations);
    downloadCsv(csv, `aggregation_compact_${new Date().toISOString().slice(0, 10)}.csv`);
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

      <TableContainer sx={{ maxHeight: 600 }}>
        <Table stickyHeader size="small">
          <TableHead>
            <TableRow>
              <TableCell sx={{ fontWeight: 600, minWidth: 80 }}>Structure</TableCell>
              <TableCell sx={{ fontWeight: 600, minWidth: 120 }}>Compound</TableCell>
              <TableCell sx={{ fontWeight: 600, minWidth: 100 }}>Target</TableCell>
              {protocols.map((protocol) => (
                aggregations.map((agg) => (
                  <TableCell
                    key={`${protocol.id}-${agg}`}
                    sx={{ fontWeight: 600, minWidth: 80 }}
                    align="right"
                  >
                    <Tooltip title={protocol.name}>
                      <span>
                        {protocol.name.length > 15
                          ? `${protocol.name.slice(0, 15)}...`
                          : protocol.name}
                        <br />
                        <Typography variant="caption" color="text.secondary">
                          {formatAggLabel(agg)}
                        </Typography>
                      </span>
                    </Tooltip>
                  </TableCell>
                ))
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {paginatedRows.map((row) => (
              <TableRow
                key={row.compound_id}
                hover
                sx={{ cursor: 'pointer' }}
                onClick={() => router.push(`/registry/compounds/${row.compound_id}`)}
              >
                <TableCell>
                  {row.smiles ? (
                    <MoleculeChip smiles={row.smiles} size={60} />
                  ) : (
                    <Box
                      sx={{
                        width: 60,
                        height: 60,
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
                <TableCell>
                  <Typography variant="body2">{row.target_name || '-'}</Typography>
                </TableCell>
                {protocols.map((protocol) => (
                  aggregations.map((agg) => {
                    const protocolData = row.protocols[protocol.id];
                    const value = protocolData?.[agg];

                    return (
                      <TableCell key={`${protocol.id}-${agg}`} align="right">
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
                        ) : (
                          <Typography variant="body2" fontFamily="monospace">
                            {agg === 'count'
                              ? value ?? '-'
                              : formatKpiValue(value as number | null)}
                          </Typography>
                        )}
                      </TableCell>
                    );
                  })
                ))}
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </TableContainer>

      <TablePagination
        component="div"
        count={rows.length}
        page={page}
        onPageChange={(_, newPage) => setPage(newPage)}
        rowsPerPage={rowsPerPage}
        onRowsPerPageChange={(e) => {
          setRowsPerPage(parseInt(e.target.value, 10));
          setPage(0);
        }}
        rowsPerPageOptions={[10, 25, 50, 100]}
      />
    </>
  );
}

/**
 * Medium table component (one row per compound-protocol pair).
 * Clicking a row opens a modal with data series details and charts.
 */
function MediumTable({
  data,
  aggregations,
}: {
  data: AggregationResponse;
  aggregations: AggregationType[];
}) {
  const [page, setPage] = useState(0);
  const [rowsPerPage, setRowsPerPage] = useState(25);

  // Modal state
  const [modalOpen, setModalOpen] = useState(false);
  const [selectedRow, setSelectedRow] = useState<MediumRow | null>(null);

  const rows = data.data as MediumRow[];

  const paginatedRows = useMemo(() => {
    const start = page * rowsPerPage;
    return rows.slice(start, start + rowsPerPage);
  }, [rows, page, rowsPerPage]);

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
            {data.meta.compound_count} compounds, {data.meta.protocol_count} protocols,{' '}
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

      <TableContainer sx={{ maxHeight: 600 }}>
        <Table stickyHeader size="small">
          <TableHead>
            <TableRow>
              <TableCell sx={{ fontWeight: 600, minWidth: 80 }}>Structure</TableCell>
              <TableCell sx={{ fontWeight: 600, minWidth: 120 }}>Compound</TableCell>
              <TableCell sx={{ fontWeight: 600, minWidth: 100 }}>Target</TableCell>
              <TableCell sx={{ fontWeight: 600, minWidth: 120 }}>Protocol</TableCell>
              {aggregations.map((agg) => (
                <TableCell key={agg} sx={{ fontWeight: 600 }} align="right">
                  {formatAggLabel(agg)}
                </TableCell>
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {paginatedRows.map((row, idx) => (
              <TableRow
                key={`${row.compound_id}-${row.protocol_id}-${idx}`}
                hover
                sx={{ cursor: 'pointer' }}
                onClick={() => handleRowClick(row)}
              >
                <TableCell>
                  {row.smiles ? (
                    <MoleculeChip smiles={row.smiles} size={60} />
                  ) : (
                    <Box
                      sx={{
                        width: 60,
                        height: 60,
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
                <TableCell>
                  <Typography variant="body2">{row.target_name || '-'}</Typography>
                </TableCell>
                <TableCell>
                  <Typography variant="body2">{row.protocol_name}</Typography>
                </TableCell>
                {aggregations.map((agg) => {
                  const value = row[agg as keyof MediumRow];
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
                      ) : (
                        <Typography variant="body2" fontFamily="monospace">
                          {agg === 'count'
                            ? value ?? '-'
                            : formatKpiValue(value as number | null)}
                        </Typography>
                      )}
                    </TableCell>
                  );
                })}
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </TableContainer>

      <TablePagination
        component="div"
        count={rows.length}
        page={page}
        onPageChange={(_, newPage) => setPage(newPage)}
        rowsPerPage={rowsPerPage}
        onRowsPerPageChange={(e) => {
          setRowsPerPage(parseInt(e.target.value, 10));
          setPage(0);
        }}
        rowsPerPageOptions={[10, 25, 50, 100]}
      />

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

/**
 * Long table component (one row per measurement).
 */
function LongTable({ data }: { data: AggregationResponse }) {
  const router = useRouter();
  const [page, setPage] = useState(0);
  const [rowsPerPage, setRowsPerPage] = useState(25);

  const rows = data.data as LongRow[];

  const paginatedRows = useMemo(() => {
    const start = page * rowsPerPage;
    return rows.slice(start, start + rowsPerPage);
  }, [rows, page, rowsPerPage]);

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

      <TableContainer sx={{ maxHeight: 600 }}>
        <Table stickyHeader size="small">
          <TableHead>
            <TableRow>
              <TableCell sx={{ fontWeight: 600, minWidth: 80 }}>Structure</TableCell>
              <TableCell sx={{ fontWeight: 600, minWidth: 120 }}>Compound</TableCell>
              <TableCell sx={{ fontWeight: 600, minWidth: 100 }}>Target</TableCell>
              <TableCell sx={{ fontWeight: 600, minWidth: 120 }}>Protocol</TableCell>
              <TableCell sx={{ fontWeight: 600 }}>Date</TableCell>
              <TableCell sx={{ fontWeight: 600 }} align="right">KPI</TableCell>
              <TableCell sx={{ fontWeight: 600 }}>Status</TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {paginatedRows.map((row) => (
              <TableRow
                key={row.data_series_id}
                hover
                sx={{ cursor: 'pointer' }}
                onClick={() => row.compound_id && router.push(`/registry/compounds/${row.compound_id}`)}
              >
                <TableCell>
                  {row.smiles ? (
                    <MoleculeChip smiles={row.smiles} size={60} />
                  ) : (
                    <Box
                      sx={{
                        width: 60,
                        height: 60,
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
                  <Typography variant="body2" fontFamily="monospace" fontWeight={500}>
                    {formatKpiValue(row.kpi_value)}
                  </Typography>
                </TableCell>
                <TableCell>
                  <Chip
                    label={row.status || 'unknown'}
                    size="small"
                    color={row.status === 'valid' ? 'success' : row.status === 'invalid' ? 'error' : 'default'}
                  />
                </TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </TableContainer>

      <TablePagination
        component="div"
        count={rows.length}
        page={page}
        onPageChange={(_, newPage) => setPage(newPage)}
        rowsPerPage={rowsPerPage}
        onRowsPerPageChange={(e) => {
          setRowsPerPage(parseInt(e.target.value, 10));
          setPage(0);
        }}
        rowsPerPageOptions={[10, 25, 50, 100]}
      />
    </>
  );
}

/**
 * Main aggregation table component.
 * Renders either compact or long format based on the response data.
 */
export function AggregationTable({ data, loading, aggregations }: AggregationTableProps) {
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
          Configure predicates above and click Generate Report to see results.
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
      return <CompactTable data={data} aggregations={aggregations} />;
    } else if (isMediumResponse(data)) {
      return <MediumTable data={data} aggregations={aggregations} />;
    } else {
      return <LongTable data={data} />;
    }
  };

  return (
    <Paper sx={{ width: '100%', overflow: 'hidden' }}>
      <Box sx={{ p: 2 }}>
        {renderTable()}
      </Box>
    </Paper>
  );
}
