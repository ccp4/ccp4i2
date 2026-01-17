'use client';

import { useState, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Box,
  Typography,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Chip,
  CircularProgress,
  IconButton,
  Alert,
  Paper,
} from '@mui/material';
import { Close, Medication, Science, OpenInNew } from '@mui/icons-material';
import Link from 'next/link';
import { MoleculeChip } from './MoleculeView';
import { DoseResponseChart, DoseResponseThumb, FitParameters, DoseResponseData } from './DoseResponseChart';
import { formatKpiValue } from '@/lib/compounds/aggregation-api';

// =============================================================================
// Authentication Integration
// =============================================================================

// For Docker integration: Try to import auth helpers from ccp4i2 client's auth-token
// Falls back to no-op for standalone development
let getAccessToken: () => Promise<string | null>;
let getUserEmail: () => string | null;

try {
  const authModule = require('../../utils/auth-token');
  getAccessToken = authModule.getAccessToken;
  getUserEmail = authModule.getUserEmail || (() => null);
} catch {
  getAccessToken = async () => null;
  getUserEmail = () => null;
}

/**
 * Authenticated fetch wrapper
 */
async function authFetch(url: string, options: RequestInit = {}): Promise<Response> {
  const headers: Record<string, string> = {};

  const token = await getAccessToken();
  if (token) {
    headers['Authorization'] = `Bearer ${token}`;
  }

  const email = getUserEmail();
  if (email) {
    headers['X-User-Email'] = email;
  }

  if (options.headers) {
    if (options.headers instanceof Headers) {
      options.headers.forEach((value, key) => {
        headers[key] = value;
      });
    } else if (Array.isArray(options.headers)) {
      options.headers.forEach(([key, value]) => {
        headers[key] = value;
      });
    } else {
      Object.assign(headers, options.headers);
    }
  }

  return fetch(url, { ...options, headers });
}

interface DataSeriesItem {
  id: string;
  assay: string;
  compound: string | null;
  compound_formatted_id: string | null;
  compound_name: string | null;
  row: number | null;
  start_column: number | null;
  end_column: number | null;
  dilution_series: {
    id: string;
    concentrations: number[];
    unit: string;
    display_name: string;
  } | null;
  extracted_data: number[] | null;
  analysis: {
    id: string;
    status: string;
    results: Record<string, any> | null;
    kpi_value: number | null;
  } | null;
  analysis_status: string | null;
  analysis_kpi: number | null;
}

interface DataSeriesResponse {
  compound: {
    id: string;
    formatted_id: string;
    smiles: string | null;
  } | null;
  protocol: {
    id: string;
    name: string;
  } | null;
  count: number;
  data_series: DataSeriesItem[];
}

interface DataSeriesDetailModalProps {
  open: boolean;
  onClose: () => void;
  compoundId: string;
  protocolId: string;
  compoundName?: string;
  protocolName?: string;
}

/**
 * Extract fit parameters from analysis results for DoseResponseChart
 */
function extractFitParams(analysis: DataSeriesItem['analysis']): FitParameters | undefined {
  if (!analysis?.results) return undefined;

  const results = analysis.results;

  return {
    ec50: results.EC50 ?? results.ec50 ?? null,
    hill: results.Hill ?? results.hill ?? results.slope ?? null,
    minVal: results.minVal ?? results.Min ?? results.min ?? results.bottom ?? null,
    maxVal: results.maxVal ?? results.Max ?? results.max ?? results.top ?? null,
    status: analysis.status,
  };
}

/**
 * Extract dose-response data from a data series - requires dilution_series
 */
function extractDoseResponseData(series: DataSeriesItem): DoseResponseData | null {
  if (!series.dilution_series?.concentrations || !series.extracted_data) {
    return null;
  }

  const concentrations = series.dilution_series.concentrations;
  let responses = series.extracted_data;
  const unit = series.dilution_series.unit || 'nM';

  // Detect format: if extracted_data has 2 more elements than concentrations,
  // it has embedded controls at first and last positions
  if (responses.length === concentrations.length + 2) {
    responses = responses.slice(1, -1);
  }

  return { concentrations, responses, unit };
}

export function DataSeriesDetailModal({
  open,
  onClose,
  compoundId,
  protocolId,
  compoundName,
  protocolName,
}: DataSeriesDetailModalProps) {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [data, setData] = useState<DataSeriesResponse | null>(null);
  const [selectedSeries, setSelectedSeries] = useState<DataSeriesItem | null>(null);

  useEffect(() => {
    if (!open || !compoundId || !protocolId) {
      return;
    }

    const fetchData = async () => {
      setLoading(true);
      setError(null);
      setSelectedSeries(null);

      try {
        const params = new URLSearchParams({
          compound: compoundId,
          protocol: protocolId,
        });

        const response = await authFetch(`/api/proxy/compounds/aggregations/data_series/?${params}`);

        if (!response.ok) {
          throw new Error(`Failed to fetch data series: ${response.statusText}`);
        }

        const result = await response.json();
        setData(result);

        // Auto-select the first series for chart preview
        if (result.data_series?.length > 0) {
          setSelectedSeries(result.data_series[0]);
        }
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Failed to fetch data');
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, [open, compoundId, protocolId]);

  const displayCompoundName = data?.compound?.formatted_id || compoundName || 'Unknown Compound';
  const displayProtocolName = data?.protocol?.name || protocolName || 'Unknown Protocol';

  return (
    <Dialog
      open={open}
      onClose={onClose}
      maxWidth="lg"
      fullWidth
      PaperProps={{ sx: { minHeight: '70vh' } }}
    >
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
        <Box sx={{ flex: 1, display: 'flex', alignItems: 'center', gap: 2 }}>
          <Medication color="secondary" />
          <Typography variant="h6" component="span">
            {displayCompoundName}
          </Typography>
          <Typography variant="body1" color="text.secondary">
            &mdash;
          </Typography>
          <Science color="primary" />
          <Typography variant="h6" component="span">
            {displayProtocolName}
          </Typography>
        </Box>
        <IconButton onClick={onClose} size="small">
          <Close />
        </IconButton>
      </DialogTitle>

      <DialogContent dividers>
        {loading && (
          <Box sx={{ display: 'flex', justifyContent: 'center', py: 6 }}>
            <CircularProgress />
          </Box>
        )}

        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        {!loading && !error && data && (
          <Box sx={{ display: 'flex', gap: 3, flexDirection: { xs: 'column', md: 'row' } }}>
            {/* Left: Data series table */}
            <Box sx={{ flex: 1, minWidth: 0 }}>
              <Typography variant="subtitle2" gutterBottom>
                {data.count} Data Series
              </Typography>

              {data.compound?.smiles && (
                <Box sx={{ mb: 2 }}>
                  <MoleculeChip smiles={data.compound.smiles} size={80} />
                </Box>
              )}

              <TableContainer component={Paper} variant="outlined" sx={{ maxHeight: 400 }}>
                <Table size="small" stickyHeader>
                  <TableHead>
                    <TableRow>
                      <TableCell sx={{ fontWeight: 600, width: 100 }}>Chart</TableCell>
                      <TableCell sx={{ fontWeight: 600 }}>KPI</TableCell>
                      <TableCell sx={{ fontWeight: 600 }}>Status</TableCell>
                      <TableCell sx={{ fontWeight: 600 }}>Assay</TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {data.data_series.map((series) => {
                      const doseData = extractDoseResponseData(series);
                      const fitParams = extractFitParams(series.analysis);
                      const isSelected = selectedSeries?.id === series.id;

                      return (
                        <TableRow
                          key={series.id}
                          hover
                          selected={isSelected}
                          onClick={() => setSelectedSeries(series)}
                          sx={{ cursor: 'pointer' }}
                        >
                          <TableCell>
                            {doseData ? (
                              <DoseResponseThumb data={doseData} fit={fitParams} size={80} />
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
                                <Typography variant="caption" color="text.secondary">
                                  No data
                                </Typography>
                              </Box>
                            )}
                          </TableCell>
                          <TableCell>
                            <Typography variant="body2" fontFamily="monospace" fontWeight={500}>
                              {formatKpiValue(series.analysis_kpi)}
                            </Typography>
                          </TableCell>
                          <TableCell>
                            <Chip
                              label={series.analysis_status || 'unknown'}
                              size="small"
                              color={
                                series.analysis_status === 'valid'
                                  ? 'success'
                                  : series.analysis_status === 'invalid'
                                  ? 'error'
                                  : 'default'
                              }
                            />
                          </TableCell>
                          <TableCell>
                            <Link
                              href={`/assays/${series.assay}`}
                              onClick={(e) => e.stopPropagation()}
                            >
                              <Typography
                                variant="caption"
                                color="primary"
                                sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}
                              >
                                View Assay <OpenInNew fontSize="inherit" />
                              </Typography>
                            </Link>
                          </TableCell>
                        </TableRow>
                      );
                    })}
                  </TableBody>
                </Table>
              </TableContainer>

              {data.data_series.length === 0 && (
                <Typography color="text.secondary" sx={{ mt: 2, textAlign: 'center' }}>
                  No data series found for this compound-protocol pair.
                </Typography>
              )}
            </Box>

            {/* Right: Selected series chart */}
            <Box sx={{ width: { xs: '100%', md: 450 }, flexShrink: 0 }}>
              <Typography variant="subtitle2" gutterBottom>
                Dose-Response Curve
              </Typography>

              {selectedSeries ? (
                (() => {
                  const doseData = extractDoseResponseData(selectedSeries);
                  const fitParams = extractFitParams(selectedSeries.analysis);

                  if (!doseData) {
                    return (
                      <Paper
                        variant="outlined"
                        sx={{
                          height: 350,
                          display: 'flex',
                          alignItems: 'center',
                          justifyContent: 'center',
                        }}
                      >
                        <Typography color="text.secondary">
                          No dose-response data available
                        </Typography>
                      </Paper>
                    );
                  }

                  return (
                    <Box>
                      <DoseResponseChart
                        data={doseData}
                        fit={fitParams}
                        title={selectedSeries.compound_name || displayCompoundName}
                        width={430}
                        height={350}
                        showLegend={true}
                      />

                      {/* Fit parameters */}
                      {fitParams && (
                        <Paper variant="outlined" sx={{ p: 2, mt: 2 }}>
                          <Typography variant="subtitle2" gutterBottom>
                            Fit Parameters
                          </Typography>
                          <Box
                            sx={{
                              display: 'grid',
                              gridTemplateColumns: 'repeat(2, 1fr)',
                              gap: 1,
                            }}
                          >
                            <Box>
                              <Typography variant="caption" color="text.secondary">
                                EC50
                              </Typography>
                              <Typography variant="body2" fontFamily="monospace">
                                {fitParams.ec50 != null
                                  ? fitParams.ec50.toExponential(2)
                                  : '-'}
                              </Typography>
                            </Box>
                            <Box>
                              <Typography variant="caption" color="text.secondary">
                                Hill
                              </Typography>
                              <Typography variant="body2" fontFamily="monospace">
                                {fitParams.hill != null ? fitParams.hill.toFixed(2) : '-'}
                              </Typography>
                            </Box>
                            <Box>
                              <Typography variant="caption" color="text.secondary">
                                Min
                              </Typography>
                              <Typography variant="body2" fontFamily="monospace">
                                {fitParams.minVal != null ? fitParams.minVal.toFixed(1) : '-'}
                              </Typography>
                            </Box>
                            <Box>
                              <Typography variant="caption" color="text.secondary">
                                Max
                              </Typography>
                              <Typography variant="body2" fontFamily="monospace">
                                {fitParams.maxVal != null ? fitParams.maxVal.toFixed(1) : '-'}
                              </Typography>
                            </Box>
                          </Box>
                        </Paper>
                      )}
                    </Box>
                  );
                })()
              ) : (
                <Paper
                  variant="outlined"
                  sx={{
                    height: 350,
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                  }}
                >
                  <Typography color="text.secondary">
                    Select a data series to view its chart
                  </Typography>
                </Paper>
              )}
            </Box>
          </Box>
        )}
      </DialogContent>

      <DialogActions>
        <Button
          component={Link}
          href={`/registry/compounds/${compoundId}`}
          startIcon={<Medication />}
          size="small"
        >
          View Compound
        </Button>
        <Box sx={{ flex: 1 }} />
        <Button onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
}
