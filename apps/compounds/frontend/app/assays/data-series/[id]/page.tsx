'use client';

import { use, useState } from 'react';
import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  Grid2 as Grid,
  Chip,
  Skeleton,
  Divider,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Button,
  CircularProgress,
  Tooltip,
} from '@mui/material';
import {
  Assessment,
  Medication,
  CheckCircle,
  Cancel,
  HelpOutline,
  ShowChart,
  Refresh,
  Science,
} from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DoseResponseChart } from '@/components/compounds/DoseResponseChart';
import { MoleculeView } from '@/components/compounds/MoleculeView';
import { AuthenticatedImage } from '@/components/compounds/AuthenticatedImage';
import { DilutionSeriesSelectDialog } from '@/components/compounds/DilutionSeriesSelectDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import { formatKpiUnit } from '@/lib/compounds/aggregation-api';
import { DataSeries, Assay, Compound, Protocol } from '@/types/compounds/models';

interface PageProps {
  params: Promise<{ id: string }>;
}

function InfoRow({ label, value }: { label: string; value: React.ReactNode }) {
  return (
    <Box sx={{ display: 'flex', py: 0.5 }}>
      <Typography
        color="text.secondary"
        sx={{ minWidth: 120, fontWeight: 500 }}
      >
        {label}:
      </Typography>
      <Typography component="div">{value ?? '-'}</Typography>
    </Box>
  );
}

function StatusChip({ status }: { status: string | undefined }) {
  switch (status) {
    case 'valid':
      return (
        <Chip
          icon={<CheckCircle />}
          label="Valid"
          size="small"
          color="success"
        />
      );
    case 'invalid':
      return (
        <Chip icon={<Cancel />} label="Invalid" size="small" color="error" />
      );
    default:
      return <Chip icon={<HelpOutline />} label="Unassigned" size="small" />;
  }
}

export default function DataSeriesDetailPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();
  const { canContribute } = useAuth();
  const [analysing, setAnalysing] = useState(false);
  const [dilutionDialogOpen, setDilutionDialogOpen] = useState(false);

  const { data: series, isLoading: seriesLoading, mutate: mutateSeries } = api.get<DataSeries>(
    `data-series/${id}/`
  );
  const { data: assay } = api.get<Assay>(
    series?.assay ? `assays/${series.assay}/` : null
  );
  const { data: protocol } = api.get<Protocol>(
    assay?.protocol ? `protocols/${assay.protocol}/` : null
  );
  const { data: compound } = api.get<Compound>(
    series?.compound ? `compounds/${series.compound}/` : null
  );

  const handleReanalyse = async () => {
    if (!id) return;
    setAnalysing(true);
    try {
      await api.post(`data-series/${id}/analyse/`, {});
      // Refresh the series data to show updated analysis
      mutateSeries();
    } catch (err) {
      console.error('Analysis error:', err);
    } finally {
      setAnalysing(false);
    }
  };

  // Extract chart data from series - requires dilution_series
  const chartData = series?.dilution_series?.concentrations && series?.extracted_data ? (() => {
    const concentrations = series.dilution_series.concentrations;
    let responses = Array.isArray(series.extracted_data) ? series.extracted_data : [];
    const unit = series.dilution_series.unit || 'nM';

    // Detect format: if extracted_data has 2 more elements than concentrations,
    // it has embedded controls at first and last positions
    if (responses.length === concentrations.length + 2) {
      responses = responses.slice(1, -1);
    }

    return { concentrations, responses, unit };
  })() : null;

  // Check for missing dilution series
  const missingDilutionSeries = series && !series.dilution_series;

  // Check if this is a table_of_values assay with an uploaded plot image
  const isTableOfValues = protocol?.analysis_method === 'table_of_values';
  const hasImageFile = series?.analysis?.results?.['Image File'];
  // Use the media proxy endpoint instead of direct blob URL (requires auth/SAS)
  const plotImageUrl = series?.plot_image
    ? `/api/proxy/compounds/media/data-series/${series.id}/plot/`
    : null;

  const fitParams = series?.analysis?.results ? (() => {
    const results = series.analysis.results;
    // Use KPI field to find the correct ec50 value
    const kpiKey = results.KPI || results.kpi;
    const ec50Value = kpiKey
      ? (results[kpiKey] ?? results[kpiKey.toLowerCase?.()] ?? results[kpiKey.toUpperCase?.()])
      : null;
    return {
      ec50: ec50Value ?? results.EC50 ?? results.ec50 ?? results.IC50 ?? results.ic50 ?? null,
      hill: results.Hill ?? results.hill ?? results.hill_slope ?? null,
      minVal: results.minVal ?? results.bottom ?? null,
      maxVal: results.maxVal ?? results.top ?? null,
      status: series.analysis.status,
    };
  })() : undefined;

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Assays', href: routes.assays.list(), icon: 'assay' },
          ...(assay?.protocol_name ? [{
            label: assay.protocol_name,
            href: routes.assays.protocol(assay.protocol),
            icon: 'protocol' as const,
          }] : []),
          {
            label: assay?.data_filename || 'Assay',
            href: series?.assay ? routes.assays.detail(series.assay) : undefined,
            icon: 'assay' as const,
          },
          { label: series?.compound_name || 'Data Series' },
        ]}
      />

      {/* Series header */}
      <Paper sx={{ p: 3, mb: 3 }}>
        {seriesLoading ? (
          <>
            <Skeleton variant="text" width={300} height={40} />
            <Skeleton variant="rectangular" height={300} />
          </>
        ) : series ? (
          <>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 2 }}>
              <ShowChart sx={{ fontSize: 48, color: 'secondary.main' }} />
              <Box sx={{ flex: 1 }}>
                <Typography variant="h4">
                  {series.compound_name || 'Unknown Compound'}
                </Typography>
                <Box sx={{ display: 'flex', gap: 1, mt: 0.5, flexWrap: 'wrap' }}>
                  <StatusChip status={series.analysis?.status} />
                  {series.compound_formatted_id && (
                    <Chip
                      icon={<Medication fontSize="small" />}
                      label={series.compound_formatted_id}
                      size="small"
                      onClick={() =>
                        series.compound && router.push(routes.registry.compound(series.compound))
                      }
                    />
                  )}
                  {assay && (
                    <Chip
                      icon={<Assessment fontSize="small" />}
                      label={assay.data_filename}
                      size="small"
                      variant="outlined"
                      onClick={() => router.push(routes.assays.detail(assay.id))}
                    />
                  )}
                </Box>
              </Box>
            </Box>

            <Divider sx={{ my: 2 }} />

            <Grid container spacing={3}>
              {/* Chart */}
              <Grid size={{ xs: 12, md: 7 }}>
                <Typography variant="h6" gutterBottom>
                  {isTableOfValues ? 'Plot Image' : 'Dose-Response Curve'}
                </Typography>
                {/* For table_of_values, always use plot_image - never attempt interactive chart */}
                {isTableOfValues ? (
                  hasImageFile && plotImageUrl ? (
                    <AuthenticatedImage
                      src={plotImageUrl}
                      alt={`Plot for ${series.compound_name || 'compound'}`}
                      width="100%"
                      height="auto"
                      objectFit="contain"
                      sx={{
                        maxWidth: '100%',
                        maxHeight: 400,
                        borderRadius: 1,
                        border: '1px solid',
                        borderColor: 'divider',
                      }}
                    />
                  ) : hasImageFile ? (
                    <Paper
                      sx={{
                        p: 4,
                        bgcolor: 'warning.50',
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                        height: 350,
                        flexDirection: 'column',
                        gap: 1,
                      }}
                    >
                      <Typography color="warning.main" fontWeight={500}>
                        Image not uploaded
                      </Typography>
                      <Typography variant="caption" color="text.secondary">
                        Expected: {hasImageFile}
                      </Typography>
                    </Paper>
                  ) : (
                    <Paper
                      sx={{
                        p: 4,
                        bgcolor: 'grey.50',
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                        height: 350,
                        flexDirection: 'column',
                        gap: 1,
                      }}
                    >
                      <Typography color="text.secondary">
                        No plot image available
                      </Typography>
                    </Paper>
                  )
                ) : chartData && chartData.concentrations.length > 0 ? (
                  <DoseResponseChart
                    data={chartData}
                    fit={fitParams}
                    compoundName={series.compound_name || undefined}
                    width={500}
                    height={350}
                    skipPoints={Array.isArray(series.skip_points) ? series.skip_points : []}
                  />
                ) : (
                  <Paper
                    sx={{
                      p: 4,
                      bgcolor: missingDilutionSeries ? 'error.50' : 'grey.50',
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                      height: 350,
                      flexDirection: 'column',
                      gap: 1,
                    }}
                  >
                    <Typography color={missingDilutionSeries ? 'error' : 'text.secondary'} fontWeight={missingDilutionSeries ? 500 : 400}>
                      {missingDilutionSeries
                        ? 'No dilution series configured'
                        : 'No concentration data available'}
                    </Typography>
                    {missingDilutionSeries && (
                      <Typography variant="caption" color="text.secondary">
                        Set preferred dilutions on the protocol to enable analysis
                      </Typography>
                    )}
                  </Paper>
                )}
              </Grid>

              {/* Analysis Results */}
              <Grid size={{ xs: 12, md: 5 }}>
                <Typography variant="h6" gutterBottom>
                  Analysis Results
                </Typography>
                {series.analysis ? (
                  <TableContainer component={Paper} variant="outlined">
                    <Table size="small">
                      <TableBody>
                        <TableRow>
                          <TableCell component="th">Status</TableCell>
                          <TableCell>
                            <StatusChip status={series.analysis.status} />
                          </TableCell>
                        </TableRow>
                        {series.analysis.results?.error && (
                          <TableRow>
                            <TableCell component="th">Error</TableCell>
                            <TableCell>
                              <Typography color="error" variant="body2">
                                {series.analysis.results.error}
                              </Typography>
                            </TableCell>
                          </TableRow>
                        )}
                        {series.analysis.results?.flags && series.analysis.results.flags.length > 0 && (
                          <TableRow>
                            <TableCell component="th">Flags</TableCell>
                            <TableCell>
                              {series.analysis.results.flags.map((flag: string) => (
                                <Chip
                                  key={flag}
                                  label={flag}
                                  size="small"
                                  color="warning"
                                  variant="outlined"
                                  sx={{ mr: 0.5, mb: 0.5 }}
                                />
                              ))}
                            </TableCell>
                          </TableRow>
                        )}
                        {/* KPI row - highlighted at top */}
                        {series.analysis.kpi_value !== undefined && (() => {
                          // Get unit from analysis results first, then fall back to dilution_series
                          const kpiUnit = series.analysis.results?.kpi_unit
                            || series.dilution_series?.unit
                            || null;
                          return (
                            <TableRow sx={{ bgcolor: 'primary.50' }}>
                              <TableCell component="th">
                                <Typography fontWeight={600}>
                                  KPI ({series.analysis.results?.KPI || 'EC50'})
                                </Typography>
                              </TableCell>
                              <TableCell>
                                <Typography fontFamily="monospace" fontWeight={600} color="primary">
                                  {typeof series.analysis.kpi_value === 'number'
                                    ? (series.analysis.kpi_value >= 0.01 && series.analysis.kpi_value < 10000
                                        ? series.analysis.kpi_value.toFixed(2)
                                        : series.analysis.kpi_value.toExponential(3))
                                    : series.analysis.kpi_value}
                                  {kpiUnit && (
                                    <Typography
                                      component="span"
                                      color="text.secondary"
                                      sx={{ ml: 1, fontWeight: 400 }}
                                    >
                                      {formatKpiUnit(kpiUnit)}
                                    </Typography>
                                  )}
                                </Typography>
                              </TableCell>
                            </TableRow>
                          );
                        })()}
                        {/* Dynamically render all other results fields */}
                        {series.analysis.results && Object.entries(series.analysis.results)
                          .filter(([key]) => !['error', 'flags', 'KPI', 'source', 'time_course', 'curve_points'].includes(key.toLowerCase()))
                          .map(([key, value]) => {
                            // Skip arrays of arrays (e.g., coordinate data better shown in charts)
                            if (Array.isArray(value) && value.length > 0 && Array.isArray(value[0])) {
                              return null;
                            }
                            // Format the key for display (snake_case to Title Case)
                            const label = key
                              .replace(/_/g, ' ')
                              .replace(/\b\w/g, c => c.toUpperCase())
                              .replace(/T1 2/i, 't½');

                            // Format the value based on type
                            let displayValue: React.ReactNode;
                            if (value === null || value === undefined) {
                              displayValue = '-';
                            } else if (typeof value === 'number') {
                              displayValue = (
                                <Typography fontFamily="monospace">
                                  {value >= 0.01 && value < 10000
                                    ? value.toFixed(2)
                                    : value.toExponential(3)}
                                </Typography>
                              );
                            } else if (typeof value === 'object') {
                              // For nested objects, show as formatted JSON
                              displayValue = (
                                <Typography
                                  fontFamily="monospace"
                                  variant="body2"
                                  sx={{ whiteSpace: 'pre-wrap', maxWidth: 200 }}
                                >
                                  {JSON.stringify(value, null, 2)}
                                </Typography>
                              );
                            } else {
                              displayValue = String(value);
                            }

                            return (
                              <TableRow key={key}>
                                <TableCell component="th">{label}</TableCell>
                                <TableCell>{displayValue}</TableCell>
                              </TableRow>
                            );
                          })}
                      </TableBody>
                    </Table>
                  </TableContainer>
                ) : (
                  <Paper sx={{ p: 2, bgcolor: 'grey.50' }}>
                    <Typography color="text.secondary">
                      No analysis results available
                    </Typography>
                  </Paper>
                )}

                {/* Action buttons */}
                <Box sx={{ display: 'flex', gap: 1, mt: 2, flexWrap: 'wrap' }}>
                  {/* Change Dilution Series button */}
                  <Tooltip
                    title={
                      !canContribute
                        ? 'Requires Contributor or Admin operating level'
                        : 'Select a different dilution series for this data'
                    }
                    arrow
                  >
                    <span>
                      <Button
                        variant="outlined"
                        startIcon={<Science />}
                        onClick={() => setDilutionDialogOpen(true)}
                        disabled={!canContribute}
                        size="small"
                      >
                        Change Dilutions
                      </Button>
                    </span>
                  </Tooltip>

                  {/* Re-analyse button */}
                  <Tooltip
                    title={
                      !canContribute
                        ? 'Requires Contributor or Admin operating level'
                        : missingDilutionSeries
                        ? 'Cannot analyse: No dilution series configured'
                        : 'Re-run curve fitting analysis'
                    }
                    arrow
                  >
                    <span>
                      <Button
                        variant="outlined"
                        startIcon={analysing ? <CircularProgress size={16} /> : <Refresh />}
                        onClick={handleReanalyse}
                        disabled={analysing || missingDilutionSeries || !canContribute}
                        size="small"
                      >
                        {analysing ? 'Analysing...' : 'Re-analyse'}
                      </Button>
                    </span>
                  </Tooltip>
                </Box>

                {/* Compound structure if matched */}
                {compound && (
                  <Box sx={{ mt: 3 }}>
                    <Typography variant="h6" gutterBottom>
                      Compound Structure
                    </Typography>
                    <Box sx={{ display: 'flex', gap: 2, alignItems: 'flex-start' }}>
                      <MoleculeView smiles={compound.smiles} width={150} height={150} />
                      <Box>
                        <InfoRow label="ID" value={compound.formatted_id} />
                        <InfoRow
                          label="MW"
                          value={compound.molecular_weight?.toFixed(1)}
                        />
                        <InfoRow label="Target" value={compound.target_name} />
                      </Box>
                    </Box>
                  </Box>
                )}
              </Grid>
            </Grid>

            {/* Raw data table */}
            {chartData && chartData.concentrations.length > 0 && (() => {
              // Extract control values from extracted_data if present
              const rawData = Array.isArray(series.extracted_data) ? series.extracted_data : [];
              const hasControls = rawData.length === chartData.concentrations.length + 2;
              const minControl = hasControls ? rawData[0] : null;
              const maxControl = hasControls ? rawData[rawData.length - 1] : null;

              return (
                <>
                  <Divider sx={{ my: 3 }} />
                  <Typography variant="h6" gutterBottom>
                    Raw Data
                  </Typography>
                  <TableContainer component={Paper} variant="outlined">
                    <Table size="small">
                      <TableHead>
                        <TableRow>
                          <TableCell>Point</TableCell>
                          <TableCell>
                            Concentration ({series.dilution_series?.unit || 'nM'})
                          </TableCell>
                          <TableCell>Response</TableCell>
                          <TableCell>Included</TableCell>
                        </TableRow>
                      </TableHead>
                      <TableBody>
                        {/* Min control row */}
                        {hasControls && (
                          <TableRow sx={{ bgcolor: 'info.50' }}>
                            <TableCell>
                              <Chip label="Min Ctrl" size="small" color="info" variant="outlined" />
                            </TableCell>
                            <TableCell>
                              <Typography variant="body2" color="text.secondary">
                                (low signal control)
                              </Typography>
                            </TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace" fontWeight={500}>
                                {typeof minControl === 'number' ? minControl.toFixed(1) : '-'}
                              </Typography>
                            </TableCell>
                            <TableCell>
                              <Typography variant="caption" color="text.secondary">
                                Control
                              </Typography>
                            </TableCell>
                          </TableRow>
                        )}
                        {/* Data points */}
                        {chartData.concentrations.map((conc, idx) => (
                          <TableRow
                            key={idx}
                            sx={{
                              bgcolor: (Array.isArray(series.skip_points) ? series.skip_points : []).includes(idx)
                                ? 'grey.100'
                                : undefined,
                            }}
                          >
                            <TableCell>{idx + 1}</TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace">
                                {conc.toExponential(2)}
                              </Typography>
                            </TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace">
                                {chartData.responses[idx]?.toFixed(1) ?? '-'}
                              </Typography>
                            </TableCell>
                            <TableCell>
                              {(Array.isArray(series.skip_points) ? series.skip_points : []).includes(idx) ? (
                                <Cancel fontSize="small" color="disabled" />
                              ) : (
                                <CheckCircle fontSize="small" color="success" />
                              )}
                            </TableCell>
                          </TableRow>
                        ))}
                        {/* Max control row */}
                        {hasControls && (
                          <TableRow sx={{ bgcolor: 'primary.50' }}>
                            <TableCell>
                              <Chip label="Max Ctrl" size="small" color="primary" variant="outlined" />
                            </TableCell>
                            <TableCell>
                              <Typography variant="body2" color="text.secondary">
                                (high signal control)
                              </Typography>
                            </TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace" fontWeight={500}>
                                {typeof maxControl === 'number' ? maxControl.toFixed(1) : '-'}
                              </Typography>
                            </TableCell>
                            <TableCell>
                              <Typography variant="caption" color="text.secondary">
                                Control
                              </Typography>
                            </TableCell>
                          </TableRow>
                        )}
                      </TableBody>
                    </Table>
                  </TableContainer>
                </>
              );
            })()}

            {/* Time Course Data (for ADME assays) */}
            {series.analysis?.results?.time_course && (() => {
              const timeCourse = series.analysis.results.time_course;
              const timePoints = timeCourse.time_points_min || [];
              const remainingPct = timeCourse.remaining_pct || [];

              if (timePoints.length === 0) return null;

              return (
                <>
                  <Divider sx={{ my: 3 }} />
                  <Typography variant="h6" gutterBottom>
                    Time Course Data
                  </Typography>
                  <TableContainer component={Paper} variant="outlined">
                    <Table size="small">
                      <TableHead>
                        <TableRow>
                          <TableCell>Time (min)</TableCell>
                          <TableCell>Remaining (%)</TableCell>
                        </TableRow>
                      </TableHead>
                      <TableBody>
                        {timePoints.map((time: number, idx: number) => (
                          <TableRow key={idx}>
                            <TableCell>
                              <Typography fontFamily="monospace">{time}</Typography>
                            </TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace">
                                {remainingPct[idx] !== null && remainingPct[idx] !== undefined
                                  ? typeof remainingPct[idx] === 'number'
                                    ? remainingPct[idx].toFixed(1)
                                    : remainingPct[idx]
                                  : '-'}
                              </Typography>
                            </TableCell>
                          </TableRow>
                        ))}
                      </TableBody>
                    </Table>
                  </TableContainer>
                </>
              );
            })()}

            {/* Source info */}
            <Divider sx={{ my: 3 }} />
            <Typography variant="h6" gutterBottom>
              Source Information
            </Typography>
            <Grid container spacing={2}>
              <Grid size={{ xs: 12, sm: 6 }}>
                <InfoRow label="Row" value={series.row} />
                <InfoRow
                  label="Columns"
                  value={`${series.start_column} - ${series.end_column}`}
                />
              </Grid>
              <Grid size={{ xs: 12, sm: 6 }}>
                <InfoRow label="Assay" value={assay?.data_filename} />
                <InfoRow label="Protocol" value={assay?.protocol_name} />
              </Grid>
            </Grid>
          </>
        ) : (
          <Typography color="error">Data series not found</Typography>
        )}
      </Paper>

      {/* Dilution Series Selection Dialog */}
      {series && (
        <DilutionSeriesSelectDialog
          open={dilutionDialogOpen}
          onClose={() => setDilutionDialogOpen(false)}
          dataSeriesId={series.id}
          currentDilutionSeriesId={series.dilution_series?.id}
          protocolDilutionSeriesId={protocol?.preferred_dilutions || undefined}
          onSave={() => mutateSeries()}
        />
      )}
    </Container>
  );
}
