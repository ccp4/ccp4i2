'use client';

import { use } from 'react';
import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  Grid,
  Chip,
  Skeleton,
  Divider,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
} from '@mui/material';
import {
  Assessment,
  Science,
  Medication,
  CheckCircle,
  Cancel,
  HelpOutline,
  ShowChart,
} from '@mui/icons-material';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { DoseResponseChart } from '@/components/DoseResponseChart';
import { MoleculeView } from '@/components/MoleculeView';
import { useCompoundsApi } from '@/lib/api';
import { DataSeries, Assay, Compound } from '@/types/models';

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

  const { data: series, isLoading: seriesLoading } = api.get<DataSeries>(
    `data-series/${id}/`
  );
  const { data: assay } = api.get<Assay>(
    series?.assay ? `assays/${series.assay}/` : null
  );
  const { data: compound } = api.get<Compound>(
    series?.compound ? `compounds/${series.compound}/` : null
  );

  // Extract chart data from series
  const chartData = series?.dilution_series && series?.extracted_data ? {
    concentrations: series.dilution_series.concentrations || [],
    responses: Array.isArray(series.extracted_data)
      ? series.extracted_data
      : [],
    unit: series.dilution_series.unit,
  } : null;

  const fitParams = series?.analysis ? {
    ec50: series.analysis.results?.EC50,
    hill: series.analysis.results?.Hill,
    minVal: series.analysis.results?.minVal,
    maxVal: series.analysis.results?.maxVal,
    status: series.analysis.status,
  } : undefined;

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Breadcrumbs
        items={[
          { label: 'Home', href: '/', icon: 'home' },
          { label: 'Assays', href: '/assays', icon: 'assay' },
          {
            label: assay?.data_filename || 'Assay',
            href: series?.assay ? `/assays/${series.assay}` : undefined,
            icon: 'assay',
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
                        series.compound && router.push(`/registry/compounds/${series.compound}`)
                      }
                    />
                  )}
                  {assay && (
                    <Chip
                      icon={<Assessment fontSize="small" />}
                      label={assay.data_filename}
                      size="small"
                      variant="outlined"
                      onClick={() => router.push(`/assays/${assay.id}`)}
                    />
                  )}
                </Box>
              </Box>
            </Box>

            <Divider sx={{ my: 2 }} />

            <Grid container spacing={3}>
              {/* Chart */}
              <Grid item xs={12} md={7}>
                <Typography variant="h6" gutterBottom>
                  Dose-Response Curve
                </Typography>
                {chartData && chartData.concentrations.length > 0 ? (
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
                      bgcolor: 'grey.50',
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                      height: 350,
                    }}
                  >
                    <Typography color="text.secondary">
                      No concentration data available
                    </Typography>
                  </Paper>
                )}
              </Grid>

              {/* Analysis Results */}
              <Grid item xs={12} md={5}>
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
                        {series.analysis.results?.EC50 !== undefined && (
                          <TableRow>
                            <TableCell component="th">EC50</TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace">
                                {typeof series.analysis.results.EC50 === 'number'
                                  ? series.analysis.results.EC50.toExponential(3)
                                  : series.analysis.results.EC50}{' '}
                                {series.dilution_series?.unit || 'nM'}
                              </Typography>
                            </TableCell>
                          </TableRow>
                        )}
                        {series.analysis.results?.Hill !== undefined && (
                          <TableRow>
                            <TableCell component="th">Hill Coefficient</TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace">
                                {typeof series.analysis.results.Hill === 'number'
                                  ? series.analysis.results.Hill.toFixed(2)
                                  : series.analysis.results.Hill}
                              </Typography>
                            </TableCell>
                          </TableRow>
                        )}
                        {series.analysis.results?.minVal !== undefined && (
                          <TableRow>
                            <TableCell component="th">Min Value</TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace">
                                {typeof series.analysis.results.minVal === 'number'
                                  ? series.analysis.results.minVal.toFixed(1)
                                  : series.analysis.results.minVal}
                              </Typography>
                            </TableCell>
                          </TableRow>
                        )}
                        {series.analysis.results?.maxVal !== undefined && (
                          <TableRow>
                            <TableCell component="th">Max Value</TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace">
                                {typeof series.analysis.results.maxVal === 'number'
                                  ? series.analysis.results.maxVal.toFixed(1)
                                  : series.analysis.results.maxVal}
                              </Typography>
                            </TableCell>
                          </TableRow>
                        )}
                        {series.analysis.kpi_value !== undefined && (
                          <TableRow>
                            <TableCell component="th">
                              KPI ({series.analysis.results?.KPI || 'EC50'})
                            </TableCell>
                            <TableCell>
                              <Typography fontFamily="monospace" fontWeight={600}>
                                {typeof series.analysis.kpi_value === 'number'
                                  ? series.analysis.kpi_value.toExponential(3)
                                  : series.analysis.kpi_value}
                              </Typography>
                            </TableCell>
                          </TableRow>
                        )}
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
            {chartData && chartData.concentrations.length > 0 && (
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
                    </TableBody>
                  </Table>
                </TableContainer>
              </>
            )}

            {/* Source info */}
            <Divider sx={{ my: 3 }} />
            <Typography variant="h6" gutterBottom>
              Source Information
            </Typography>
            <Grid container spacing={2}>
              <Grid item xs={12} sm={6}>
                <InfoRow label="Row" value={series.row} />
                <InfoRow
                  label="Columns"
                  value={`${series.start_column} - ${series.end_column}`}
                />
              </Grid>
              <Grid item xs={12} sm={6}>
                <InfoRow label="Assay" value={assay?.data_filename} />
                <InfoRow label="Protocol" value={assay?.protocol_name} />
              </Grid>
            </Grid>
          </>
        ) : (
          <Typography color="error">Data series not found</Typography>
        )}
      </Paper>
    </Container>
  );
}
